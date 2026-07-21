# ORFquant mirai parallel backend — Disk-Load strategy
#
# Replaces parallel::mclapply (fork) with mirai daemons (socket-based
# independent processes, NNG nanonext transport).  Eliminates GC finalizer
# crashes on FaFile/Bioconductor reference objects that plague forked
# subprocesses.
#
# Architecture (v2 — Disk-Load):
#   Main process saves genome_seq as RDS → starts mirai daemons
#   → everywhere() sets ONLY paths + scalar params (NOT large objects)
#   → Each daemon independently loads GTF_annotation from Rannot file,
#     for_ORFquant_data from for_ORFquant file, genome_seq from RDS
#   → mirai_map() dispatches gene regions to daemons
#   → Each daemon calls ORFquant() independently (no fork, no GC issues)
#
# Memory: 1 copy per daemon (~5GB each) — NOT serialized via everywhere().
#
# Requirements:
#   - mirai  (>= 0.13.0) — async daemon framework (CRAN)
#   - nanonext            — NNG transport (auto-installed with mirai)

orfquant_mirai_parallel <- function(
    genes_red,
    for_ORFquant_data,
    GTF_annotation,
    genome_seq,
    for_ORFquant_file,
    annotation_file,
    n_cores,
    canonical_start_only,
    stn.orf_find.all_starts,
    stn.orf_find.nostarts,
    stn.orf_find.start_sel_cutoff,
    stn.orf_find.start_sel_cutoff_ave,
    stn.orf_find.cutoff_fr_ave,
    stn.orf_quant.cutoff_cums,
    stn.orf_quant.cutoff_pct,
    stn.orf_quant.cutoff_P_sites,
    unique_reads_only
) {
    if (!requireNamespace("mirai", quietly = TRUE)) {
        stop(
            "mirai package is required for parallel_backend = 'mirai'.\n",
            "  Install with: install.packages('mirai', repos = 'https://cloud.r-project.org')"
        )
    }

    n_regions <- length(genes_red)
    cat(sprintf("[mirai] Preparing disk-load strategy for %d gene regions... %s\n",
        n_regions, date()))

    # ---- Step 1: Limit daemon count to avoid memory exhaustion ----
    # Each daemon loads ~5GB (human genome). Allow ~10GB overhead per daemon.
    mem_gb <- tryCatch(
        as.numeric(system("awk '/MemAvailable/{print $2/1024/1024}' /proc/meminfo",
                           intern = TRUE)),
        error = function(e) 32
    )
    max_daemons <- max(1L, as.integer(floor(mem_gb / 10)))
    n_cores <- min(n_cores, max_daemons, 16L)
    cat(sprintf("[mirai] Memory: %.0f GB available, capping at %d daemons\n",
        mem_gb, n_cores))

    # ---- Step 2: Start daemon pool ----
    cat(sprintf("[mirai] Starting %d daemons (mirai %s)... %s\n",
        n_cores, as.character(packageVersion("mirai")), date()))

    mirai::daemons(
        n          = n_cores,
        dispatcher = TRUE
    )
    on.exit(
        tryCatch(
            mirai::daemons(0),
            error = function(e) message("[mirai] cleanup: ", conditionMessage(e))
        ),
        add = TRUE
    )

    # ---- Step 3: Set up daemon environment ----
    # Packages and file paths via .expr; scalar params and genes_red as ...
    # (assigned to daemon global env, making them visible to mirai_map callbacks).
    cat(sprintf("[mirai] Broadcasting paths and packages to %d daemons... %s\n",
        n_cores, date()))

    # Pre-load annotation data in each daemon at startup.
    # <<- inside everywhere({}) pushes objects to the daemon's long-lived
    # environment, which persists across mirai_map callbacks (unlike .GlobalEnv
    # which mirai creates fresh for each callback invocation).
    mirai::everywhere({
        suppressPackageStartupMessages({
            library(GenomicRanges)
            library(GenomicFeatures)
            library(Biostrings)
            library(Rsamtools)
            library(ORFquant)
        })
        cat(sprintf("[daemon %d] Loading annotation + P-sites from disk...\n",
            Sys.getpid()))
        ORFquant::load_annotation(ANNOTATION_FILE)
        GTF_annotation <<- GTF_annotation
        genome_seq <<- genome_seq
        for_ORFquant_data <<- get(load(FOR_ORFQUANT_FILE))
        cat(sprintf("[daemon %d] Loaded: GTF=%.0fMB genome=%.0fMB pdata=%.0fMB\n",
            Sys.getpid(),
            as.numeric(object.size(GTF_annotation)) / 1e6,
            as.numeric(object.size(genome_seq)) / 1e6,
            as.numeric(object.size(for_ORFquant_data)) / 1e6))
    },
        ANNOTATION_FILE       = annotation_file,
        FOR_ORFQUANT_FILE     = for_ORFquant_file,
        genes_red             = genes_red,
        canonical_start_only  = canonical_start_only,
        unique_reads_only     = unique_reads_only,
        stn.orf_find.all_starts          = stn.orf_find.all_starts,
        stn.orf_find.nostarts            = stn.orf_find.nostarts,
        stn.orf_find.start_sel_cutoff    = stn.orf_find.start_sel_cutoff,
        stn.orf_find.start_sel_cutoff_ave = stn.orf_find.start_sel_cutoff_ave,
        stn.orf_find.cutoff_fr_ave       = stn.orf_find.cutoff_fr_ave,
        stn.orf_quant.cutoff_cums        = stn.orf_quant.cutoff_cums,
        stn.orf_quant.cutoff_pct         = stn.orf_quant.cutoff_pct,
        stn.orf_quant.cutoff_P_sites     = stn.orf_quant.cutoff_P_sites
    )

    # ---- Step 5: Parallel computation ----
    cat(sprintf("[mirai] Processing %d gene regions with %d daemons... %s\n",
        n_regions, n_cores, date()))

    results <- mirai::mirai_map(
        seq_along(genes_red),
        function(g) {
            tryCatch({

                gen_region <- genes_red[g]
                chr_name <- as.character(seqnames(gen_region))

                # Match process_gene logic: genetic_codes is a data.frame,
                # the genetic_code column has no names — use rownames() instead.
                code_id <- GTF_annotation$genetic_codes$genetic_code[
                    rownames(GTF_annotation$genetic_codes) == chr_name
                ]
                genetcd <- getGeneticCode(code_id)

                if (canonical_start_only) {
                    attributes(genetcd)$alt_init_codons <- names(
                        which(genetcd == "M")
                    )
                }

                ORFquant(
                    region = gen_region,
                    for_ORFquant = for_ORFquant_data,
                    genetic_code_region = genetcd,
                    orf_find.all_starts = stn.orf_find.all_starts,
                    orf_find.nostarts = stn.orf_find.nostarts,
                    orf_find.start_sel_cutoff = stn.orf_find.start_sel_cutoff,
                    orf_find.start_sel_cutoff_ave = stn.orf_find.start_sel_cutoff_ave,
                    orf_find.cutoff_fr_ave = stn.orf_find.cutoff_fr_ave,
                    orf_quant.cutoff_cums = stn.orf_quant.cutoff_cums,
                    orf_quant.cutoff_pct = stn.orf_quant.cutoff_pct,
                    orf_quant.cutoff_P_sites = stn.orf_quant.cutoff_P_sites,
                    unique_reads = unique_reads_only
                )
            }, error = function(e) {
                message(sprintf(
                    "\n[mirai] Gene region %d (%s) error: %s",
                    g, as.character(genes_red[g]), conditionMessage(e)
                ))
                NULL
            })
        }
    )[]

    # ---- Step 6: Filter failed regions (NULL, try-error, empty list) ----
    # Consistent with mclapply path filtering in orfquant.R
    is_invalid <- vapply(results, function(x) {
        inherits(x, "try-error") ||
            is.null(x) ||
            (is.list(x) && length(x) == 0)
    }, logical(1L))
    n_failed <- sum(is_invalid)
    if (n_failed > 0) {
        cat(sprintf(
            "\n[mirai] %d / %d gene regions failed or empty, %d succeeded\n",
            n_failed, n_regions, n_regions - n_failed
        ))
    }
    if (n_failed > 0) {
        results <- results[!is_invalid]
    }

    cat(sprintf("[mirai] Processing complete. %d regions successful. %s\n",
        length(results), date()))

    return(results)
}
