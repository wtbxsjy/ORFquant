# ORFquant mirai parallel backend — OPTIMIZED Disk-Load strategy (v3)
#
# KEY OPTIMIZATION over v2:
#   v2: each daemon calls load_annotation() → FaFile → DNAStringSet (~2.2 GB in RAM)
#   v3: daemon keeps genome as FaFile (file descriptor, ~0 MB RAM), reads sequences on demand
#        Total per-daemon memory: ~300 MB (annotation only) vs ~2.5 GB (annotation + genome)
#        Allows 16+ daemons on the same hardware without OOM.
#
# Memory model:
#   Main process: saves genome path + annotation file path → starts mirai daemons
#   everywhere(): loads GTF_annotation from Rannot, opens FaFile from path
#   mirai_map(): dispatches gene indices only (integers)
#   Per daemon: ~300 MB (GTF_annotation GRangesList) + ~0 MB (FaFile) = ~300 MB
#   16 daemons: ~5 GB total (vs ~40 GB in v2)
#
# Requirements:
#   - mirai (>= 0.13.0)
#   - nanonext
#   - ORFquant internal: .orfquant_genome_ref(), .orfquant_open_genome_ref()

orfquant_mirai_parallel_v3 <- function(
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
    cat(sprintf("[mirai v3] Preparing disk-load strategy for %d gene regions... %s\n",
        n_regions, date()))

    # ---- Step 1: Extract genome reference (file path, NOT DNAStringSet) ----
    # .orfquant_genome_ref() stores the FASTA file path without loading sequences.
    # Each daemon will open its own FaFile handle — near-zero memory overhead.
    genome_ref <- .orfquant_genome_ref(genome_seq)
    if (is.null(genome_ref)) {
        # genome_seq is not a FaFile (e.g. BSgenome, already DNAStringSet)
        # Fall back to saving as temporary RDS for daemons to load.
        genome_rds <- tempfile(fileext = ".rds")
        saveRDS(genome_seq, genome_rds, compress = FALSE)
        use_genome_ref <- FALSE
        genome_path_or_rds <- genome_rds
        genome_ref_obj <- NULL
        cat(sprintf("[mirai v3] Genome not FaFile — saved as RDS (%s)\n",
            format(file.size(genome_rds), units = "MB")))
    } else {
        use_genome_ref <- TRUE
        # Normalize to absolute path — daemons may run from a different working dir
        genome_path_or_rds <- normalizePath(genome_ref$path, mustWork = TRUE)
        genome_ref_obj <- genome_ref
        # Also update the path in genome_ref_obj so .orfquant_open_genome_ref opens correctly
        genome_ref_obj$path <- genome_path_or_rds
        cat(sprintf("[mirai v3] Genome reference: %s (daemons open FaFile, ~0MB each)\n",
            genome_ref$path))
    }

    # ---- Step 2: Limit daemon count based on per-daemon memory ----
    # Measured ~2.9 GB per daemon for genome_ref (GOTCHA #32); configurable via ORFQ_PER_DAEMON_GB
    mem_gb <- tryCatch(
        as.numeric(system("awk '/MemAvailable/{print $2/1024/1024}' /proc/meminfo",
                           intern = TRUE)),
        error = function(e) 32
    )
    env_per_daemon <- as.numeric(Sys.getenv("ORFQ_PER_DAEMON_GB", ""))
    per_daemon_gb <- if (!is.na(env_per_daemon) && env_per_daemon > 0) {
        env_per_daemon
    } else {
        if (use_genome_ref) 2.9 else 4.0
    }
    max_daemons <- max(1L, as.integer(floor(mem_gb / per_daemon_gb)))
    req_cores <- n_cores
    n_cores <- min(n_cores, max_daemons, 32L)
    if (req_cores > n_cores) {
        cat(sprintf(
            "[mirai v3] Available memory (%.1f GB) throttles daemons from %d to %d (%.1f GB/daemon)\n",
            mem_gb, req_cores, n_cores, per_daemon_gb
        ))
    } else {
        cat(sprintf("[mirai v3] Memory: %.0f GB available, capping at %d daemons (%.1f GB/daemon)\n",
            mem_gb, n_cores, per_daemon_gb))
    }

    # ---- Step 3: Start daemon pool ----
    cat(sprintf("[mirai v3] Starting %d daemons (mirai %s)... %s\n",
        n_cores, as.character(packageVersion("mirai")), date()))

    mirai::daemons(
        n          = n_cores,
        dispatcher = TRUE
    )
    on.exit(
        tryCatch({
            mirai::daemons(0)
            # Cleanup temp RDS if used
            if (!use_genome_ref && file.exists(genome_path_or_rds)) {
                unlink(genome_path_or_rds)
            }
        }, error = function(e) message("[mirai v3] cleanup: ", conditionMessage(e))),
        add = TRUE
    )

    # ---- Step 4: Set up daemon environment ----
    cat(sprintf("[mirai v3] Broadcasting paths and packages to %d daemons... %s\n",
        n_cores, date()))

    # Capture the genome loading strategy for the daemon closure
    USE_GENOME_REF <- use_genome_ref
    GENOME_PATH_OR_RDS <- genome_path_or_rds
    GENOME_REF_OBJ <- genome_ref_obj

    mirai::everywhere({
        suppressPackageStartupMessages({
            library(GenomicRanges)
            library(GenomicFeatures)
            library(Biostrings)
            library(Rsamtools)
            library(ORFquant)
        })
        cat(sprintf("[daemon %d] Loading data from disk...\n", Sys.getpid()))

        # ── Load annotation WITHOUT FaFile→DNAStringSet conversion ──
        # load_annotation() converts FaFile → DNAStringSet (2.2 GB) for fork
        # safety. mirai daemons are independent socket processes — they don't
        # fork, so we load the Rannot directly and keep the genome as FaFile.
        load_env <- new.env(parent = emptyenv())
        ann_raw <- get(load(ANNOTATION_FILE, envir = load_env), envir = load_env)

        # Build annotation without the genome conversion
        GTF_annotation <<- ann_raw
        # Also assign to .GlobalEnv so ORFquant() internal functions find it
        assign("GTF_annotation", ann_raw, envir = .GlobalEnv)

        # ── Load genome — disk-backed FaFile (~0 MB RAM) ──
        if (USE_GENOME_REF) {
            genome_seq <<- ORFquant:::.orfquant_open_genome_ref(GENOME_REF_OBJ)
            assign("genome_seq", genome_seq, envir = .GlobalEnv)
            cat(sprintf("[daemon %d] Genome: FaFile (on-disk path, ~0MB RAM)\n",
                Sys.getpid()))
        } else {
            genome_ds <- readRDS(GENOME_PATH_OR_RDS)
            genome_seq <<- genome_ds
            assign("genome_seq", genome_ds, envir = .GlobalEnv)
            cat(sprintf("[daemon %d] Genome: DNAStringSet from RDS (%.0fMB RAM)\n",
                Sys.getpid(), as.numeric(object.size(genome_ds)) / 1e6))
        }

        # ── Load P-sites data ──
        for_ORFquant_data <<- get(load(FOR_ORFQUANT_FILE))
        cat(sprintf("[daemon %d] P-sites loaded: class=%s fields=%s psites_len=%d\n",
            Sys.getpid(), class(for_ORFquant_data)[1],
            paste(names(for_ORFquant_data), collapse=","),
            length(for_ORFquant_data$P_sites_all)))

        cat(sprintf("[daemon %d] Ready: annot=%.0fMB genome=%.0fMB psites=%.0fMB (total %.0fMB)\n",
            Sys.getpid(),
            as.numeric(object.size(GTF_annotation)) / 1e6,
            as.numeric(object.size(genome_seq)) / 1e6,
            as.numeric(object.size(for_ORFquant_data)) / 1e6,
            (as.numeric(object.size(GTF_annotation)) +
             as.numeric(object.size(genome_seq)) +
             as.numeric(object.size(for_ORFquant_data))) / 1e6
        ))
    },
        ANNOTATION_FILE       = annotation_file,
        FOR_ORFQUANT_FILE     = for_ORFquant_file,
        USE_GENOME_REF        = USE_GENOME_REF,
        GENOME_REF_OBJ        = GENOME_REF_OBJ,
        GENOME_PATH_OR_RDS    = GENOME_PATH_OR_RDS,
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

    # ---- Step 5: Parallel computation with DISK-BACKED streaming ----
    # Each daemon writes its result to an RDS file immediately.
    # The main process only collects file paths (~100 bytes each), NOT GRanges objects.
    # This avoids the ~150 GB memory explosion from collecting all ORF results in RAM.
    cat(sprintf("[mirai v3] Processing %d gene regions with %d daemons (disk-backed)... %s\n",
        n_regions, n_cores, date()))

    result_dir <- tempfile("orfq_mirai_results_")
    dir.create(result_dir)
    cat(sprintf("[mirai v3] Results directory: %s\n", result_dir))

    result_paths <- mirai::mirai_map(
        seq_along(genes_red),
        function(g) {
            # Resolve from daemon .GlobalEnv — closure's lexical scope is main process,
            # but everywhere() created these in the daemon's global env.
            GTF_annotation <- get("GTF_annotation", envir = .GlobalEnv)
            for_ORFquant_data <- get("for_ORFquant_data", envir = .GlobalEnv)
            genome_seq <- get("genome_seq", envir = .GlobalEnv)

            t_start <- as.numeric(Sys.time())
            pid <- Sys.getpid()
            outfile <- file.path(RESULT_DIR, sprintf("gene_%06d.rds", g))
            tryCatch({
                gen_region <- genes_red[g]
                chr_name <- as.character(seqnames(gen_region))

                code_id <- GTF_annotation$genetic_codes$genetic_code[
                    rownames(GTF_annotation$genetic_codes) == chr_name
                ]
                genetcd <- getGeneticCode(code_id)

                if (canonical_start_only) {
                    attributes(genetcd)$alt_init_codons <- names(
                        which(genetcd == "M")
                    )
                }

                res <- ORFquant(
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
                saveRDS(res, outfile, compress = FALSE)
                t_end <- as.numeric(Sys.time())
                is_ok <- !is.null(res) && !(is.list(res) && length(res) == 0)
                list(
                    path = outfile,
                    ok = is_ok,
                    n = if (is_ok) length(res) else 0L,
                    gene_idx = g,
                    pid = pid,
                    t_start = t_start,
                    t_end = t_end,
                    elapsed_s = round(t_end - t_start, 3)
                )
            }, error = function(e) {
                t_end <- as.numeric(Sys.time())
                message(sprintf(
                    "\n[mirai v3] Gene region %d (%s) error: %s",
                    g, as.character(genes_red[g]), conditionMessage(e)
                ))
                # Write error marker
                saveRDS(structure(list(error=conditionMessage(e)), class="orfquant_gene_error"), outfile, compress=FALSE)
                list(
                    path = outfile,
                    ok = FALSE,
                    n = 0L,
                    gene_idx = g,
                    pid = pid,
                    t_start = t_start,
                    t_end = t_end,
                    elapsed_s = round(t_end - t_start, 3),
                    error = conditionMessage(e)
                )
            })
        },
        RESULT_DIR = result_dir
    )[]

    cat(sprintf("[mirai v3] mirai_map returned. %s\n", date()))

    # ---- Step 6: Filter results without full deserialization (Task 01) ----
    paths <- vapply(result_paths, function(x) x$path, character(1L))
    is_valid <- vapply(result_paths, function(x) isTRUE(x$ok), logical(1L))

    # Task 03: if ORFQ_GANTT_FILE or ORFQ_KEEP_RESULT_DIR is set, output daemon gantt profile
    gantt_file <- Sys.getenv("ORFQ_GANTT_FILE", "")
    if (nzchar(gantt_file) || nzchar(Sys.getenv("ORFQ_KEEP_RESULT_DIR", ""))) {
        target_gantt <- if (nzchar(gantt_file)) gantt_file else file.path(result_dir, "daemon_gantt.tsv")
        tryCatch({
            gantt_df <- data.frame(
                gene_idx = vapply(result_paths, function(x) if (!is.null(x$gene_idx)) as.integer(x$gene_idx) else NA_integer_, integer(1L)),
                pid = vapply(result_paths, function(x) if (!is.null(x$pid)) as.integer(x$pid) else NA_integer_, integer(1L)),
                t_start = vapply(result_paths, function(x) if (!is.null(x$t_start)) as.numeric(x$t_start) else NA_real_, numeric(1L)),
                t_end = vapply(result_paths, function(x) if (!is.null(x$t_end)) as.numeric(x$t_end) else NA_real_, numeric(1L)),
                elapsed_s = vapply(result_paths, function(x) if (!is.null(x$elapsed_s)) as.numeric(x$elapsed_s) else NA_real_, numeric(1L)),
                n_orfs = vapply(result_paths, function(x) if (!is.null(x$n)) as.integer(x$n) else 0L, integer(1L)),
                ok = is_valid,
                stringsAsFactors = FALSE
            )
            write.table(gantt_df, target_gantt, sep = "\t", quote = FALSE, row.names = FALSE)
            cat(sprintf("[mirai v3] Gantt profile written to: %s\n", target_gantt))
        }, error = function(e) {
            warning("[mirai v3] Failed to write daemon gantt TSV: ", conditionMessage(e))
        })
    }

    n_failed <- sum(!is_valid)
    n_succeeded <- sum(is_valid)
    if (n_failed > 0) {
        cat(sprintf(
            "\n[mirai v3] %d / %d gene regions failed or empty, %d succeeded\n",
            n_failed, n_regions, n_succeeded
        ))
    }

    cat(sprintf("[mirai v3] Processing complete. %d regions → %d RDS files in %s. %s\n",
        n_regions, n_succeeded, result_dir, date()))

    # Return validated file paths + result_dir for later cleanup
    result_paths <- paths
    attr(result_paths, "result_dir") <- result_dir
    attr(result_paths, "is_valid") <- is_valid
    result_paths <- result_paths[is_valid]
    # Subsetting drops attributes; re-attach
    attr(result_paths, "result_dir") <- result_dir
    attr(result_paths, "is_valid") <- is_valid[is_valid]
    return(result_paths)
}
# DEBUG PATCH - remove after testing
