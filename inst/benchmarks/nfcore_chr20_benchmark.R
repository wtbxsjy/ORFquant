#!/usr/bin/env Rscript
# =============================================================================
# ORFquant benchmark on the nf-core/riboseq test data (human chr20)
# =============================================================================
#
# Data (nf-core/test-datasets, branch "modules",
# data/genomics/homo_sapiens/riboseq_expression):
#   Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz
#   Homo_sapiens.GRCh38.111_chr20.gtf
#   aligned_reads/SRX11780887_chr20.bam(.bai)   (Ribo-seq, nf-core/riboseq)
#
#   git clone --depth 1 --filter=blob:none --sparse --branch modules \
#       https://github.com/nf-core/test-datasets.git
#   git -C test-datasets sparse-checkout set \
#       data/genomics/homo_sapiens/riboseq_expression
#
# Usage:
#   Rscript nfcore_chr20_benchmark.R prepare <data_dir> <work_dir>
#   Rscript nfcore_chr20_benchmark.R run <work_dir> <lib> <tag> <A|B> \
#       [end_pos|Inf] [n_cores] [backend]
#   Rscript nfcore_chr20_benchmark.R compare <res1.rds> <res2.rds>
#
# Benchmark A: real P-sites from the (down-sampled) nf-core BAM.
# Benchmark B: same annotation, simulated high-depth P-sites (3-nt periodic
#              signal on CDSs of ~70% of protein-coding genes, 1-2 isoforms
#              each, plus background and junction reads).  Exercises the
#              quantification / annotation code paths that the sparse test
#              BAM rarely reaches.
# =============================================================================

args <- commandArgs(TRUE)
mode <- args[1]

prepare <- function(data_dir, work_dir) {
    suppressPackageStartupMessages({
        library(ORFquant)
        library(Rsamtools)
    })
    dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
    owd <- setwd(work_dir)
    on.exit(setwd(owd))
    fa <- "chr20.fa"
    if (!file.exists(fa)) {
        con <- gzfile(file.path(
            data_dir, "Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz"
        ))
        writeLines(readLines(con), fa)
        close(con)
    }
    indexFa(fa)
    prepare_annotation_files(
        annotation_directory = "annot",
        gtf_file = file.path(data_dir, "Homo_sapiens.GRCh38.111_chr20.gtf"),
        scientific_name = "Homo.sapiens",
        annotation_name = "chr20",
        genome_seq = fa,
        export_bed_tables_TxDb = FALSE
    )
    annot <- list.files("annot", pattern = "Rannot$", full.names = TRUE)
    rl <- data.frame(read_length = 26:34, cutoff = 12, comp = "nucl")
    write.table(rl, "rl_cutoff.txt", sep = "\t", quote = FALSE,
        row.names = FALSE)
    prepare_for_ORFquant(
        annotation_file = annot,
        bam_file = file.path(data_dir, "aligned_reads",
            "SRX11780887_chr20.bam"),
        path_to_rl_cutoff_file = "rl_cutoff.txt",
        dest_name = "A"
    )
    simulate_psites(annot, "B_for_ORFquant")
}

simulate_psites <- function(annot, out) {
    set.seed(1)
    load_annotation(annot)
    tr <- GTF_annotation$trann
    cds <- GTF_annotation$cds_txs
    pc <- names(cds)[tr$transcript_biotype[
        match(names(cds), tr$transcript_id)
    ] == "protein_coding"]
    genes <- tr$gene_id[match(pc, tr$transcript_id)]
    sel <- unlist(lapply(split(pc, genes), function(x) {
        head(sample(x), sample(1:2, 1))
    }))
    sel <- sel[runif(length(sel)) < .7]
    ps <- list()
    for (tx in sel) {
        cd <- cds[[tx]]
        ncod <- sum(width(cd)) %/% 3
        if (ncod < 20) next
        lam <- rlnorm(1, log(1.5), 1)
        cnt <- rpois(ncod * 3, lam * rep(c(.75, .12, .13), ncod))
        pos <- which(cnt > 0)
        g <- mapFromTranscripts(
            GRanges(tx, IRanges(pos, width = 1)),
            GRangesList(setNames(list(cd), tx))
        )
        mcols(g) <- DataFrame(score = cnt[pos])
        ps[[tx]] <- g
        ex <- GTF_annotation$exons_txs[[tx]]
        bg <- ex[sample(length(ex), 1)]
        nb <- rpois(1, 5)
        if (nb > 0) {
            b <- GRanges(seqnames(bg),
                IRanges(sample(start(bg):end(bg), nb, TRUE), width = 1),
                strand(bg), score = 1L)
            ps[[paste0(tx, "_bg")]] <- b
        }
    }
    P <- sort(unlist(GRangesList(unname(ps))))
    key <- paste(start(P), strand(P))
    keep <- !duplicated(key)
    P <- GRanges(seqnames(P)[keep], ranges(P)[keep], strand(P)[keep],
        score = as.integer(tapply(P$score,
            factor(key, levels = unique(key)), sum)))
    seqlevels(P) <- seqlevels(GTF_annotation$seqinfo)
    seqinfo(P) <- GTF_annotation$seqinfo
    Pu <- P
    Pu$score <- rbinom(length(P), P$score, .9)
    Pu <- Pu[Pu$score > 0]
    Pmm <- Pu[runif(length(Pu)) < .05]
    j <- GTF_annotation$junctions
    jt <- vapply(j$tx_name, function(x) any(x %in% sel), logical(1))
    j$reads <- 0
    j$reads[jt] <- rpois(sum(jt), 8)
    j$unique_reads <- rbinom(length(j), j$reads, .9)
    for_ORFquant <- list(P_sites_all = P, P_sites_uniq = Pu,
        P_sites_uniq_mm = Pmm, junctions = j)
    save(for_ORFquant, file = out)
}

run <- function(work_dir, lib, tag, which, end_pos = Inf, n_cores = 1,
    backend = "mclapply") {
    suppressPackageStartupMessages(library(ORFquant, lib.loc = lib))
    pf <- file.path(work_dir, if (which == "A") "A_for_ORFquant" else
        "B_for_ORFquant")
    annot <- list.files(file.path(work_dir, "annot"), pattern = "Rannot$",
        full.names = TRUE)
    rg <- if (is.finite(end_pos)) {
        GenomicRanges::GRanges("20", IRanges::IRanges(1, end_pos))
    } else {
        NA
    }
    t <- system.time(res <- run_ORFquant(
        for_ORFquant_file = pf, annotation_file = annot,
        genomic_region = rg, n_cores = n_cores, parallel_backend = backend,
        prefix = file.path(tempdir(), tag), write_temp_files = FALSE,
        interactive = FALSE
    ))
    res$psite_data_file <- NULL
    saveRDS(res, paste0("res_", tag, ".rds"))
    cat(sprintf("%s elapsed=%.1fs ORFs_tx=%d ORFs_gen=%d selected_txs=%d\n",
        tag, t[["elapsed"]], length(res$ORFs_tx), length(res$ORFs_gen),
        length(res$selected_txs)))
}

compare <- function(f1, f2) {
    suppressPackageStartupMessages(library(GenomicRanges))
    a <- readRDS(f1)
    b <- readRDS(f2)
    for (nm in union(names(a), names(b))) {
        r <- all.equal(a[[nm]], b[[nm]])
        cat(sprintf("%-24s %5d vs %5d  all.equal=%s\n", nm,
            length(a[[nm]]), length(b[[nm]]),
            if (isTRUE(r)) "TRUE" else paste(head(r, 2), collapse = " | ")))
    }
}

switch(mode,
    prepare = prepare(args[2], args[3]),
    run = run(args[2], args[3], args[4], args[5],
        if (length(args) > 5) as.numeric(args[6]) else Inf,
        if (length(args) > 6) as.integer(args[7]) else 1L,
        if (length(args) > 7) args[8] else "mclapply"),
    compare = compare(args[2], args[3]),
    stop("mode must be one of prepare / run / compare")
)
