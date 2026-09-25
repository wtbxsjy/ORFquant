library(testthat)
suppressPackageStartupMessages({
    library(GenomicRanges)
    library(Biostrings)
})

# Previous get_orfs() implementation (one translate() per frame, linear scan
# of the stop codons for every start), kept here as the reference.
legacy_get_orfs <- function(tx_name, sequence, genetic_code_table) {
    list_frames <- list()
    for (u in 0:2) {
        pept <- unlist(strsplit(
            as.character(suppressWarnings(translate(
                subseq(sequence, start = u + 1),
                genetic.code = genetic_code_table,
                if.fuzzy.codon = "solve"
            ))),
            split = ""
        ))
        start_pos <- ((seq_along(pept))[pept == "M"]) * 3
        start_pos <- if (length(start_pos) > 0) start_pos + u - 2 else NA
        stop_pos <- ((seq_along(pept))[pept == "*"]) * 3 - 3
        stop_pos <- if (length(stop_pos) > 0) stop_pos + u else NA
        st2vect <- c()
        for (h in seq_along(start_pos)) {
            diff <- stop_pos - start_pos[h]
            diff <- diff[diff > 0]
            st2vect[h] <- if (length(diff) > 0) start_pos[h] + min(diff) else NA
        }
        st_st <- data.frame(cbind(start_pos, st2vect))
        st_st <- st_st[!is.na(st_st[, 1]), ]
        st_st <- st_st[!is.na(st_st[, 2]), ]
        if (nrow(st_st) == 0) {
            list_frames[[paste("frame", u, sep = "_")]] <- GRanges()
            next
        }
        gra_orf <- GRanges(
            seqnames = paste(tx_name, "frame", u, sep = "_"),
            strand = "+",
            ranges = IRanges(start = st_st[, 1], end = st_st[, 2])
        )
        gra_orf$type <- "ORF"
        gra_orf$score <- 1
        list_frames[[paste("frame", u, sep = "_")]] <- gra_orf
    }
    GRangesList(list_frames)
}

test_that("get_orfs matches the previous implementation", {
    set.seed(11)
    gc_std <- getGeneticCode("1")
    gc_can <- gc_std
    attributes(gc_can)$alt_init_codons <- names(which(gc_can == "M"))
    for (len in c(2, 3, 5, 60, 301, 1500)) {
        for (i in 1:3) {
            s <- DNAString(paste(
                sample(c("A", "C", "G", "T"), len, replace = TRUE),
                collapse = ""
            ))
            for (gc in list(gc_std, gc_can)) {
                expect_identical(
                    get_orfs("tx1", s, genetic_code_table = gc),
                    legacy_get_orfs("tx1", s, gc)
                )
            }
        }
    }
})

test_that("get_orfs handles sequences without start or stop codons", {
    gc <- getGeneticCode("1")
    no_start <- DNAString(strrep("CCC", 20))
    expect_identical(
        get_orfs("tx", no_start, genetic_code_table = gc),
        legacy_get_orfs("tx", no_start, gc)
    )
    no_stop <- DNAString(paste0("ATG", strrep("CCC", 20)))
    expect_identical(
        get_orfs("tx", no_stop, genetic_code_table = gc),
        legacy_get_orfs("tx", no_stop, gc)
    )
})

test_that("region index reproduces x[x %over% region] for every region", {
    set.seed(3)
    ps <- GRanges(
        rep(c("1", "2"), each = 300),
        IRanges(sample(1:5000, 600, replace = TRUE), width = 1),
        strand = sample(c("+", "-"), 600, replace = TRUE),
        score = sample(1:20, 600, replace = TRUE)
    )
    junc <- GRanges(
        "1",
        IRanges(c(100, 2000, 4000), width = 300),
        strand = "+",
        reads = c(1, 2, 3)
    )
    fo <- list(
        P_sites_all = ps,
        P_sites_uniq = ps[ps$score > 5],
        P_sites_uniq_mm = ps[0],
        junctions = junc
    )
    regions <- GRanges(
        c("1", "1", "2", "2"),
        IRanges(c(1, 3000, 10, 6000), c(2500, 4999, 1000, 7000)),
        strand = c("+", "-", "+", "+")
    )
    idx <- ORFquant:::.orfquant_region_index(fo, regions)
    for (g in seq_along(regions)) {
        sub <- ORFquant:::.orfquant_region_data(fo, idx, g)
        for (nm in names(fo)) {
            expect_identical(
                sub[[nm]],
                fo[[nm]][fo[[nm]] %over% regions[g]]
            )
        }
    }
})

test_that("select_start frame statistics are computed per ORF", {
    # in-frame signal on codon position 1 of an ORF at 1-30
    cov <- Rle(c(rep(c(5L, 0L, 0L), 10), rep(0L, 30)))
    orfs <- GRanges("tx", IRanges(c(1, 4), c(30, 30)), "+")
    res <- select_start(orfs, cov, cutoff = NA, cutoff_ave = .5)
    expect_s4_class(res, "GRanges")
    expect_equal(length(res), 1L)
    expect_equal(start(res), 1L)
    expect_equal(res$ave_pct_fr, 100)
    expect_equal(res$pct_fr, 100)
})

test_that("select_txs accepts exonicParts() bins with extra columns", {
    # exonicParts() (Bioc >= 3.20) adds tx_id/exon_id/exon_name/exon_rank,
    # which previously broke the column alignment with the junction counts.
    bins <- GRanges(
        "1",
        IRanges(c(1, 201), c(100, 300)),
        "+",
        tx_id = IntegerList(1, 1),
        tx_name = CharacterList("t1", "t1"),
        gene_id = CharacterList("g1", "g1"),
        exon_id = IntegerList(1, 2),
        exon_name = CharacterList("e1", "e2"),
        exon_rank = IntegerList(1, 2)
    )
    junctions <- GRanges(
        "1", IRanges(101, 200), "+",
        type = "J",
        tx_name = CharacterList("t1"),
        gene_id = CharacterList("g1"),
        reads = 3,
        unique_reads = 3
    )
    ps <- GRanges("1", IRanges(c(10, 13, 16, 220, 223), width = 1), "+",
        score = 2L)
    res <- select_txs(
        region = GRanges("1", IRanges(1, 300), "+"),
        annotation = list(exons_bins = bins),
        P_sites = ps,
        P_sites_uniq = ps,
        junction_counts = junctions
    )
    expect_s4_class(res, "GRanges")
    expect_equal(sort(unique(res$type)), c("E", "J"))
    expect_equal(sum(res$reads[res$type == "E"]), 10)
})
