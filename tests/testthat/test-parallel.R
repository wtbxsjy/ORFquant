library(testthat)

# Helper: find test_data relative to the package root
test_data_path <- function(file) {
    pkg_root <- tryCatch(
        rprojroot::find_package_root_file(),
        error = function(e) "."
    )
    file.path(pkg_root, "test_data", file)
}

# Helper: cheaply confirm an annotation file exists and load it
load_test_annot <- function(name) {
    p <- test_data_path(name)
    skip_if_not(file.exists(p), paste0("test_data/", name, " not found"))
    load_env <- new.env(parent = emptyenv())
    load(p, envir = load_env)
    get(ls(load_env)[1], envir = load_env)
}

# ─────────────────────────────────────────────────────────────────────────────
# 1. .orfquant_genome_ref: extracts path and circular ranges from FaFile_Circ
# ─────────────────────────────────────────────────────────────────────────────
test_that(".orfquant_genome_ref returns NULL for non-FaFile", {
    skip_if_not_installed("ORFquant")
    ref <- ORFquant:::.orfquant_genome_ref(DNAStringSet(c(chr1 = "ACGTACGT")))
    expect_null(ref)
})

test_that(".orfquant_genome_ref returns ORFquantGenomeRef from FaFile", {
    skip_if_not_installed("Rsamtools")
    fa <- withr::local_tempfile(fileext = ".fa")
    Biostrings::writeXStringSet(
        Biostrings::DNAStringSet(c(chr1 = "ACGTACGT", chrM = "AAAAGGGG")),
        filepath = fa
    )
    Rsamtools::indexFa(fa)
    fa_obj <- ORFquant:::FaFile_Circ(
        Rsamtools::FaFile(fa),
        circularRanges = "chrM"
    )
    ref <- ORFquant:::.orfquant_genome_ref(fa_obj)
    expect_s3_class(ref, "ORFquantGenomeRef")
    expect_equal(ref$type, "fasta")
    expect_true(file.exists(ref$path))
    expect_equal(ref$circularRanges, "chrM")
})

# ─────────────────────────────────────────────────────────────────────────────
# 2. .orfquant_open_genome_ref: reconstructs a FaFile_Circ from a ref
# ─────────────────────────────────────────────────────────────────────────────
test_that(".orfquant_open_genome_ref creates usable FaFile_Circ", {
    skip_if_not_installed("Rsamtools")
    fa <- withr::local_tempfile(fileext = ".fa")
    Biostrings::writeXStringSet(
        Biostrings::DNAStringSet(c(chr1 = "ACGTACGT", chrM = "AAAAGGGG")),
        filepath = fa
    )
    Rsamtools::indexFa(fa)
    ref <- structure(
        list(type = "fasta", path = normalizePath(fa), circularRanges = "chrM"),
        class = "ORFquantGenomeRef"
    )
    fa_obj <- ORFquant:::.orfquant_open_genome_ref(ref)
    expect_true(inherits(fa_obj, "FaFile_Circ"))
    seq <- Biostrings::getSeq(fa_obj, Rsamtools::scanFaIndex(fa_obj))
    expect_equal(length(seq), 2L)
    ORFquant:::.orfquant_close_genome(fa_obj)
})

# ─────────────────────────────────────────────────────────────────────────────
# 3. load_annotation: FaFile_Circ annotation creates genome_ref, keeps FaFile
# ─────────────────────────────────────────────────────────────────────────────
test_that("load_annotation with FaFile annotation sets genome_ref", {
    skip_if_not_installed("ORFquant")
    annot_file <- test_data_path("mini_Rannot_fafile")
    skip_if_not(file.exists(annot_file), "mini_Rannot_fafile not found")

    # load_annotation uses <<- which super-assigns into the calling env's parent
    # chain. We call it directly (not inside local()) so <<- finds this frame's
    # bindings and updates them.
    GTF_annotation <- NULL
    genome_seq     <- NULL
    load_annotation(annot_file)
    ann <- GTF_annotation
    seq <- genome_seq

    # genome_ref must now be present when annotation contained a FaFile genome
    expect_false(is.null(ann$genome_ref),
        info = "genome_ref should be set for FaFile annotation")
    expect_equal(ann$genome_ref$type, "fasta")

    # genome_seq is still the live FaFile (each worker opens its own copy)
    expect_true(inherits(seq, "FaFile"))
})

# ─────────────────────────────────────────────────────────────────────────────
# 4. serial run on mini FaFile annotation produces results without errors
#    Uses gene_name filter to restrict to the first few genes in the test data
# ─────────────────────────────────────────────────────────────────────────────
test_that("run_ORFquant serial finishes on FaFile annotation", {
    skip_if_not_installed("ORFquant")
    annot_file <- test_data_path("mini_Rannot_fafile")
    psites_file <- test_data_path("mini_for_ORFquant")
    skip_if_not(
        file.exists(annot_file) && file.exists(psites_file),
        "test_data files not found"
    )
    skip_on_ci()

    # Load annotation to get available gene ids on chr1
    load_env <- new.env(parent = emptyenv())
    load(annot_file, envir = load_env)
    ann <- get(ls(load_env)[1], envir = load_env)
    chr1_gene_ids <- names(ann$genes)[
        as.character(GenomicRanges::seqnames(ann$genes)) == "1"
    ]
    skip_if(length(chr1_gene_ids) == 0, "No chr1 genes found in test data")

    withr::with_tempdir({
        msgs  <- character()
        warns <- character()
        result <- tryCatch(
            withCallingHandlers(
                withCallingHandlers(
                    run_ORFquant(
                        for_ORFquant_file   = psites_file,
                        annotation_file     = annot_file,
                        n_cores             = 1,
                        prefix              = "test_serial",
                        gene_id             = chr1_gene_ids[seq_len(min(3, length(chr1_gene_ids)))],
                        write_temp_files    = FALSE,
                        write_GTF_file      = FALSE,
                        write_protein_fasta = FALSE,
                        interactive         = FALSE,
                        parallel_backend    = "serial"
                    ),
                    message = function(m) {
                        msgs <<- c(msgs, conditionMessage(m))
                        invokeRestart("muffleMessage")
                    }
                ),
                warning = function(w) {
                    warns <<- c(warns, conditionMessage(w))
                    invokeRestart("muffleWarning")
                }
            ),
            error = function(e) {
                # If no ORFs found it will error; treat as a skip, not failure
                msg <- conditionMessage(e)
                if (grepl("Not enough P_sites|Incorrect gene", msg)) {
                    skip(paste("Insufficient data in test_data:", msg))
                }
                stop(e)
            }
        )
        # No FaFile finalizer errors in warnings
        finalizer_err <- any(grepl("finalize|non-function", warns, ignore.case = TRUE))
        expect_false(finalizer_err,
            info = paste("Unexpected finalizer warnings:", paste(warns, collapse = "; ")))
        expect_true(is.list(result))
    })
})
