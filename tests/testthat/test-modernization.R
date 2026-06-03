# =============================================================================
# ORFquant v1.3.2 Modernization Tests
# =============================================================================
#
# These tests verify that the codebase follows modern R conventions:
# - nrow()/ncol() instead of dim()[1]/dim()[2]
# - accessor functions (ranges(), strand(), width()) instead of @ slot access
# - type-stable vapply() instead of sapply() where appropriate
# - pre-allocation instead of c() accumulation in loops
# - dplyr::if_else() instead of base ifelse()
# - No redundant stringsAsFactors = FALSE in data.frame() calls
# - R >= 4.0.0 dependency
# =============================================================================

library(testthat)

# ---------------------------------------------------------------------------
# Helper: find the package source file
# ---------------------------------------------------------------------------
pkg_src <- function() {
    r_dir <- system.file("R", package = "ORFquant")
    if (r_dir == "") {
        # During development, look relative to working directory
        r_dir <- "../../R"
    }
    file.path(r_dir, "orfquant.R")
}

# ---------------------------------------------------------------------------
# Test 1: R version requirement
# ---------------------------------------------------------------------------
test_that("R version requirement is >= 4.0.0", {
    desc <- readLines(system.file("DESCRIPTION", package = "ORFquant"))
    dep_line <- grep("Depends:.*R", desc, value = TRUE)
    expect_true(length(dep_line) > 0, info = "No R dependency line in DESCRIPTION")
    expect_true(grepl("4\\.0", dep_line),
        info = sprintf("R dependency should be >= 4.0.0, got: %s", dep_line))
})

# ---------------------------------------------------------------------------
# Test 2: Package version is 1.3.2+
# ---------------------------------------------------------------------------
test_that("Package version is 1.3.2 or later", {
    ver <- as.character(packageVersion("ORFquant"))
    expect_true(compareVersion(ver, "1.3.2") >= 0,
        info = sprintf("Expected version >= 1.3.2, got %s", ver))
})

# ---------------------------------------------------------------------------
# Test 3: dim()[] patterns are replaced with nrow()/ncol()
# ---------------------------------------------------------------------------
test_that("Source code uses nrow()/ncol() instead of dim()[1]/dim()[2]", {
    src <- readLines(pkg_src())

    # These patterns should NOT exist
    dim1 <- grep("dim\\([^)]+\\)\\[1\\]", src, value = TRUE)
    dim2 <- grep("dim\\([^)]+\\)\\[2\\]", src, value = TRUE)

    expect_equal(length(dim1), 0,
        info = sprintf("dim()[1] patterns found: %s", paste(dim1, collapse = "; ")))
    expect_equal(length(dim2), 0,
        info = sprintf("dim()[2] patterns found: %s", paste(dim2, collapse = "; ")))
})

# ---------------------------------------------------------------------------
# Test 4: No stringsAsFactors in data.frame() calls
# ---------------------------------------------------------------------------
test_that("data.frame() calls don't use redundant stringsAsFactors = FALSE", {
    src <- readLines(pkg_src())

    saf_lines <- grep("stringsAsFactors", src, value = TRUE)
    saf_in_df <- saf_lines[grepl("data\\.frame", saf_lines)]

    expect_equal(length(saf_in_df), 0,
        info = sprintf("data.frame() with stringsAsFactors: %s",
            paste(saf_in_df, collapse = "; ")))
})

# ---------------------------------------------------------------------------
# Test 5: S4 @ranges replaced with ranges() accessor
# ---------------------------------------------------------------------------
test_that("Source code uses ranges() instead of @ranges slot access", {
    src <- readLines(pkg_src())

    at_ranges <- grep("@ranges", src, value = TRUE)
    at_ranges_code <- at_ranges[!grepl("^\\s*#", at_ranges)]

    expect_equal(length(at_ranges_code), 0,
        info = sprintf("@ranges still used in: %s",
            paste(at_ranges_code, collapse = "; ")))
})

# ---------------------------------------------------------------------------
# Test 6: S4 @strand replaced with strand() accessor
# ---------------------------------------------------------------------------
test_that("Source code uses strand() instead of @strand slot access", {
    src <- readLines(pkg_src())

    at_strand <- grep("@strand", src, value = TRUE)
    at_strand_code <- at_strand[!grepl("^\\s*#", at_strand)]

    expect_equal(length(at_strand_code), 0,
        info = sprintf("@strand still used in: %s",
            paste(at_strand_code, collapse = "; ")))
})

# ---------------------------------------------------------------------------
# Test 7: S4 @width replaced with width() accessor
# ---------------------------------------------------------------------------
test_that("Source code uses width() instead of @width slot access", {
    src <- readLines(pkg_src())

    at_width <- grep("@width", src, value = TRUE)
    at_width_code <- at_width[!grepl("^\\s*#", at_width)]

    expect_equal(length(at_width_code), 0,
        info = sprintf("@width still used in: %s",
            paste(at_width_code, collapse = "; ")))
})

# ---------------------------------------------------------------------------
# Test 8: S4 @is_circular replaced with isCircular() accessor
# ---------------------------------------------------------------------------
test_that("Source code uses isCircular() instead of @is_circular slot access", {
    src <- readLines(pkg_src())

    at_ic <- grep("@is_circular", src, value = TRUE)
    at_ic_code <- at_ic[!grepl("^\\s*#", at_ic)]

    expect_equal(length(at_ic_code), 0,
        info = sprintf("@is_circular still used in: %s",
            paste(at_ic_code, collapse = "; ")))
})

# ---------------------------------------------------------------------------
# Test 9: S4 @values replaced with runValue() for Rle objects
# ---------------------------------------------------------------------------
test_that("Source code uses runValue() instead of @values for Rle objects", {
    src <- readLines(pkg_src())

    at_values <- grep("@values", src, value = TRUE)
    at_values_code <- at_values[!grepl("^\\s*#", at_values)]

    expect_equal(length(at_values_code), 0,
        info = sprintf("@values still used in: %s",
            paste(at_values_code, collapse = "; ")))
})

# ---------------------------------------------------------------------------
# Test 10: c() accumulation in simple loops replaced with pre-allocation
# ---------------------------------------------------------------------------
test_that("Simple c() accumulation loops use pre-allocation", {
    src <- readLines(pkg_src())

    # shft <- c() pattern (was O(n^2) growth)
    shft_init <- grep("shft <- c\\(\\)", src)
    shft_grow <- grep("shft <- c\\(shft,", src)
    expect_equal(length(shft_init), 0,
        info = sprintf("shft <- c() should be numeric(length(rangok))"))
    expect_equal(length(shft_grow), 0,
        info = sprintf("shft <- c(shft, ...) should be shft[i] <- ..."))

    # stok <- c() pattern
    stok_init <- grep("stok <- c\\(\\)", src)
    stok_grow <- grep("stok <- c\\(stok,", src)
    expect_equal(length(stok_init), 0,
        info = sprintf("stok <- c() should be numeric(length(shft))"))
    expect_equal(length(stok_grow), 0,
        info = sprintf("stok <- c(stok, ...) should be stok[i] <- ..."))

    # nmss <- c() replaced with rep()
    nmss_init <- grep("nmss <- c\\(\\)", src)
    expect_equal(length(nmss_init), 0,
        info = "nmss <- c() pattern should be replaced with rep()")
})

# ---------------------------------------------------------------------------
# Test 11: vapply() used for type-stable mapping
# ---------------------------------------------------------------------------
test_that("vapply() is used alongside sapply() for type-stable mapping", {
    src <- readLines(pkg_src())

    n_vapply <- sum(grepl("vapply\\(", src))
    expect_gt(n_vapply, 0,
        info = "vapply() should be used for type-stable mapping")
})

# ---------------------------------------------------------------------------
# Test 12: dplyr::if_else() replaces base ifelse()
# ---------------------------------------------------------------------------
test_that("dplyr::if_else() is used instead of base ifelse()", {
    src <- readLines(pkg_src())

    n_ifelse <- sum(grepl("ifelse\\(", src))
    n_if_else <- sum(grepl("if_else\\(", src))

    expect_equal(n_ifelse, 0,
        info = sprintf("%d base ifelse() calls remaining, should use dplyr::if_else()",
            n_ifelse))
    expect_gt(n_if_else, 0,
        info = "dplyr::if_else() should be used at least once")
})

# ---------------------------------------------------------------------------
# Test 13: DPSS caching infrastructure exists
# ---------------------------------------------------------------------------
test_that("DPSS caching environment and functions exist", {
    expect_true(exists(".dpss_cache", envir = asNamespace("ORFquant"),
        inherits = FALSE),
        info = ".dpss_cache environment missing from package namespace")
})

test_that("get_cached_dpss function is available", {
    # Check it exists in the package namespace
    ns <- asNamespace("ORFquant")
    expect_true(exists("get_cached_dpss", envir = ns, inherits = FALSE),
        info = "get_cached_dpss() missing from package")
})

test_that("clear_dpss_cache is exported", {
    expect_true("clear_dpss_cache" %in% getNamespaceExports("ORFquant"),
        info = "clear_dpss_cache() should be exported")
})

# ---------------------------------------------------------------------------
# Test 14: Source file parses without error
# ---------------------------------------------------------------------------
test_that("Package source file parses without syntax errors", {
    src_file <- pkg_src()
    if (file.exists(src_file)) {
        expect_silent(parse(file = src_file))
    } else {
        skip(sprintf("Source file not found at %s (expected in dev mode)", src_file))
    }
})

# ---------------------------------------------------------------------------
# Test 15: Exported functions have correct signatures
# ---------------------------------------------------------------------------
test_that("Core exported functions exist and are callable", {
    exported <- getNamespaceExports("ORFquant")
    core_funcs <- c("ORFquant", "run_ORFquant", "detect_translated_orfs",
        "select_txs", "calc_orf_pval", "get_orfs", "select_start",
        "prepare_annotation_files", "prepare_for_ORFquant",
        "plot_ORFquant_results", "annotate_ORFs", "clear_dpss_cache")

    for (fn in core_funcs) {
        expect_true(fn %in% exported,
            info = sprintf("'%s' should be exported", fn))
    }
})
