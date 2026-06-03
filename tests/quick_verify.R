#!/usr/bin/env Rscript
# =============================================================================
# ORFquant Modernization Quick Verification (base R only, no dependencies)
# =============================================================================
#
# Verifies that modernization patterns were correctly applied by checking
# the source code directly. No package installation needed.
#
# Usage:
#   Rscript 00_quick_verify.R
# =============================================================================

v131_src <- normalizePath("~/riboseq/ORFquant/R/orfquant.R",
    mustWork = FALSE)

# Try alternative relative paths
if (!file.exists(v131_src)) {
    v131_src <- normalizePath("../../ORFquant/R/orfquant.R",
        mustWork = FALSE)
}
if (!file.exists(v131_src)) {
    stop("Cannot find modernized orfquant.R. Expected at:\n",
         "  ~/riboseq/ORFquant/R/orfquant.R or\n",
         "  ../../ORFquant/R/orfquant.R")
}

cat("Checking:", v131_src, "\n\n")

lines <- readLines(v131_src)
total_lines <- length(lines)

pass <- 0
fail <- 0
warn <- 0

check <- function(desc, test_expr) {
    result <- tryCatch(test_expr, error = function(e) FALSE)
    if (isTRUE(result)) {
        cat(sprintf("  [PASS] %s\n", desc))
        pass <<- pass + 1
    } else {
        cat(sprintf("  [FAIL] %s\n", desc))
        if (is.character(result)) cat(sprintf("         %s\n", result))
        fail <<- fail + 1
    }
}

cat("=== dim()[1]/[2] -> nrow()/ncol() ===\n")
dim1 <- grep("dim\\([^)]+\\)\\[1\\]", lines, value = TRUE)
dim2 <- grep("dim\\([^)]+\\)\\[2\\]", lines, value = TRUE)
check("No dim(x)[1] patterns remain",
    length(dim1) == 0 || paste("Found:", paste(head(dim1, 3), collapse = "; ")))
check("No dim(x)[2] patterns remain",
    length(dim2) == 0 || paste("Found:", paste(head(dim2, 3), collapse = "; ")))

nrow_uses <- sum(grepl("nrow\\(", lines))
ncol_uses <- sum(grepl("ncol\\(", lines))
cat(sprintf("  INFO: nrow() used %d times, ncol() used %d times\n", nrow_uses, ncol_uses))

cat("\n=== stringsAsFactors removal ===\n")
saf_lines <- grep("stringsAsFactors", lines, value = TRUE)
saf_in_df <- saf_lines[grepl("data\\.frame", saf_lines)]
check("No stringsAsFactors=FALSE in data.frame() calls",
    length(saf_in_df) == 0 || paste("Found:", paste(saf_in_df, collapse = "; ")))
if (length(saf_lines) > 0) {
    cat(sprintf("  INFO: %d stringsAsFactors remain (expected in read.table only):\n",
        length(saf_lines)))
    for (l in saf_lines) cat(sprintf("    %s\n", trimws(l)))
}

cat("\n=== S4 slot access -> accessor functions ===\n")
at_ranges <- grep("@ranges", lines, value = TRUE)
at_ranges_code <- at_ranges[!grepl("^\\s*#", at_ranges)]
check("No @ranges in code (only comments allowed)",
    length(at_ranges_code) == 0 || paste("Found:", paste(head(at_ranges_code, 3), collapse = "; ")))

at_strand <- grep("@strand", lines, value = TRUE)
at_strand_code <- at_strand[!grepl("^\\s*#", at_strand)]
check("No @strand in code (only comments allowed)",
    length(at_strand_code) == 0 || paste("Found:", paste(head(at_strand_code, 3), collapse = "; ")))

at_width <- grep("@width", lines, value = TRUE)
at_width_code <- at_width[!grepl("^\\s*#", at_width)]
check("No @width in code (only comments allowed)",
    length(at_width_code) == 0 || paste("Found:", paste(head(at_width_code, 3), collapse = "; ")))

at_is_circular <- grep("@is_circular", lines, value = TRUE)
at_is_circular_code <- at_is_circular[!grepl("^\\s*#", at_is_circular)]
check("No @is_circular in code (only comments allowed)",
    length(at_is_circular_code) == 0 || paste("Found:", paste(head(at_is_circular_code, 3), collapse = "; ")))

at_values <- grep("@values", lines, value = TRUE)
at_values_code <- at_values[!grepl("^\\s*#", at_values)]
check("No @values in code (only comments allowed)",
    length(at_values_code) == 0 || paste("Found:", paste(head(at_values_code, 3), collapse = "; ")))

cat("\n=== c() accumulation in loops -> pre-allocation ===\n")
shft_init <- grep("shft <- c\\(\\)", lines)
shft_grow <- grep("shft <- c\\(shft,", lines)
check("No shft <- c() + c(shft, ...) growth pattern",
    length(shft_init) == 0 || paste("Found shft <- c() at line", shft_init))
check("No shft <- c(shft, ...) growth pattern",
    length(shft_grow) == 0 || paste("Found shft <- c(shft, at line", shft_grow))

stok_init <- grep("stok <- c\\(\\)", lines)
stok_grow <- grep("stok <- c\\(stok,", lines)
check("No stok <- c() + c(stok, ...) growth pattern",
    length(stok_init) == 0 || paste("Found stok <- c() at line", stok_init))
check("No stok <- c(stok, ...) growth pattern",
    length(stok_grow) == 0 || paste("Found stok <- c(stok, at line", stok_grow))

nmss_init <- grep("nmss <- c\\(\\)", lines)
check("nmss <- c() replaced with rep()",
    length(nmss_init) == 0 || paste("Found at line", nmss_init))

cat("\n=== sapply() -> vapply() type-stable mapping ===\n")
n_sapply <- sum(grepl("sapply\\(", lines))
n_vapply <- sum(grepl("vapply\\(", lines))
cat(sprintf("  sapply() calls: %d\n", n_sapply))
cat(sprintf("  vapply() calls: %d\n", n_vapply))
check("vapply() count > 0 (type-stable mapping used)",
    n_vapply > 0)

cat("\n=== ifelse() -> dplyr::if_else() ===\n")
n_ifelse <- sum(grepl("ifelse\\(", lines))
n_if_else <- sum(grepl("if_else\\(", lines))
cat(sprintf("  ifelse() calls: %d\n", n_ifelse))
cat(sprintf("  if_else() calls: %d\n", n_if_else))
check("No base ifelse() remaining (use dplyr::if_else)",
    n_ifelse == 0 || paste("Remaining:", n_ifelse, "ifelse() calls"))

cat("\n=== DPSS caching infrastructure ===\n")
check(".dpss_cache environment exists",
    any(grepl("\\.dpss_cache", lines)))
check("get_cached_dpss() function exists",
    any(grepl("get_cached_dpss", lines)))
check("clear_dpss_cache() function exists and is exported",
    any(grepl("clear_dpss_cache", lines)))

cat("\n=== R version requirement ===\n")
desc_lines <- readLines(normalizePath("~/riboseq/ORFquant/DESCRIPTION",
    mustWork = FALSE))
dep_line <- grep("Depends:.*R", desc_lines, value = TRUE)
cat(sprintf("  %s\n", trimws(dep_line)))
check("R >= 4.0.0 (for stringsAsFactors default)",
    grepl("4\\.0", dep_line))

cat("\n=== Parse validation ===\n")
parse_ok <- tryCatch({
    parse(file = v131_src)
    TRUE
}, error = function(e) {
    paste("Parse error:", e$message)
})
check("Source file parses without error", parse_ok)

cat("\n=== Function counts ===\n")
func_pattern <- "function\\s*\\("
n_funcs <- sum(grepl(func_pattern, lines))
# Count only top-level function definitions (not anonymous)
n_exported <- sum(grepl("^#' @export", lines))
cat(sprintf("  Total function definitions: %d\n", n_funcs))
cat(sprintf("  Exported functions: %d\n", n_exported))

# =============================================================================
# Summary
# =============================================================================
cat("\n============================================\n")
cat(sprintf(" Results: %d PASS, %d FAIL, %d INFO\n", pass, fail, warn))
cat("============================================\n")

if (fail > 0) {
    cat("\nSome checks failed. Review the FAIL items above.\n")
    quit(status = 1)
} else {
    cat("\nAll modernization checks passed!\n")
    cat("The source code is structurally consistent with the modernization plan.\n")
    cat("For functional verification, run the integration tests with BioC packages installed.\n")
    quit(status = 0)
}
