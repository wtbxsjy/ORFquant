#!/usr/bin/env Rscript
# Benchmark script to compare ORFquant v1.2 (baseline) vs v1.3 (optimized)
# Usage: Rscript benchmark_orfquant.R

message("=== ORFquant Performance Benchmark ===")
message("Comparing v1.2.0 (baseline) vs v1.3.0 (optimized)")
message("")

# Load required packages
suppressPackageStartupMessages({
    library(microbenchmark)
    library(multitaper)
})

# Test 1: DPSS Caching Performance
message("Test 1: DPSS Caching Performance")
message("---")

# Simulate different ORF lengths (common scenario: many ORFs with similar lengths)
set.seed(42)
orf_lengths <- sample(50:500, 100, replace = TRUE)

# Function to simulate old behavior (no caching)
benchmark_dpss_uncached <- function(lengths, k = 24, nw = 12) {
    results <- list()
    for (len in lengths) {
        results[[length(results) + 1]] <- dpss(n = len, k = k, nw = nw)
    }
    invisible(results)
}

# Try to load v1.3 caching function
tryCatch(
    {
        # Source the get_cached_dpss function
        # In actual use, this comes from the ORFquant package
        .dpss_cache <- new.env(hash = TRUE, parent = emptyenv())

        get_cached_dpss <- function(n, k, nw) {
            key <- paste0(n, "_", k, "_", nw)
            if (!exists(key, envir = .dpss_cache, inherits = FALSE)) {
                assign(key, dpss(n = n, k = k, nw = nw), envir = .dpss_cache)
            }
            get(key, envir = .dpss_cache, inherits = FALSE)
        }

        benchmark_dpss_cached <- function(lengths, k = 24, nw = 12) {
            results <- list()
            for (len in lengths) {
                results[[length(results) + 1]] <- get_cached_dpss(
                    n = len,
                    k = k,
                    nw = nw
                )
            }
            invisible(results)
        }

        # Run benchmark
        message("Running timing comparison (100 ORF lengths)...")
        message("")

        # First run: both methods start fresh
        rm(list = ls(.dpss_cache), envir = .dpss_cache)

        # Time uncached method
        time_uncached <- system.time({
            for (i in 1:3) {
                benchmark_dpss_uncached(orf_lengths)
            }
        })["elapsed"]

        # Time cached method (includes cache population on first pass)
        rm(list = ls(.dpss_cache), envir = .dpss_cache)
        time_cached <- system.time({
            for (i in 1:3) {
                benchmark_dpss_cached(orf_lengths)
            }
        })["elapsed"]

        message(sprintf(
            "Uncached (v1.2 behavior): %.3f seconds",
            time_uncached
        ))
        message(sprintf("Cached (v1.3 optimized):  %.3f seconds", time_cached))
        message(sprintf("Speedup factor: %.2fx", time_uncached / time_cached))
        message("")

        # Test with repeated lengths (simulates real-world scenario)
        message("Test 2: Repeated Lengths (simulates many similar ORFs)")
        message("---")

        repeated_lengths <- rep(c(100, 150, 200, 250, 300), each = 20)

        rm(list = ls(.dpss_cache), envir = .dpss_cache)

        time_uncached_rep <- system.time({
            for (i in 1:3) {
                benchmark_dpss_uncached(repeated_lengths)
            }
        })["elapsed"]

        rm(list = ls(.dpss_cache), envir = .dpss_cache)
        time_cached_rep <- system.time({
            for (i in 1:3) {
                benchmark_dpss_cached(repeated_lengths)
            }
        })["elapsed"]

        message(sprintf(
            "Uncached (v1.2 behavior): %.3f seconds",
            time_uncached_rep
        ))
        message(sprintf(
            "Cached (v1.3 optimized):  %.3f seconds",
            time_cached_rep
        ))
        message(sprintf(
            "Speedup factor: %.2fx",
            time_uncached_rep / time_cached_rep
        ))
    },
    error = function(e) {
        message("Error running benchmark: ", conditionMessage(e))
    }
)

message("")
message("=== Benchmark Complete ===")
