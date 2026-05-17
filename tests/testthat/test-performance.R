# Test for DPSS caching performance and correctness
library(testthat)

context("Performance optimizations")

test_that("DPSS cache returns correct values", {
    skip_if_not_installed("multitaper")
    library(multitaper)

    # Test that cached values match non-cached computation
    n <- 100
    k <- 24
    nw <- 12

    # Direct computation
    direct_result <- dpss(n = n, k = k, nw = nw)

    # Cached computation (run twice to test caching)
    cached_result1 <- get_cached_dpss(n, k, nw)
    cached_result2 <- get_cached_dpss(n, k, nw)

    # Results should be identical
    expect_equal(cached_result1$v, direct_result$v, tolerance = 1e-10)
    expect_equal(cached_result2$v, cached_result1$v, tolerance = 1e-10)
})

test_that("DPSS cache improves performance", {
    skip_if_not_installed("multitaper")
    library(multitaper)

    n_values <- c(50, 100, 150, 200, 250)
    k <- 24
    nw <- 12

    # Clear cache first
    if (exists(".dpss_cache", envir = .GlobalEnv)) {
        rm(list = ls(envir = .dpss_cache), envir = .dpss_cache)
    }

    # Time uncached (first call to each length)
    time_uncached <- system.time({
        for (n in n_values) {
            get_cached_dpss(n, k, nw)
        }
    })

    # Time cached (second call to same lengths)
    time_cached <- system.time({
        for (n in n_values) {
            get_cached_dpss(n, k, nw)
        }
    })

    # Cached should be faster
    expect_true(time_cached["elapsed"] < time_uncached["elapsed"])
})

test_that("calc_orf_pval produces valid output structure", {
    # This test requires test data
    skip("Test data loading not implemented yet")
})
