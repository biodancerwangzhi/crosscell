# Extracted from test-error-handling.R:22

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "crosscell", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
library(testthat)
library(crosscell)

# test -------------------------------------------------------------------------
tmp <- tempfile(fileext = ".h5ad")
file.create(tmp)
on.exit(unlink(tmp))
expect_error(
    read_rds_fast(tmp),
    "Expected .rds"
  )
