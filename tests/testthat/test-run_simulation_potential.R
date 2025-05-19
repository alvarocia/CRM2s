test_that("run_simulation_potential returns valid summary data frame", {
  result <- run_simulation_potential(num_rep = 50, seed = 1234)

  # Check the result is a data.frame with 2 rows and expected columns
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_equal(sort(result$method), c("2stage", "3+3"))

  # Check for expected column names
  expected_cols <- c(
    "method", "mean_pat", "median_pat",
    "mean_mtd", "var_mtd", "median_mtd", "min_mtd", "q1_mtd", "q3_mtd", "max_mtd", "iqr_mtd",
    "mean_tox", "var_tox", "median_tox", "min_tox", "q1_tox", "q3_tox", "max_tox", "iqr_tox"
  )
  expect_named(result, expected_cols)

  # Basic checks: all numeric columns should be non-negative
  numeric_cols <- setdiff(expected_cols, "method")
  for (col in numeric_cols) {
    expect_true(all(result[[col]] >= 0), info = paste("Column", col, "should be >= 0"))
  }

})
