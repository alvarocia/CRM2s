test_that("logistic_3_3 returns correct structure and values", {
  res <- logistic_3_3(seed = 1234)

  expect_type(res, "list")
  expect_named(res, c("n_toxicities", "mtd_estimated", "n_patients", "x", "y"))
  expect_length(res$x, res$n_patients)
  expect_length(res$y, res$n_patients)
  expect_true(all(res$y %in% c(0, 1)))
  expect_gt(res$mtd_estimated, 0)
})
