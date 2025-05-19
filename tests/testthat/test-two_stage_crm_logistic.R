test_that("two_stage_crm_logistic returns expected structure", {
  res <- two_stage_crm_logistic(seed = 1234)

  expect_type(res, "list")
  expect_named(res, c("n_toxicities", "mtd_estimated", "mle_theta", "x", "y"))
  expect_length(res$x, length(res$y))
  expect_true(all(res$y %in% c(0, 1)))
  expect_gt(res$mtd_estimated, 0)
  expect_true(is.numeric(res$mle_theta))
  expect_length(res$mle_theta, 2)
})
