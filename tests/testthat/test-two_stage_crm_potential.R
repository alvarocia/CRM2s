test_that("two_stage_crm_potential runs and returns expected output", {
  result <- two_stage_crm_potential(seed = 1234)

  expect_type(result, "list")
  expect_named(result, c("n_toxicities", "mtd_estimated", "mle_theta", "x", "y"))

  expect_length(result$x, length(result$y))
  expect_true(all(result$y %in% c(0, 1)))
  expect_true(result$mtd_estimated > 0)
  expect_true(result$mle_theta > 0)
})
