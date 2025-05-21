test_that("power_3_3 returns correct structure and values", {
  result <- power_3_3(seed = 1234)

  expect_type(result, "list")
  expect_named(result, c("n_toxicities", "mtd_estimated", "n_patients", "x", "y"))

  expect_length(result$x, result$n_patients)
  expect_length(result$y, result$n_patients)
  expect_true(all(result$y %in% c(0, 1)))
  expect_true(result$mtd_estimated > 0)
})
