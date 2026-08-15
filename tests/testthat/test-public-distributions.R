test_that("trajectory ensembles produce a stable distribution object", {
  constructor <- public_function("make_trajectory_objects")
  reducer <- public_function("make_trajectory_distributions")

  expect_true(is.function(constructor))
  expect_true(is.function(reducer))
  if (!is.function(constructor) || !is.function(reducer)) {
    return(invisible())
  }

  tree <- small_simmap()
  samples <- constructor(list(tree, tree), include_tt = FALSE)
  result <- reducer(samples, include = c("P", "T", "SS"), backend = "R")

  expect_s3_class(result, "trajectory_distribution")
  expect_named(result, c("metadata", "P", "T", "SS", "TT"))
  expect_equal(result$metadata$n_samples, 2L)
})
