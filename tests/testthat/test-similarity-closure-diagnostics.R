test_that("all-state closure errors expose the failing vector", {
  metadata <- data.frame(
    matrix_name = c(
      "A_Shared", "Divergent",
      "A_Shared", "Divergent"
    ),
    state = c("A", NA, "A", NA),
    phase = c("Shared", "Independent", "Shared", "Independent"),
    kind = c("state", "residual", "state", "residual"),
    requested_maximum_transition_count = c(NA, NA, 0L, 0L),
    actual_transition_count_used = c(NA, NA, 0L, 0L),
    view = c("sync", "sync", "state_only", "state_only"),
    stringsAsFactors = FALSE
  )
  result <- list(summaries = list(
    components = metadata,
    sync = list(tree_wide_proportions = c(A_Shared = 0.4, Divergent = 0.5)),
    async = list(state_only = list(
      tree_wide_proportions = c(A_Shared = 0.4, Divergent = 0.6)
    ))
  ))

  error <- tryCatch(
    {
      TrajectoryTrees:::.scenario_mats_v3_add_similarity_vectors(result)
      NULL
    },
    error = function(error) error
  )

  expect_s3_class(error, "trajectorytrees_similarity_closure_error")
  expect_identical(error$view, "sync")
  expect_equal(error$total, 0.9, tolerance = 1e-12)
  expect_equal(error$tolerance, 1e-8, tolerance = 0)
  expect_equal(
    error$values,
    c(
      shared = 0.4,
      conserved = 0,
      conserved_homoplasy = 0,
      convergent = 0,
      divergent = 0.5
    ),
    tolerance = 1e-12
  )
})
