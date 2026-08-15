test_that("trajectory_similarity is a public package-native calculation", {
  similarity_function <- public_function("trajectory_similarity")

  expect_true(is.function(similarity_function))
  if (!is.function(similarity_function)) {
    return(invisible())
  }

  trajectory <- make_trajectory_objects(
    small_simmap(),
    include_tt = FALSE,
    include_ss = FALSE,
    include_scenario_mats = TRUE,
    include_path_maps = TRUE
  )
  result <- similarity_function(
    trajectory,
    maximum_transition_counts = 0:3,
    async_metric = "bhattacharyya",
    return_mode = "matrices_and_summaries"
  )

  expect_named(result, c("matrices", "summaries"))
  expect_true(is.loaded("_TrajectoryTrees_scenario_mats_v3_calculate_cpp"))
  expect_true(all(c(
    "sync", "async", "similarity_vectors"
  ) %in% names(result$summaries)))
})

test_that("trajectory_similarity validates its public arguments", {
  similarity_function <- public_function("trajectory_similarity")

  expect_true(is.function(similarity_function))
  if (!is.function(similarity_function)) {
    return(invisible())
  }

  trajectory <- make_trajectory_objects(
    small_simmap(),
    include_tt = FALSE,
    include_ss = FALSE
  )

  expect_error(
    similarity_function(trajectory, maximum_transition_counts = 1:3),
    "include 0"
  )
  expect_error(
    similarity_function(trajectory, async_metric = "unsupported"),
    "bhattacharyya"
  )
  expect_error(
    similarity_function(trajectory, return_mode = "unsupported"),
    "matrices_and_summaries"
  )
})
