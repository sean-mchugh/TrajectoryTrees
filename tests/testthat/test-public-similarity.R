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

test_that("trajectory_similarity supports summary-only output", {
  trajectory <- make_trajectory_objects(
    small_simmap(),
    include_tt = FALSE,
    include_ss = FALSE,
    include_scenario_mats = TRUE,
    include_path_maps = TRUE
  )
  result <- trajectory_similarity(
    trajectory,
    maximum_transition_counts = 0:3,
    async_metric = "bhattacharyya",
    return_mode = "summaries_only"
  )

  expect_named(result, "summaries")
  expect_false("matrices" %in% names(result))
  expect_true(all(c("sync", "async", "components", "weighting") %in%
                    names(result$summaries)))
  expect_false("synchronous_available_similarity_by_pair" %in%
                 names(result$summaries$weighting))
  expect_false("asynchronous_available_similarity_by_pair" %in%
                 names(result$summaries$weighting))
  expect_true(is.list(result$summaries$similarity_vectors))

  reference <- trajectory_similarity(
    trajectory,
    maximum_transition_counts = 0:3,
    async_metric = "bhattacharyya",
    return_mode = "matrices_and_summaries"
  )
  expect_equal(result$summaries$sync, reference$summaries$sync)
  expect_equal(result$summaries$async, reference$summaries$async)
  expect_equal(result$summaries$components, reference$summaries$components)
  expect_equal(result$summaries$similarity_vectors,
               reference$summaries$similarity_vectors)
  expect_equal(
    result$summaries$weighting[c(
      "tree_is_ultrametric", "available_similarity_weighting_applied",
      "tip_history_times"
    )],
    reference$summaries$weighting[c(
      "tree_is_ultrametric", "available_similarity_weighting_applied",
      "tip_history_times"
    )]
  )
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
