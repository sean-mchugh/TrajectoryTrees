test_that("V2.1 literal history distinguishes all synchronous classes", {
  states <- rbind(
    first_tip = c("B", "B", "B", "B", "C", "A"),
    second_tip = c("B", "B", "C", "B", "C", "A")
  )
  phylogeny_edges <- rbind(
    first_tip = c(1L, 2L, 2L, 2L, 2L, 2L),
    second_tip = c(1L, 3L, 3L, 3L, 3L, 3L)
  )
  intervals <- make_literal_similarity_intervals(states, phylogeny_edges)

  history <- build_pair_history_v2_1(
    intervals,
    "first_tip",
    "second_tip"
  )
  expect_identical(
    history$intervals$sync_class,
    c(
      "Shared", "Conserved", "Divergent",
      "Conserved_Homoplasy", "Parallel", "Parallel"
    )
  )
  expect_identical(
    history$intervals$actual_transition_count_used,
    c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, 1L, 2L)
  )

  similarity <- calculate_pair_similarity_v2_1(
    intervals,
    "first_tip",
    "second_tip",
    async_metric = "bhattacharyya"
  )
  expected <- c(
    B_Shared = 1 / 6,
    B_Conserved = 1 / 6,
    B_Conserved_Homoplasy = 1 / 6,
    C_Parallel_1Transition = 1 / 6,
    A_Parallel_2Transitions = 1 / 6,
    Divergent = 1 / 6
  )
  expect_equal(
    unname(similarity$sync[names(expected)]),
    unname(expected),
    tolerance = 1e-12
  )
  expect_equal(sum(similarity$sync), 1, tolerance = 1e-12)
})

test_that("package similarity equals V2.1 on all 33 reference histories", {
  case_specs <- three_tip_reference_case_specs()
  simmaps <- build_three_tip_reference_simmaps()

  common_v21_surface <- function(result) {
    result$summaries$similarity_vectors <- NULL
    if (!is.null(result$pair_details$interval_classification)) {
      rownames(result$pair_details$interval_classification) <- NULL
    }
    result
  }

  expect_equal(nrow(case_specs), 33L)
  expect_identical(names(simmaps), case_specs$case_id)

  for (case_index in seq_len(nrow(case_specs))) {
    case_id <- case_specs$case_id[[case_index]]
    trajectory <- make_trajectory_objects(
      simmaps[[case_id]],
      root_state = case_specs$root_state[[case_index]],
      include_tt = FALSE,
      include_ss = FALSE,
      include_scenario_mats = TRUE,
      include_path_maps = TRUE,
      time_tolerance = 1e-6
    )

    for (metric in c("bhattacharyya", "minimum")) {
      complete_reference <- trajectory_similarity_v2_1(
        trajectory,
        maximum_transition_counts = 0:3,
        async_metric = metric,
        return_mode = "complete",
        time_tolerance = 1e-6
      )
      complete_candidate <- trajectory_similarity(
        trajectory,
        maximum_transition_counts = 0:3,
        async_metric = metric,
        return_mode = "complete",
        time_tolerance = 1e-6
      )
      compact_reference <- trajectory_similarity_v2_1(
        trajectory,
        maximum_transition_counts = 0:3,
        async_metric = metric,
        return_mode = "matrices_and_summaries",
        time_tolerance = 1e-6
      )
      compact_candidate <- trajectory_similarity(
        trajectory,
        maximum_transition_counts = 0:3,
        async_metric = metric,
        return_mode = "matrices_and_summaries",
        time_tolerance = 1e-6
      )

      expect_equal(
        common_v21_surface(complete_candidate),
        common_v21_surface(complete_reference),
        tolerance = 1e-10,
        info = paste(case_id, metric, "complete")
      )
      expect_equal(
        common_v21_surface(compact_candidate),
        common_v21_surface(compact_reference),
        tolerance = 1e-10,
        info = paste(case_id, metric, "compact")
      )
      expect_equal(
        complete_candidate[c("matrices", "summaries")],
        compact_candidate,
        tolerance = 1e-10,
        info = paste(case_id, metric, "complete versus compact")
      )
    }
  }
})

test_that("reference corpus covers ultrametric and non-ultrametric weighting", {
  case_specs <- three_tip_reference_case_specs()
  simmaps <- build_three_tip_reference_simmaps()

  for (ultrametric in c(TRUE, FALSE)) {
    case_index <- which(case_specs$ultrametric == ultrametric)[[1L]]
    trajectory <- make_trajectory_objects(
      simmaps[[case_specs$case_id[[case_index]]]],
      root_state = case_specs$root_state[[case_index]],
      include_tt = FALSE,
      include_ss = FALSE,
      include_scenario_mats = TRUE,
      include_path_maps = TRUE,
      time_tolerance = 1e-6
    )
    result <- trajectory_similarity(
      trajectory,
      maximum_transition_counts = 0:3,
      return_mode = "complete",
      time_tolerance = 1e-6
    )

    expect_identical(
      result$summaries$weighting$tree_is_ultrametric,
      ultrametric
    )
    expect_identical(
      result$summaries$weighting$available_similarity_weighting_applied,
      !ultrametric
    )
    expect_equal(
      sum(result$summaries$sync$tree_wide_proportions),
      1,
      tolerance = 1e-10
    )
    for (view in result$summaries$async) {
      expect_equal(
        sum(view$tree_wide_proportions),
        1,
        tolerance = 1e-10
      )
    }
  }
})

test_that("package vectors retain the hand-derived tree-wide aggregation", {
  case_specs <- three_tip_reference_case_specs()
  case_id <- "core__root_both_inherit__internal_both_inherit"
  case_index <- match(case_id, case_specs$case_id)
  simmap <- build_three_tip_reference_simmaps()[[case_id]]
  trajectory <- make_trajectory_objects(
    simmap,
    root_state = case_specs$root_state[[case_index]],
    include_tt = FALSE,
    include_ss = FALSE,
    include_scenario_mats = TRUE,
    include_path_maps = TRUE,
    time_tolerance = 1e-6
  )
  vectors <- trajectory_similarity(
    trajectory,
    maximum_transition_counts = 0:3,
    time_tolerance = 1e-6
  )$summaries$similarity_vectors

  expected <- c(
    shared = 2 / 15,
    conserved = 13 / 15,
    conserved_homoplasy = 0,
    convergent = 0,
    parallel_1 = 0,
    parallel_2 = 0,
    parallel_3 = 0,
    divergent = 0
  )
  expect_true(vectors$available)
  expect_equal(
    vectors$sync$all_states$mutually_exclusive,
    expected,
    tolerance = 1e-10
  )
  expect_equal(
    vectors$async_combined$all_states$mutually_exclusive,
    stats::setNames(
      unname(expected),
      sub("^parallel_3$", "parallel_3_plus", names(expected))
    ),
    tolerance = 1e-10
  )
  expect_equal(
    vectors$sync$all_states$cumulative,
    c(
      all_parallel = 0,
      convergent_plus_parallel = 0,
      conserved_homoplasy_plus_convergent_plus_parallel = 0
    ),
    tolerance = 1e-10
  )
})
