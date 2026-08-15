expect_summary_reconstructed_from_matrices <- function(
    matrix_family,
    summary,
    available_similarity_by_pair = NULL,
    tolerance = 1e-10) {
  expect_true(length(matrix_family) > 0L)
  component_names <- names(matrix_family)
  upper_triangle <- upper.tri(matrix_family[[1L]])

  expected_totals <- vapply(
    matrix_family,
    function(component_matrix) {
      values <- component_matrix[upper_triangle]
      if (!is.null(available_similarity_by_pair)) {
        values <- values * available_similarity_by_pair[upper_triangle]
      }
      sum(values)
    },
    numeric(1)
  )
  names(expected_totals) <- component_names
  expected_total <- sum(expected_totals)
  expected_proportions <- expected_totals / expected_total

  expect_equal(summary$totals, expected_totals, tolerance = tolerance)
  expect_equal(
    summary$tree_wide_proportions,
    expected_proportions,
    tolerance = tolerance
  )
  expect_equal(summary$means, expected_proportions, tolerance = tolerance)
  expect_equal(
    summary$total_across_matrices,
    expected_total,
    tolerance = tolerance
  )
  if (!is.null(available_similarity_by_pair)) {
    expect_equal(
      summary$available_similarity_total,
      sum(available_similarity_by_pair[upper_triangle]),
      tolerance = tolerance
    )
  }
}

test_that("hand-derived three-tip matrices have the correct pair values", {
  case_id <- "core__root_both_inherit__internal_both_inherit"
  case_specs <- three_tip_reference_case_specs()
  simmap <- build_three_tip_reference_simmaps()[[case_id]]
  trajectory <- make_trajectory_objects(
    simmap,
    root_state = case_specs$root_state[match(case_id, case_specs$case_id)],
    include_tt = FALSE,
    include_ss = FALSE,
    include_scenario_mats = TRUE,
    include_path_maps = TRUE,
    time_tolerance = 1e-6
  )
  result <- trajectory_similarity(
    trajectory,
    maximum_transition_counts = 0:3,
    async_metric = "bhattacharyya",
    return_mode = "complete",
    time_tolerance = 1e-6
  )

  expected_shared <- matrix(
    c(0, 0, 0, 0, 0, 0.4, 0, 0.4, 0),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(trajectory$phylo$tip.label, trajectory$phylo$tip.label)
  )
  expected_conserved <- matrix(
    c(0, 1, 1, 1, 0, 0.6, 1, 0.6, 0),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(trajectory$phylo$tip.label, trajectory$phylo$tip.label)
  )

  expect_equal(result$matrices$sync$A_Shared, expected_shared, tolerance = 1e-10)
  expect_equal(
    result$matrices$sync$A_Conserved,
    expected_conserved,
    tolerance = 1e-10
  )
  remaining_components <- setdiff(
    names(result$matrices$sync),
    c("A_Shared", "A_Conserved")
  )
  expect_true(all(vapply(
    result$matrices$sync[remaining_components],
    function(component_matrix) all(abs(component_matrix) <= 1e-10),
    logical(1)
  )))
})

test_that("reported summaries are independently recovered from matrices", {
  case_specs <- three_tip_reference_case_specs()
  case_ids <- c(
    "paths__unit_cycles_both_both",
    "nonultra__unit_cycles_both_both__short_tip_1_by_1"
  )
  simmaps <- build_three_tip_reference_simmaps()

  for (case_id in case_ids) {
    trajectory <- make_trajectory_objects(
      simmaps[[case_id]],
      root_state = case_specs$root_state[match(case_id, case_specs$case_id)],
      include_tt = FALSE,
      include_ss = FALSE,
      include_scenario_mats = TRUE,
      include_path_maps = TRUE,
      time_tolerance = 1e-6
    )
    for (metric in c("bhattacharyya", "minimum")) {
      result <- trajectory_similarity(
        trajectory,
        maximum_transition_counts = 0:3,
        async_metric = metric,
        return_mode = "complete",
        time_tolerance = 1e-6
      )
      weighted <- isTRUE(
        result$summaries$weighting$available_similarity_weighting_applied
      )
      sync_weights <- if (weighted) {
        result$summaries$weighting$synchronous_available_similarity_by_pair
      } else {
        NULL
      }
      async_weights <- if (weighted) {
        result$summaries$weighting$asynchronous_available_similarity_by_pair
      } else {
        NULL
      }

      expect_summary_reconstructed_from_matrices(
        result$matrices$sync,
        result$summaries$sync,
        sync_weights
      )
      for (view_name in names(result$matrices$async)) {
        expect_summary_reconstructed_from_matrices(
          result$matrices$async[[view_name]],
          result$summaries$async[[view_name]],
          async_weights
        )
      }
    }
  }
})
