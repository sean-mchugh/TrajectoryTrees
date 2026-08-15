source(
  testthat::test_path("reference", "ScenarioMatsV2_1.R"),
  local = environment()
)
source(
  testthat::test_path(
    "reference",
    "construct_three_tip_reference_trees.R"
  ),
  local = environment()
)

make_literal_similarity_intervals <- function(
    states,
    phylogeny_edge_ids,
    root_anchor = "B",
    duration = rep(1, ncol(states)),
    maximum_transition_counts = 0:3) {
  tip_labels <- rownames(states)
  column_count <- ncol(states)
  list(
    tip_labels = tip_labels,
    time = data.frame(
      interval_id = seq_len(column_count),
      source_column = seq_len(column_count),
      start = c(0, head(cumsum(duration), -1L)),
      end = cumsum(duration),
      duration = duration
    ),
    states = states,
    paths = matrix(
      paste0("|", as.character(states)),
      nrow(states),
      ncol(states),
      dimnames = dimnames(states)
    ),
    phylo_edge_ids = phylogeny_edge_ids,
    scenario_edge_ids = matrix(
      1L, nrow(states), column_count,
      dimnames = dimnames(states)
    ),
    scenario_edge_step_ids = matrix(
      1L, nrow(states), column_count,
      dimnames = dimnames(states)
    ),
    active = matrix(
      TRUE, nrow(states), column_count,
      dimnames = dimnames(states)
    ),
    duration = duration,
    root_anchor = root_anchor,
    root_policy = list(
      anchor = root_anchor,
      provenance = "explicit",
      synthetic = FALSE
    ),
    maximum_transition_counts = as.integer(maximum_transition_counts),
    state_levels = sort(unique(as.character(states))),
    time_tolerance = 1e-10
  )
}
