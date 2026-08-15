#ifndef PST_V34_TREE_OUTPUT_HPP
#define PST_V34_TREE_OUTPUT_HPP

#include "pst_v34_types.hpp"

#include <string>

/**
 * Materialize one public topology edge per traversal-owned stable scenario ID.
 */
Rcpp::List pst_v34_scenario_tree_initial_impl(Rcpp::List summary);

/** Validate and expose traversal-owned topology without semantic repair. */
Rcpp::List pst_v34_canonicalize_scenario_tree_impl(
  Rcpp::List tree,
  Rcpp::List summary,
  double time_tolerance);

/** Return an unchanged summary after validating direct stable-ID provenance. */
Rcpp::List pst_v34_canonicalize_summary_scenario_ids_impl(
  Rcpp::List summary,
  Rcpp::List scenario_tree);

/**
 * Rebuild only the scenario-edge id matrix against the finalized scenario tree.
 * All inputs remain read only; the return is a new integer matrix with the same
 * dimensions as `existing_ids` and finalized `scenario_tree$edge` row ids.
 */
Rcpp::IntegerMatrix pst_v34_final_scenario_edge_ids_impl(
  Rcpp::List scenario_tree,
  Rcpp::CharacterVector phylo_tip_labels,
  Rcpp::NumericVector time_vec,
  Rcpp::IntegerMatrix existing_ids,
  double time_tolerance);

/**
 * Build exact tree-map coordinates for every populated scenario-matrix cell.
 * The result contains finalized scenario-edge ids plus scenario and phylogeny
 * step ids; optional character matrices expose the resolved path labels.
 */
Rcpp::List pst_v34_scenario_matrix_coordinates_impl(
  Rcpp::List scenario_tree,
  Rcpp::List phylo_tree,
  Rcpp::NumericVector time_vec,
  Rcpp::IntegerMatrix existing_scenario_ids,
  Rcpp::IntegerMatrix existing_phylo_ids,
  bool include_paths,
  double time_tolerance);

/** Build the transition tree from C++ transition-history support. */
Rcpp::List pst_v34_transition_tree_from_support_impl(
  Rcpp::List transition_support,
  Rcpp::CharacterVector tip_labels,
  double time_tolerance);

/** Apply terminal STOPPED-map correction without changing biological paths. */
Rcpp::List pst_v34_correct_stopped_branches_impl(
  Rcpp::List trans_tree,
  std::string stop_char = "STOPPED",
  bool shrink_bl = true);

/** Strip V34 size/path suffixes and rebuild state-only mapped-edge exposure. */
Rcpp::List pst_v34_state_only_tree_impl(
  Rcpp::List tree,
  std::string summary_suffix = "",
  Rcpp::Nullable<Rcpp::CharacterVector> state_levels = R_NilValue,
  bool include_edge_boundaries = false,
  bool source_labels_are_state_only = false);

/** Relabel a state tree with canonical root-to-state path identities. */
Rcpp::List pst_v34_pathify_state_tree_impl(
  Rcpp::List tree,
  Rcpp::List root_policy,
  std::string sep = "|");

/** Materialize transition-tree size maps from committed online path support. */
Rcpp::List pst_v34_transition_size_maps_impl(
  Rcpp::List state_tree,
  Rcpp::List path_tree,
  Rcpp::NumericMatrix cumulative_lineage_trans,
  Rcpp::NumericMatrix cumulative_scenario_trans,
  Rcpp::CharacterVector path_labels_by_id,
  Rcpp::IntegerVector terminal_lineage_counts_by_id,
  Rcpp::IntegerVector terminal_scenario_counts_by_id,
  Rcpp::IntegerVector terminal_trajectory_group_id_by_tip);

#endif
