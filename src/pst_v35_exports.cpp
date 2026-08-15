#include <Rcpp.h>

#include "pst_v35_batch.hpp"
#include "pst_v35_debug.hpp"
#include "pst_v35_parse.hpp"
#include "pst_v35_ss_output.hpp"
#include "pst_v35_traversal.hpp"
#include "pst_v35_tree_output.hpp"
#include "pst_v35_tt_output.hpp"

#include <chrono>
#include <string>

/**
 * Round-trip one compact input for bridge and parser diagnostics.
 *
 * Input ownership remains with Rcpp for this call. Output copies every
 * edge/map vector plus structural counts; parse failure throws before a partial
 * diagnostic object is created.
 */
// [[Rcpp::export]]
Rcpp::List pst_v35_cpp_roundtrip(Rcpp::List payload) {
  V35Input input = pst_v35_parse_input(payload);
  return Rcpp::List::create(
    Rcpp::Named("summary") = pst_v35_traversal_smoke_summary(input),
    Rcpp::Named("edge_parent") = input.edge_parent,
    Rcpp::Named("edge_child") = input.edge_child,
    Rcpp::Named("edge_length") = input.edge_length,
    Rcpp::Named("map_edge_id") = input.map_edge_id,
    Rcpp::Named("map_state_id") = input.map_state_id,
    Rcpp::Named("map_duration") = input.map_duration,
    Rcpp::Named("map_edge_offset") = input.map_edge_offset,
    Rcpp::Named("tip_label") = input.tip_label
  );
}

/**
 * Convert public numeric times through the exact integer-key implementation
 * used by V35 event batching.
 */
// [[Rcpp::export]]
Rcpp::NumericVector pst_v35_cpp_batch_time_keys(
    Rcpp::NumericVector times,
    double time_tolerance) {
  double time_scale = pst_v35_batch_time_scale(time_tolerance);
  Rcpp::NumericVector keys(times.size());
  for (int index = 0; index < times.size(); ++index) {
    if (!R_finite(times[index])) {
      Rcpp::stop("V35 batch-time key input contains a nonfinite time");
    }
    keys[index] = static_cast<double>(pst_v35_batch_time_key(
      times[index],
      time_scale
    ));
  }
  return keys;
}

/**
 * Parse a simmap directly and execute V35 without R-side map flattening.
 *
 * Parse time is recorded separately from traversal stages. Biological ids stay
 * 1-based across the boundary and the resolved root policy is immutable.
 */
// [[Rcpp::export]]
Rcpp::List pst_v35_cpp_active_summary_from_tree(
    Rcpp::List tree,
    Rcpp::DataFrame state_lookup,
    Rcpp::IntegerVector root_anchor_state_ids,
    bool root_is_synthetic,
    bool recover_debug_failure = false,
    int debug_fail_batch_id = NA_INTEGER,
    double time_tolerance = 1e-6) {
  typedef std::chrono::steady_clock V35Clock;
  V35Clock::time_point parse_start = V35Clock::now();
  V35Input input = pst_v35_parse_tree_input(
    tree,
    state_lookup,
    root_anchor_state_ids,
    root_is_synthetic
  );
  V35Clock::time_point parse_end = V35Clock::now();
  Rcpp::List out = pst_v35_extract_active_summary(
    input,
    recover_debug_failure,
    debug_fail_batch_id,
    time_tolerance
  );
  Rcpp::List metadata = out["metadata"];
  Rcpp::List stage_timings = metadata["stage_timings"];
  stage_timings["parse_schedule"] = std::chrono::duration<double>(
    parse_end - parse_start
  ).count();
  metadata["stage_timings"] = stage_timings;
  out["metadata"] = metadata;
  return out;
}

/**
 * Clip one mapped PST projection at an atomically committed debug time.
 *
 * `tree` and `cutoff_time` are read-only R inputs. The returned simmap preserves
 * map/size-map alignment through the cutoff and uses deterministic provisional
 * tips for active frontiers. The debug module owns naming and node remapping;
 * the final tree is never mutated. Invalid time or malformed tree fields throw
 * before a partial snapshot crosses the boundary.
 */
// [[Rcpp::export]]
Rcpp::List pst_v35_cpp_debug_clip_tree(
    Rcpp::List tree,
    double cutoff_time,
    double time_tolerance) {
  return pst_v35_debug_clip_tree_cpp(tree, cutoff_time, time_tolerance);
}

/**
 * Format one online event-record family on the shared committed TT axis.
 *
 * `canonical_trans` is the traversal-owned numerical source. `records` only
 * preserve the ordered legacy event vocabulary, `summary` supplies labels,
 * and `time_vec` is the committed batch grid. This boundary performs no
 * scoring and never mutates inputs; inconsistent identities raise an R error.
 */
// [[Rcpp::export]]
Rcpp::List pst_v35_cpp_tt_transition_views(
    Rcpp::DataFrame records,
    Rcpp::NumericMatrix canonical_trans,
    Rcpp::List summary,
    Rcpp::NumericVector time_vec,
    double time_tolerance) {
  return pst_v35_transition_views_cpp(
    records,
    canonical_trans,
    summary,
    time_vec,
    time_tolerance
  );
}

/**
 * Derive dependent TT total/state views from one canonical path matrix.
 *
 * Inputs are read only. The path matrix remains the sole numeric source;
 * lookups provide semantic aggregation labels and output order.
 */
// [[Rcpp::export]]
Rcpp::List pst_v35_cpp_tt_path_counter_views(
    Rcpp::NumericMatrix path_matrix,
    Rcpp::DataFrame path_lookup,
    Rcpp::DataFrame state_lookup) {
  return pst_v35_path_counter_views_cpp(
    path_matrix,
    path_lookup,
    state_lookup
  );
}

/**
 * Assemble SS from compact P/S/T metrics and traversal-owned TT/path support.
 *
 * `tree_parts` contains materialization-owned compatibility metrics,
 * `path_size_maps` contains online through/stopped support, and optional `tt`
 * contains online cumulative counters. The return matches the V30 SS field
 * names and shapes without rescoring events or scanning public trees. Inputs
 * are read only; incomplete required parts fail rather than being reconstructed.
 */
// [[Rcpp::export]]
Rcpp::List pst_v35_cpp_summary_stats_from_parts(
    Rcpp::List tree_parts,
    SEXP path_size_maps,
    Rcpp::Nullable<Rcpp::List> tt = R_NilValue) {
  return pst_v35_build_summary_stats_from_parts_cpp(
    tree_parts,
    path_size_maps,
    tt
  );
}

/**
 * Materialize the initial scenario tree from traversal-owned segment records.
 *
 * `summary` is read only and must contain aligned state, path, size, duration,
 * stable-parent, and tip-membership records. The returned simmap has exactly
 * one edge per stable traversal scenario.
 */
// [[Rcpp::export]]
Rcpp::List pst_v35_cpp_scenario_tree_initial(Rcpp::List summary) {
  return pst_v35_scenario_tree_initial_impl(summary);
}

/**
 * Validate a traversal-owned scenario tree for the public PST contract.
 *
 * The return preserves every topology edge and segment boundary while exposing
 * state-only maps, path maps, lineage support, exact edge types, and terminal
 * events. Inputs are not mutated; invalid stable topology raises an R error.
 */
// [[Rcpp::export]]
Rcpp::List pst_v35_cpp_canonicalize_scenario_tree(
    Rcpp::List tree,
    Rcpp::List summary,
    double time_tolerance) {
  return pst_v35_canonicalize_scenario_tree_impl(
    tree,
    summary,
    time_tolerance
  );
}

/**
 * Preserve traversal scenario IDs after direct tree serialization.
 *
 * `summary` is cloned and `scenario_tree` is read only. Formatting has no
 * authority to rewrite scenario matrices, terminal records, or LDIF.
 */
// [[Rcpp::export]]
Rcpp::List pst_v35_cpp_canonicalize_summary_scenario_ids(
    Rcpp::List summary,
    Rcpp::List scenario_tree) {
  return pst_v35_canonicalize_summary_scenario_ids_impl(
    summary,
    scenario_tree
  );
}

/**
 * Build finalized edge/step coordinates for the complete scenario matrix.
 *
 * All inputs are read only. The returned list always contains finalized
 * scenario-edge ids and both step-id matrices. When `include_paths` is true it
 * also contains the scenario and phylogeny path labels resolved by those exact
 * coordinates.
 */
// [[Rcpp::export]]
Rcpp::List pst_v35_cpp_scenario_matrix_coordinates(
    Rcpp::List scenario_tree,
    Rcpp::List phylo_tree,
    Rcpp::NumericVector time_vec,
    Rcpp::IntegerMatrix existing_scenario_ids,
    Rcpp::IntegerMatrix existing_phylo_ids,
    Rcpp::IntegerMatrix existing_path_ids,
    Rcpp::DataFrame path_lookup,
    bool include_paths,
    double time_tolerance) {
  return pst_v35_scenario_matrix_coordinates_impl(
    scenario_tree,
    phylo_tree,
    time_vec,
    existing_scenario_ids,
    existing_phylo_ids,
    existing_path_ids,
    path_lookup,
    include_paths,
    time_tolerance
  );
}

/**
 * Build the public transition tree from committed transition-history support.
 *
 * `transition_support` contains online state/edge histories and the shared time
 * grid; `tip_labels` supplies biological identities. The returned simmap groups
 * equal histories, splits at first divergence, and carries deterministic
 * topology/maps. Inputs remain unchanged; malformed support fails before
 * serialization.
 */
// [[Rcpp::export]]
Rcpp::List pst_v35_cpp_transition_tree_from_support(
    Rcpp::List transition_support,
    Rcpp::CharacterVector tip_labels,
    double time_tolerance) {
  return pst_v35_transition_tree_from_support_impl(
    transition_support,
    tip_labels,
    time_tolerance
  );
}

/**
 * Build a state-only public simmap from a V35 state/size-labelled tree.
 *
 * `tree` is cloned in C++ before relabelling. The return strips V35 size
 * suffixes, merges adjacent equal biological states, and rebuilds mapped-edge
 * exposure under biological state columns. Optional summary arguments attach
 * compact SS metric parts during native materialization. Inputs remain
 * unchanged; malformed maps or topology raise an R error before any partial
 * tree is returned.
 */
// [[Rcpp::export]]
Rcpp::List pst_v35_cpp_state_only_tree(
    Rcpp::List tree,
    std::string summary_suffix = "",
    Rcpp::Nullable<Rcpp::CharacterVector> state_levels = R_NilValue,
    bool include_edge_boundaries = false,
    bool source_labels_are_state_only = false) {
  return pst_v35_state_only_tree_impl(
    tree,
    summary_suffix,
    state_levels,
    include_edge_boundaries,
    source_labels_are_state_only
  );
}

/**
 * Relabel one state tree with canonical root-to-state path identities.
 *
 * `tree` is cloned, `root_policy` seeds the initial path, and `sep` is the
 * escaped public path delimiter. The return includes path-labelled maps,
 * mapped-edge exposure, path ids/components, and a canonical lookup. Parent
 * paths must exist before daughters in root-first order; broken topology fails
 * explicitly.
 */
// [[Rcpp::export]]
Rcpp::List pst_v35_cpp_pathify_state_tree(
    Rcpp::List tree,
    Rcpp::List root_policy,
    std::string sep = "|") {
  return pst_v35_pathify_state_tree_impl(tree, root_policy, sep);
}

/**
 * Materialize transition-tree size maps from traversal-owned online counters.
 *
 * `state_tree` and `path_tree` are parallel native tree projections;
 * `cumulative_lineage_trans` and `cumulative_scenario_trans` own internal
 * through sizes; traversal-committed terminal lineage and scenario ledgers own
 * tip sizes for each exact canonical path. The return contains four
 * edge-map families and an ownership marker. Inputs are read only; malformed
 * topology, labels, or matrix dimensions fail before any public tree is modified.
 */
// [[Rcpp::export]]
Rcpp::List pst_v35_cpp_transition_size_maps(
    Rcpp::List state_tree,
    Rcpp::List path_tree,
    Rcpp::NumericMatrix cumulative_lineage_trans,
    Rcpp::NumericMatrix cumulative_scenario_trans,
    Rcpp::CharacterVector path_labels_by_id,
    Rcpp::IntegerVector terminal_lineage_counts_by_id,
    Rcpp::IntegerVector terminal_scenario_counts_by_id,
    Rcpp::IntegerVector terminal_trajectory_group_id_by_tip) {
  return pst_v35_transition_size_maps_impl(
    state_tree,
    path_tree,
    cumulative_lineage_trans,
    cumulative_scenario_trans,
    path_labels_by_id,
    terminal_lineage_counts_by_id,
    terminal_scenario_counts_by_id,
    terminal_trajectory_group_id_by_tip
  );
}
