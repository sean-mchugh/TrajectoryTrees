# V35 C++ bridge.
#
# Role:
#   Expose small R helpers around the native routines registered when the
#   TrajectoryTrees package is loaded.
#
# Shape note:
#   The round-trip result is an internal debug list with scalar counts and
#   vectors copied back from the C++ boundary. Debug artifact example:
#   roundtrip <- readRDS("testing/correctness/reference/v35_intermediates/five_tip_cpp_example/cpp_roundtrip.rds")

#' Verify that the package's V35 native routines are loaded.
#'
#' Native compilation happens once during package installation. Normal package
#' use never invokes a compiler and never calls `Rcpp::sourceCpp()`.
#' @return Invisibly returns `TRUE`; throws if the package DLL is unavailable.
#' @noRd
pst_v35_load_cpp <- function() {
  if (!is.loaded("_TrajectoryTrees_pst_v35_cpp_roundtrip")) {
    stop(
      "TrajectoryTrees native routines are not loaded; reinstall the package",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Send a compact encoded input payload through the C++ smoke-test boundary.
#'
#' A C++ round-trip diagnostic list.
#' Loads the V35 C++ translation unit before crossing the R/C++ boundary.
#' @param payload R list payload crossing the V35 R/C++ boundary.
#' @return A C++ round-trip diagnostic list.
#' @noRd
pst_v35_cpp_roundtrip_payload <- function(payload) {
  # Load the C++ symbol immediately before the first bridge call.
  pst_v35_load_cpp()

  # Call the C++ smoke kernel once with the compact payload.
  pst_v35_cpp_roundtrip(payload)
}

#' Convert numeric event times to native V35 batch keys.
#'
#' @param times Numeric vector of root-relative event times.
#' @param time_tolerance Positive V35 event-key resolution.
#' @return Numeric vector containing exact integer-valued native keys.
#' @noRd
pst_v35_cpp_batch_time_keys_payload <- function(times, time_tolerance) {
  pst_v35_load_cpp()
  pst_v35_cpp_batch_time_keys(
    as.numeric(times),
    as.numeric(time_tolerance)
  )
}

#' Build active traversal summary matrices directly from a simmap tree through C++.
#'
#' A V35 active-summary list aligned to the supplied state lookup and root anchor ids.
#' Loads the V35 C++ translation unit before crossing the R/C++ boundary.
#' @param tree Simmap/phylo tree being formatted, relabelled, or canonicalized.
#' @param state_lookup Data frame mapping V35 integer state ids to state labels.
#' @param root_anchor_state_ids Integer state id or ids used as the synthetic root path anchor.
#' @param root_is_synthetic Logical flag distinguishing a bookkeeping root from a biological root state.
#' @param recover_debug_failure Whether a failed batch returns its pre-batch checkpoint.
#' @param debug_fail_batch_id Internal deterministic failure injection used by the debug contract test.
#' @param time_tolerance Positive numeric resolution used for both V35 event
#'   batches and public scenario-tree cleanup.
#' @return A V35 active-summary list aligned to the supplied state lookup and root anchor ids.
#' @noRd
pst_v35_cpp_active_summary_tree_payload <- function(
    tree,
    state_lookup,
    root_anchor_state_ids,
    root_is_synthetic,
    recover_debug_failure = FALSE,
    debug_fail_batch_id = NA_integer_,
    time_tolerance = 1e-6) {
  # Load the C++ symbol immediately before active traversal extraction.
  pst_v35_load_cpp()

  # Return active traversal matrices directly from a simmap tree.
  pst_v35_cpp_active_summary_from_tree(
    tree,
    state_lookup,
    root_anchor_state_ids,
    isTRUE(root_is_synthetic),
    isTRUE(recover_debug_failure),
    as.integer(debug_fail_batch_id),
    as.numeric(time_tolerance)
  )
}

#' Clip one V35 mapped tree at a committed debug time in C++.
#'
#' @param tree Public V35 phylogeny, scenario tree, or transition tree.
#' @param cutoff_time Absolute root-relative committed time.
#' @param time_tolerance Shared V35 resolution for clipping comparisons.
#' @return Clipped mapped tree with deterministic provisional frontier tips.
#' @noRd
pst_v35_cpp_debug_clip_tree_payload <- function(
    tree,
    cutoff_time,
    time_tolerance) {
  pst_v35_load_cpp()
  pst_v35_cpp_debug_clip_tree(tree, cutoff_time, time_tolerance)
}

#' Build dependent transition-event views from canonical live TT rows.
#'
#' A list of event, arrival, leaving, trans, state, and total TT matrices.
#' Loads the V35 C++ translation unit before crossing the R/C++ boundary.
#' @param records Data frame or C++ record table of transition/path events.
#' @param canonical_trans Traversal-owned canonical path matrix.
#' @param summary Canonical V35 traversal summary with lookup tables, scenario matrices, and edge records.
#' @param time_vec Numeric V35 shared time axis.
#' @param time_tolerance Shared V35 resolution for time-axis comparisons.
#' @return A list of event, arrival, leaving, path, and total TT matrices.
#' @noRd
pst_v35_cpp_tt_transition_views_payload <- function(
    records,
    canonical_trans,
    summary,
    time_vec,
    time_tolerance) {
  # Load the C++ symbol immediately before native TT event-view construction.
  pst_v35_load_cpp()

  # Derive every numerical view from traversal-owned canonical path rows.
  pst_v35_cpp_tt_transition_views(
    records,
    canonical_trans,
    summary,
    time_vec,
    time_tolerance
  )
}

#' Derive public TT total/state views from one canonical path matrix.
#'
#' @param path_matrix Traversal-owned cumulative matrix keyed by canonical path.
#' @param path_lookup Canonical path ids, labels, parents, and terminal states.
#' @param state_lookup Canonical state ids and labels defining output order.
#' @return A list containing dependent `tot`, `state`, and unchanged `trans` matrices.
#' @noRd
pst_v35_cpp_tt_path_counter_views_payload <- function(
    path_matrix,
    path_lookup,
    state_lookup) {
  pst_v35_load_cpp()
  pst_v35_cpp_tt_path_counter_views(path_matrix, path_lookup, state_lookup)
}

#' Format the initial scenario tree from a V35 summary through C++.
#'
#' A labelled simmap scenario tree before canonical cleanup passes.
#' Loads the V35 C++ translation unit before crossing the R/C++ boundary.
#' @param summary Canonical V35 traversal summary with lookup tables, scenario matrices, and edge records.
#' @return A labelled simmap scenario tree before canonical cleanup passes.
#' @noRd
pst_v35_cpp_scenario_tree_initial_payload <- function(summary) {
  # Load the C++ symbol immediately before native scenario-tree formatting.
  pst_v35_load_cpp()

  # Return the initial labelled scenario simmap before cleanup canonicalization.
  pst_v35_cpp_scenario_tree_initial(summary)
}

#' Run C++ scenario-tree canonicalization and compaction passes.
#'
#' A canonical labelled simmap scenario tree.
#' Loads the V35 C++ translation unit before crossing the R/C++ boundary.
#' @param tree Simmap/phylo tree being formatted, relabelled, or canonicalized.
#' @param summary Canonical V35 traversal summary with lookup tables, scenario matrices, and edge records.
#' @param time_tolerance Positive numeric resolution used for scenario-tree cleanup.
#' @return A canonical labelled simmap scenario tree.
#' @noRd
pst_v35_cpp_canonicalize_scenario_tree_payload <- function(
    tree,
    summary,
    time_tolerance) {
  # Load the C++ symbol immediately before scenario-tree cleanup.
  pst_v35_load_cpp()

  # Return the canonical scenario simmap after C++ cleanup passes.
  pst_v35_cpp_canonicalize_scenario_tree(tree, summary, time_tolerance)
}

#' Rewrite summary scenario ids to match merged public scenario-tree tips.
#'
#' A V35 summary whose scenario-id matrix uses canonical public scenario identities.
#' Loads the V35 C++ translation unit before crossing the R/C++ boundary.
#' @param summary Canonical V35 traversal summary with lookup tables, scenario matrices, and edge records.
#' @param scenario_tree Optional scenario tree used to align TT or size-support identities to public tree ids.
#' @return A V35 summary whose scenario-id matrix uses canonical public scenario identities.
#' @noRd
pst_v35_cpp_canonicalize_summary_scenario_ids_payload <- function(summary, scenario_tree) {
  # Load the C++ symbol immediately before scenario-id cleanup.
  pst_v35_load_cpp()

  # Return the summary with scenario ids canonicalized to merged public scenario
  # tips. This keeps STT and scenario_mats in the same scenario identity space
  # as the canonical public scenario tree.
  pst_v35_cpp_canonicalize_summary_scenario_ids(summary, scenario_tree)
}

#' Build exact scenario-matrix edge/step coordinates through C++.
#'
#' The native implementation rebuilds scenario-edge ids in the finalized tree
#' namespace and adds one step id for each scenario and phylogeny edge id. When
#' requested, it also resolves the character path label at both coordinates.
#' Every input is read only; the return is a list of new matrices sharing the
#' dimensions and dimnames of `existing_scenario_ids`.
#' @param scenario_tree Final state tree with optional edge-parallel path maps.
#' @param phylo_tree Input state tree with optional edge-parallel path maps.
#' @param time_vec Committed scenario-matrix column times.
#' @param existing_scenario_ids Traversal scenario ids used for layout only.
#' @param existing_phylo_ids Final phylogeny-edge ids for every matrix cell.
#' @param include_paths Whether to return resolved scenario/phylogeny paths.
#' @param time_tolerance Positive V35 integer-time-key resolution.
#' @return List containing finalized edge ids, both step-id matrices, and
#'   optional path matrices.
#' @noRd
pst_v35_cpp_scenario_matrix_coordinates_payload <- function(
    scenario_tree,
    phylo_tree,
    time_vec,
    existing_scenario_ids,
    existing_phylo_ids,
    existing_path_ids,
    path_lookup,
    include_paths,
    time_tolerance) {
  # The coordinate builder lives in the established V35 translation unit, so
  # ordinary calls reuse the already loaded native symbols.
  pst_v35_load_cpp()

  pst_v35_cpp_scenario_matrix_coordinates(
    scenario_tree,
    phylo_tree,
    time_vec,
    existing_scenario_ids,
    existing_phylo_ids,
    existing_path_ids,
    path_lookup,
    isTRUE(include_paths),
    as.numeric(time_tolerance)
  )
}

#' Format a transition tree from transition-support maps through C++.
#'
#' A labelled transition simmap tree.
#' Loads the V35 C++ translation unit before crossing the R/C++ boundary.
#' @param transition_support Transition-support list containing edge maps, labels, and tree assembly metadata.
#' @param tip_labels Tip labels aligned to tree tip ids.
#' @param time_tolerance Positive numeric resolution shared with V35 traversal.
#' @return A labelled transition simmap tree.
#' @noRd
pst_v35_cpp_transition_tree_from_support_payload <- function(
    transition_support,
    tip_labels,
    time_tolerance) {
  # Load the C++ symbol immediately before native transition-tree formatting.
  pst_v35_load_cpp()

  # Return the labelled transition simmap from compact transition support.
  pst_v35_cpp_transition_tree_from_support(
    transition_support,
    tip_labels,
    time_tolerance
  )
}

#' Strip V35 size suffixes from a mapped tree through C++.
#'
#' A state-labelled simmap with native mapped-edge exposure.
#' Loads the V35 C++ translation unit before crossing the R/C++ boundary.
#' @param tree Simmap/phylo tree being formatted, relabelled, or canonicalized.
#' @param summary_suffix Optional `S` or `T` suffix used to attach compact SS parts during materialization.
#' @param state_levels Optional biological state domain for attached summary parts.
#' @param include_edge_boundaries Whether topology boundaries count as transitions in attached summary parts.
#' @param source_labels_are_state_only Whether input map labels are already biological states and must not be decoded as scenario state-size labels.
#' @return A state-labelled simmap with native mapped-edge exposure.
#' @noRd
pst_v35_cpp_state_only_tree_payload <- function(
    tree,
    summary_suffix = "",
    state_levels = NULL,
    include_edge_boundaries = FALSE,
    source_labels_are_state_only = FALSE) {
  # Load the C++ symbol immediately before native state-only formatting.
  pst_v35_load_cpp()

  # Return a state-labelled simmap without R-side label stripping or mapped-edge
  # rebuilding. Optional summary parts are attached by C++ while the tree is
  # being materialized, so final SS assembly does not rescan derived PST trees.
  pst_v35_cpp_state_only_tree(
    tree,
    summary_suffix,
    state_levels,
    isTRUE(include_edge_boundaries),
    isTRUE(source_labels_are_state_only)
  )
}

#' Relabel a state tree as root-to-state paths through C++.
#'
#' A path-labelled simmap preserving the input topology and mapped durations.
#' Loads the V35 C++ translation unit before crossing the R/C++ boundary.
#' @param tree Simmap/phylo tree being formatted, relabelled, or canonicalized.
#' @param root_policy Resolved root-anchor policy list.
#' @param sep Separator used when joining state labels into path labels.
#' @return A path-labelled simmap preserving the input topology and mapped durations.
#' @noRd
pst_v35_cpp_pathify_state_tree_payload <- function(tree, root_policy, sep = "|") {
  # Load the C++ symbol immediately before path-labeled tree formatting.
  pst_v35_load_cpp()

  # Return a path-labeled simmap equivalent to the legacy R pathify helper.
  pst_v35_cpp_pathify_state_tree(tree, root_policy, sep)
}

#' Materialize transition size-map families from final cumulative TT rows.
#'
#' @param state_tree Native state-labelled transition tree.
#' @param path_tree Parallel native path-labelled transition tree.
#' @param cumulative_lineage_trans Traversal-owned canonical CLTT path matrix.
#' @param cumulative_scenario_trans Traversal-owned canonical CSTT path matrix.
#' @param path_labels_by_id C++-owned canonical path labels indexed by path id.
#' @param terminal_lineage_counts_by_id Traversal-committed terminal lineage
#'   contributions keyed by canonical path id.
#' @param terminal_scenario_counts_by_id Traversal-committed terminal scenario
#'   contributions keyed by canonical path id.
#' @return Four edge-parallel size-map lists plus a native ownership marker.
#' @noRd
pst_v35_cpp_transition_size_maps_payload <- function(
    state_tree,
    path_tree,
    cumulative_lineage_trans,
    cumulative_scenario_trans,
    path_labels_by_id,
    terminal_lineage_counts_by_id,
    terminal_scenario_counts_by_id,
    terminal_trajectory_group_id_by_tip) {
  # Load the native tree-output symbol before size maps cross the C++ boundary.
  pst_v35_load_cpp()

  # C++ aligns public path labels to final canonical TT columns and assigns
  # every segment; R performs no biological size accounting.
  pst_v35_cpp_transition_size_maps(
    state_tree,
    path_tree,
    cumulative_lineage_trans,
    cumulative_scenario_trans,
    path_labels_by_id,
    terminal_lineage_counts_by_id,
    terminal_scenario_counts_by_id,
    terminal_trajectory_group_id_by_tip
  )
}

#' Assemble SS from C++ tree-materialization parts and online counters.
#'
#' @param tree_parts Named `P`, `S`, and `T` compact metric payloads.
#' @param path_size_maps Traversal-backed transition path support.
#' @param tt Complete or debug-prefix V35 TT object.
#' @return Public V35 SS object.
#' @noRd
pst_v35_cpp_summary_stats_from_parts_payload <- function(
    tree_parts,
    path_size_maps,
    tt = NULL) {
  pst_v35_load_cpp()
  pst_v35_cpp_summary_stats_from_parts(tree_parts, path_size_maps, tt)
}
