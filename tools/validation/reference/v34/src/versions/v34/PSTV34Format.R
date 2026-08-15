# V34 public-output formatting helpers.
#
# C++ traversal data are authoritative. These helpers assign labels, materialize
# public simmap objects, attach traversal-owned size support, and package TT/SS.
# They must never reconstruct biological state from completed matrices.

#' A `phylo` tree with mapped-state metadata removed or simplified.
#' @param phylo_charmap Input `phylo`/`simmap` object whose `maps` and `mapped.edge` fields describe mapped character histories.
#' @return A `phylo` tree with mapped-state metadata removed or simplified.
#' @noRd
pst_v34_fast_char_tree <- function(phylo_charmap) {
  # Preserve simmap char maps directly; V34 fast mode keeps sizes out of labels.
  tree <- list(
    edge = phylo_charmap$edge,
    edge.length = phylo_charmap$edge.length,
    Nnode = phylo_charmap$Nnode,
    tip.label = phylo_charmap$tip.label,
    maps = phylo_charmap$maps,
    mapped.edge = phylo_charmap$mapped.edge
  )
  class(tree) <- class(phylo_charmap)
  tree
}

#' Build the public path lookup from the traversal-owned path registry.
#' @return A path lookup data frame keyed by integer path id.
#' @noRd
pst_v34_path_lookup_from_active <- function(payload, active_summary) {
  labels <- as.character(active_summary$online_size_support$labels_by_id)
  terminal_state_id <- as.integer(active_summary$path_terminal_state_id)
  # C++ owns label generation and terminal-state identity. Both vectors must
  # remain exactly parallel to the canonical path registry.
  if (length(labels) != length(terminal_state_id)) {
    stop("V34 C++ path labels are not aligned to terminal states", call. = FALSE)
  }
  data.frame(
    path_id = seq_along(labels),
    label = labels,
    parent_path_id = as.integer(active_summary$path_parent),
    added_state_id = as.integer(active_summary$path_added_state_id),
    terminal_state_id = terminal_state_id,
    stringsAsFactors = FALSE
  )
}

#' Repackage active C++ traversal output as the canonical R V34 summary shape.
#'
#' A `pst_v34_summary`-compatible list.
#' @param phylo_charmap Input `phylo`/`simmap` object whose `maps` and `mapped.edge` fields describe mapped character histories.
#' @param payload R list payload crossing the V34 R/C++ boundary.
#' @param active_summary C++ active traversal summary to repackage into R-facing structures.
#' @return A `pst_v34_summary`-compatible list.
#' @noRd
pst_v34_summary_from_active <- function(phylo_charmap, payload, active_summary) {
  path_lookup <- pst_v34_path_lookup_from_active(payload, active_summary)
  terminal_records <- active_summary$terminal_records
  # Label committed terminal rows when this trajectory has reached any tips.
  if (nrow(terminal_records)) {
    terminal_records$tip_label <- phylo_charmap$tip.label[terminal_records$tip_id]
    terminal_records <- terminal_records[
      ,
      c("tip_id", "tip_label", "final_phylo_edge_id", "final_scenario_edge_id",
        "stable_scenario_id", "final_path_id", "terminal_time"),
      drop = FALSE
    ]
  }

  # These fields are mandatory traversal products. Refuse to reconstruct them
  # from denser matrices when the C++ ownership contract is broken.
  if (is.null(payload$observed_state_ids) ||
      is.null(active_summary$scenario_edge_records)) {
    stop("V34 traversal summary omitted required owned fields", call. = FALSE)
  }

  out <- list(
    root_policy = payload$root_policy,
    state_lookup = payload$state_lookup,
    observed_state_ids = as.integer(payload$observed_state_ids),
    tip_labels = phylo_charmap$tip.label,
    path_lookup = path_lookup,
    path_terminal_state_id = as.integer(active_summary$path_terminal_state_id),
    final_frontier = active_summary$final_frontier,
    terminal_trajectory_group_id_by_tip =
      active_summary$terminal_trajectory_group_id_by_tip,
    scenario_summary = active_summary$scenario_summary,
    scenario_edge_records = active_summary$scenario_edge_records,
    traversal_size_support = active_summary$traversal_size_support,
    terminal_records = terminal_records,
    # Carry traversal-owned counter rows and event classifications directly to
    # TT serialization. This payload is not rebuilt from scenario matrices.
    online_tt_core = active_summary$online_tt_core,
    # Transition runs are compressed when each batch commits. Tree formatting
    # consumes this payload directly and never scans `state_id_matrix`.
    online_transition_support = active_summary$online_transition_support,
    # Tail, through, terminal, and transition size counters share the same
    # atomic commit boundary as paths and scenarios.
    online_size_support = active_summary$online_size_support,
    # Snapshot mode receives compact C++ event identities aligned one-to-one
    # with committed TT rows. Normal traversal carries `NULL` and allocates no
    # diagnostic event copies.
    debug_step_metadata = active_summary$debug_step_metadata,
    # Preserve C++ validation counts and stage timings as diagnostic evidence;
    # R does not infer or repair these traversal-owned values.
    metadata = active_summary$metadata
  )
  class(out) <- c("pst_v34_active_summary", "pst_v17_summary", class(out))
  out
}

#' Create legacy scenario matrices from an active or canonical V34 summary without recomputing traversal.
#'
#' A list of scenario, edge, id, time, and event matrices.
#' @param summary Canonical V34 traversal summary with lookup tables, scenario matrices, and edge records.
#' @return A list of scenario, edge, id, time, and event matrices.
#' @noRd
pst_v34_scenario_mats_from_summary_fast <- function(summary) {
  state_labels <- character(max(summary$state_lookup$state_id))
  state_labels[summary$state_lookup$state_id] <- as.character(summary$state_lookup$state)
  state_ids <- summary$scenario_summary$state_id_matrix
  scenarios <- matrix(
    state_labels[as.integer(state_ids)],
    nrow = nrow(state_ids),
    ncol = ncol(state_ids),
    dimnames = dimnames(state_ids)
  )
  edge_node_ids <- summary$scenario_summary$edge_node_ids
  edges <- matrix(
    as.character(edge_node_ids),
    nrow = nrow(edge_node_ids),
    ncol = ncol(edge_node_ids),
    dimnames = dimnames(edge_node_ids)
  )

  list(
    scenarios = scenarios,
    edges = edges,
    phylo_edge_ids = summary$scenario_summary$phylo_edge_ids,
    scenario_edge_ids = summary$scenario_summary$scenario_edge_ids,
    time_vec = summary$scenario_summary$time_vec,
    event_vec = summary$scenario_summary$event_vec,
    event_tips = summary$scenario_summary$event_tips
  )
}

#' Capture the traversal-owned scenario ledger before public formatting.
#'
#' @param summary Canonical V34 traversal summary containing scenario records.
#' @return Diagnostic ledger keyed by stable traversal scenario IDs.
#' @noRd
pst_v34_traversal_scenario_ledger <- function(summary) {
  records <- summary$scenario_edge_records
  required_fields <- c(
    "edge_durations",
    "edge_state_ids",
    "edge_path_ids",
    "edge_sizes",
    "edge_parent_scenario_edge_id",
    "edge_type"
  )
  if (!is.list(records) ||
      !all(required_fields %in% names(records))) {
    stop(
      "V34 traversal scenario records cannot form the diagnostic ledger",
      call. = FALSE
    )
  }
  record_count <- length(records$edge_durations)
  parallel_lengths <- vapply(
    records[required_fields],
    length,
    integer(1)
  )
  if (any(parallel_lengths != record_count)) {
    stop(
      "V34 traversal scenario records are not edge-parallel",
      call. = FALSE
    )
  }
  ldif <- summary$online_tt_core$lineage_differentiating_records
  if (!is.data.frame(ldif) ||
      !all(c(
        "parent_scenario_id",
        "child_scenario_id"
      ) %in% names(ldif))) {
    stop(
      "V34 traversal LDIF records lack exact stable scenario identities",
      call. = FALSE
    )
  }

  list(
    scenario_edge_id = seq_len(record_count),
    parent_scenario_edge_id =
      as.integer(records$edge_parent_scenario_edge_id),
    edge_type = as.character(records$edge_type),
    state_ids = records$edge_state_ids,
    path_ids = records$edge_path_ids,
    sizes = records$edge_sizes,
    durations = records$edge_durations,
    ldif = ldif,
    semantic_repairs_applied = 0L
  )
}

#' Resolve one exact traversal stable ID for each formatted public edge.
#'
#' @param provenance Native raw scenario-edge provenance.
#' @param edge_count Number of public scenario edges.
#' @return Integer stable-ID vector; non-bijective cleanup results are `NA`.
#' @noRd
pst_v34_stable_ids_from_provenance <- function(provenance, edge_count) {
  if (is.integer(provenance) && length(provenance) == edge_count) {
    return(provenance)
  }
  if (!is.list(provenance) || length(provenance) != edge_count) {
    stop(
      "V34 formatted scenario provenance is not edge-parallel",
      call. = FALSE
    )
  }
  vapply(
    provenance,
    function(raw_ids) {
      raw_ids <- as.integer(raw_ids)
      if (length(raw_ids) == 1L) raw_ids[[1L]] else NA_integer_
    },
    integer(1)
  )
}

#' Materialize and validate the scenario tree from committed segments.
#'
#' A simmap tree with character labels and attached size-map metadata.
#' @param summary Canonical V34 traversal summary.
#' @param time_tolerance Positive numeric resolution shared with V34 traversal.
#' @return A direct state-and-size-labelled scenario tree.
#' @noRd
pst_v34_build_scenario_tree_fast <- function(summary, time_tolerance) {
  tree <- pst_v34_cpp_scenario_tree_initial_payload(summary)
  pst_v34_cpp_canonicalize_scenario_tree_payload(
    tree,
    summary,
    time_tolerance
  )
}

#' Materialize the transition tree from committed transition runs.
#'
#' A character vector of expanded transition labels.
#' @param phylo_tree Input phylogeny supplying final tip labels.
#' @param summary Canonical V34 traversal summary.
#' @param shrink_stopped_branches Whether terminal STOPPED runs are removed.
#' @return A finalized transition-labelled simmap tree.
#' @noRd
pst_build_transition_tree_v34_fast <- function(
    phylo_tree,
    summary,
    time_tolerance,
    shrink_stopped_branches = TRUE) {
  transition_support <- summary$online_transition_support
  # Every active V34 traversal must return committed transition-run support;
  # silently rebuilding it here would violate traversal ownership.
  if (is.null(transition_support)) {
    stop("V34 traversal did not return online transition support", call. = FALSE)
  }
  trans_tree <- pst_v34_cpp_transition_tree_from_support_payload(
    transition_support,
    phylo_tree$tip.label,
    time_tolerance
  )

  pst_v34_cpp_correct_stopped_branches_payload(
    trans_tree = trans_tree,
    shrink_bl = shrink_stopped_branches
  )
}

#' Attach path-size maps to V34 public trees using the C++ size-support path.
#'
#' A PST tree bundle with size maps attached to scenario and transition trees.
#' @param phylo_tree Input phylogeny used as the physical scaffold for derived trees.
#' @param scenario_tree Optional scenario tree used to align TT or size-support identities to public tree ids.
#' @param trans_tree Transition-labelled simmap tree.
#' @param summary Canonical V34 traversal summary with lookup tables, scenario matrices, and edge records.
#' @return A PST tree bundle with size maps attached to scenario and transition trees.
#' @noRd
pst_v34_attach_size_maps_fast <- function(phylo_tree, scenario_tree, trans_tree, summary) {
  # Compatibility SS tables share the original phylogeny's biological state
  # domain. C++ receives this domain while materializing the derived state-only
  # trees so S/T summary parts are attached before the trees return to R.
  state_levels <- sort(unique(colnames(phylo_tree$mapped.edge)))
  trans_char_tree <- pst_v34_cpp_state_only_tree_payload(
    trans_tree,
    summary_suffix = "T",
    state_levels = state_levels,
    include_edge_boundaries = TRUE
  )
  trans_path_tree <- pst_v34_cpp_pathify_state_tree_payload(
    trans_char_tree,
    root_policy = summary$root_policy
  )
  tt_core <- summary$online_tt_core
  # Internal size maps consume the final row of the traversal-owned matrices
  # later exposed as CLTT/CSTT. Terminal segments consume the terminal lineage
  # and scenario contributions committed with active LTT/STT batches.
  if (is.null(tt_core$cumulative_lin_path) ||
      is.null(tt_core$cumulative_scn_path) ||
      is.null(summary$online_size_support$labels_by_id) ||
      is.null(summary$online_size_support$path_terminal_lineage_counts_by_id) ||
      is.null(summary$online_size_support$path_terminal_scenario_counts_by_id)) {
    stop("V34 traversal did not return cumulative and terminal TT size support", call. = FALSE)
  }
  # Native tree output owns all segment-level support assignment. This call
  # consumes canonical cumulative path matrices and returns edge-parallel public
  # map families; R neither interprets path ancestry nor recalculates sizes.
  transition_size_maps <- pst_v34_cpp_transition_size_maps_payload(
    state_tree = trans_char_tree,
    path_tree = trans_path_tree,
    cumulative_lineage_trans = tt_core$cumulative_lin_path,
    cumulative_scenario_trans = tt_core$cumulative_scn_path,
    path_labels_by_id = summary$online_size_support$labels_by_id,
    terminal_lineage_counts_by_id =
      summary$online_size_support$path_terminal_lineage_counts_by_id,
    terminal_scenario_counts_by_id =
      summary$online_size_support$path_terminal_scenario_counts_by_id,
    terminal_trajectory_group_id_by_tip =
      summary$terminal_trajectory_group_id_by_tip
  )
  trans_tree$size_maps <- list(
    scenario_size_maps = transition_size_maps$scenario_size_maps,
    lineage_size_maps = transition_size_maps$lineage_size_maps
  )
  trans_path_tree$size_maps <- list(
    scenario_size_path_maps = transition_size_maps$scenario_size_path_maps,
    lineage_size_path_maps = transition_size_maps$lineage_size_path_maps
  )
  scenario_char_tree <- pst_v34_cpp_state_only_tree_payload(
    scenario_tree,
    summary_suffix = "S",
    state_levels = state_levels,
    include_edge_boundaries = TRUE,
    source_labels_are_state_only = TRUE
  )
  scenario_char_tree$size_maps <- NULL
  list(
    phylo = phylo_tree,
    scenario = scenario_tree,
    trans = trans_tree,
    path_trees = list(
      trans = trans_path_tree
    ),
    path_dictionary = NULL,
    char_only_trees = list(
      scenario = scenario_char_tree,
      trans = trans_char_tree
    ),
    transition_size_map_construction = transition_size_maps$construction_phase,
    root_policy = summary$root_policy
  )
}

#' Assemble the C++-first V34 public output from active traversal summary and requested options.
#'
#' The public V34 PST result list with trees, scenario matrices, optional TT, optional SS, and metadata.
#' This is a construction boundary: it preserves public V34 object shapes while allowing internal traversal ids to remain compact integers.
#' @param phylo_charmap Input `phylo`/`simmap` object whose `maps` and `mapped.edge` fields describe mapped character histories.
#' @param payload R list payload crossing the V34 R/C++ boundary.
#' @param active_summary C++ active traversal summary to repackage into R-facing structures.
#' @param tt_mats Logical flag requesting through-time matrix construction.
#' @param ss Logical flag requesting empirical summary-statistics construction.
#' @return The public V34 PST result list with trees, scenario matrices, optional TT, optional SS, and metadata.
#' @noRd
pst_v34_build_cpp_traversal_output <- function(
    phylo_charmap,
    payload,
    active_summary,
    tt_mats = TRUE,
    ss = TRUE) {
  output_started <- proc.time()[["elapsed"]]
  # Assemble public surfaces from the dedicated C++ traversal summary. The C++
  # side owns integer PST/TT/SS support; this R boundary only converts integer
  # ids to state/path names and packages simmap-compatible objects.
  phylo_tree <- pst_v34_fast_char_tree(phylo_charmap)
  summary <- pst_v34_summary_from_active(phylo_charmap, payload, active_summary)
  state_lookup <- summary$state_lookup
  path_lookup <- summary$path_lookup
  traversal_scenario_ledger <-
    pst_v34_traversal_scenario_ledger(summary)

  scenario_state_size_tree <- pst_v34_build_scenario_tree_fast(
    summary,
    time_tolerance = payload$options$time_tolerance
  )
  summary <- pst_v34_cpp_canonicalize_summary_scenario_ids_payload(
    summary,
    scenario_state_size_tree
  )
  trans_state_size_tree <- pst_build_transition_tree_v34_fast(
    phylo_tree = phylo_tree,
    summary = summary,
    time_tolerance = payload$options$time_tolerance
  )
  size_map_started <- proc.time()[["elapsed"]]
  size_attached <- pst_v34_attach_size_maps_fast(
    phylo_tree = phylo_tree,
    scenario_tree = scenario_state_size_tree,
    trans_tree = trans_state_size_tree,
    summary = summary
  )
  size_map_seconds <- proc.time()[["elapsed"]] - size_map_started
  scenario_state_size_tree <- size_attached$scenario
  trans_state_size_tree <- size_attached$trans
  scenario_tree <- size_attached$char_only_trees$scenario
  scenario_tree$size_maps <- scenario_state_size_tree$size_maps
  scenario_tree$terminal_events <- scenario_state_size_tree$terminal_events
  scenario_edge_provenance <- scenario_state_size_tree$raw_scenario_edge_ids
  if (is.null(scenario_edge_provenance)) {
    scenario_edge_provenance <- attr(
      scenario_state_size_tree,
      "pst_v34_raw_scenario_edge_ids",
      exact = TRUE
    )
  }
  attr(scenario_tree, "pst_v34_raw_scenario_edge_ids") <-
    scenario_edge_provenance
  attr(scenario_tree, "pst_v34_stable_scenario_ids") <-
    pst_v34_stable_ids_from_provenance(
      scenario_edge_provenance,
      nrow(scenario_tree$edge)
    )
  scenario_tree$raw_scenario_edge_ids <- NULL
  trans_tree <- size_attached$char_only_trees$trans
  trans_tree$size_maps <- trans_state_size_tree$size_maps
  trans_tree$path_size_maps <- size_attached$path_trees$trans$size_maps
  # Native materialization retains raw scenario-edge provenance in the
  # V34-named diagnostic attribute. Do not also expose the edge-parallel list
  # as a new public transition-tree field.
  trans_tree$raw_scenario_edge_ids <- NULL
  scenario_mats <- if (
      length(attr(scenario_state_size_tree, "pst_v34_terminal_trajectory_merged_labels")) > 0L
  ) {
    # Terminal-trajectory merges change public scenario ids, so regenerate only
    # the label/id formatting from the already committed summary matrices.
    pst_v34_scenario_mats_from_summary_fast(summary)
  } else {
    # Without terminal merging, use the C++ matrices directly; absence is an
    # ownership failure rather than a reason to reconstruct traversal state.
    if (is.null(active_summary$scenario_mats)) {
      stop("V34 traversal omitted owned scenario matrices", call. = FALSE)
    }
    active_summary$scenario_mats
  }

  # Capture exact P/S/T compatibility metrics while all three C++-owned PST
  # projections are finalized. S/T parts are attached by native derived-tree
  # materialization; only P is captured here because the original phylogeny is
  # caller input rather than a V34-derived PST projection.
  # Legacy SS transition matrices use the sorted state domain represented by
  # the phylogeny maps. Synthetic root labels that never occupy a phylogeny
  # segment must not create zero-only compatibility rows or columns.
  state_levels <- sort(unique(colnames(phylo_tree$mapped.edge)))
  # The original phylogeny is already state-only, but passing it through the
  # native state formatter lets P compatibility metrics be captured inside the
  # same C++ materialization boundary as S and T. No completed public tree is
  # handed to a separate post-traversal summary scanner.
  phylo_tree <- pst_v34_cpp_state_only_tree_payload(
    phylo_tree,
    summary_suffix = "P",
    state_levels = state_levels,
    include_edge_boundaries = FALSE,
    source_labels_are_state_only = TRUE
  )
  phylo_summary_parts <- attr(phylo_tree, "pst_v34_tree_summary_parts")
  scenario_summary_parts <- attr(scenario_tree, "pst_v34_tree_summary_parts")
  trans_summary_parts <- attr(trans_tree, "pst_v34_tree_summary_parts")
  # Missing derived-tree metric parts mean the native materialization boundary
  # failed to attach its dependent SS payload; do not repair by rescanning here.
  if (is.null(phylo_summary_parts) ||
      is.null(scenario_summary_parts) ||
      is.null(trans_summary_parts)) {
    stop("V34 derived tree materialization omitted summary parts", call. = FALSE)
  }
  tree_summary_parts <- list(
    P = phylo_summary_parts,
    S = scenario_summary_parts,
    T = trans_summary_parts
  )
  # SS receives only compact dependency metrics captured as each native tree
  # was formatted. Later SS assembly may name tables and calculate ratios, but
  # it cannot inspect P/S/T maps or scenario matrices.
  summary$online_ss_support <- list(
    tree_parts = tree_summary_parts,
    construction_phase = "native_dependency_formatting"
  )
  pst_materialized <- proc.time()[["elapsed"]]

  out <- list(
    phylo = phylo_tree,
    scenario = scenario_tree,
    trans = trans_tree,
    scenario_mats = scenario_mats,
    tt = NULL,
    ss = NULL,
    root_policy = size_attached$root_policy
  )

  attr(out, "pst_v34_cpp_traversal_summary") <- summary
  attr(out, "pst_v34_traversal_scenario_ledger") <-
    traversal_scenario_ledger
  attr(out, "pst_v34_state_lookup") <- state_lookup
  attr(out, "pst_v34_path_lookup") <- path_lookup
  attr(out, "pst_v34_trans_path_size_maps") <- out$trans$path_size_maps
  attr(out, "pst_v34_tree_summary_parts") <- tree_summary_parts
  attr(out, "pst_v34_transition_size_map_construction") <-
    size_attached$transition_size_map_construction
  tt_seconds <- 0
  ss_seconds <- 0
  # Materialize TT directly from traversal-owned online rows when requested.
  if (isTRUE(tt_mats)) {
    tt_started <- proc.time()[["elapsed"]]
    out$tt <- pst_v34_build_tt_canonical(
      summary,
      time_tolerance = payload$options$time_tolerance
    )
    tt_seconds <- proc.time()[["elapsed"]] - tt_started
  }
  # Build optional summary statistics only when requested by public API flags.
  if (isTRUE(ss)) {
    ss_started <- proc.time()[["elapsed"]]
    tt_in_ss_seconds <- 0
    # SS requires TT counters; materialize them here only when the caller did
    # not request the public TT surface separately.
    if (is.null(out$tt)) {
      tt_started <- proc.time()[["elapsed"]]
      out$tt <- pst_v34_build_tt_canonical(
        summary,
        time_tolerance = payload$options$time_tolerance
      )
      tt_in_ss_seconds <- proc.time()[["elapsed"]] - tt_started
      tt_seconds <- tt_seconds + tt_in_ss_seconds
    }
    out$ss <- pst_v34_build_summary_stats_canonical(
      tt = out$tt,
      tree_parts = tree_summary_parts,
      path_size_maps = out$trans$path_size_maps
    )
    ss_seconds <- proc.time()[["elapsed"]] - ss_started - tt_in_ss_seconds
  }

  output_finished <- proc.time()[["elapsed"]]
  attr(out, "pst_v34_materialization_timings") <- list(
    tree_output = max(
      pst_materialized - output_started - size_map_seconds,
      0
    ),
    size_map = max(size_map_seconds, 0),
    tt_serialization = tt_seconds,
    ss = max(ss_seconds, 0),
    pst_materialization = pst_materialized - output_started,
    tt_materialization = tt_seconds,
    ss_materialization = max(ss_seconds, 0),
    r_formatting_total = output_finished - output_started
  )

  attr(out, "pst_v34_cpp_traversal_contract") <- list(
    source = "cpp_active_traversal_summary",
    maps = "char",
    size_maps = "numeric_sizes",
    state_labels = "state_lookup",
    path_labels = "path_lookup",
    tt_source = "online_traversal_counters",
    batch_semantics = "atomic_same_time"
  )
  out
}
