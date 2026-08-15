#' Build the V34 public PST trajectory object from a simmap tree.
#'
#' A list of V34 public PST trees plus optional TT and summary-statistics payloads.
#' This is a construction boundary: it preserves public V34 object shapes while allowing internal traversal ids to remain compact integers.
#' @param phylo_charmap Input `phylo`/`simmap` object whose `maps` and `mapped.edge` fields describe mapped character histories.
#' @param root_state Optional explicit root-anchor state; when `NULL`, V34 infers or synthesizes one from root-descending edges.
#' @param tt_mats Logical flag requesting through-time matrix construction.
#' @param ss Logical flag requesting empirical summary-statistics construction.
#' @param include_tt Compatibility flag overriding `tt_mats` when supplied.
#' @param include_ss Compatibility flag overriding `ss` when supplied.
#' @param include_legacy_summary Logical flag exposing the non-authoritative
#'   legacy `ss$summary` compatibility payload. It is omitted by default.
#' @param include_scenario_mats Logical flag controlling whether the completed
#'   scenario matrices are retained in the returned trajectory. The default
#'   `TRUE` includes finalized scenario-edge ids and both exact step-id matrices.
#' @param include_path_maps Logical flag controlling whether edge-parallel path
#'   maps and dense scenario/phylogeny path matrices are retained. The default
#'   `TRUE` makes every matrix cell's resolved path directly inspectable.
#' @param fast_labelled Compatibility flag accepted by V34; the C++-first path has no separate fast-labelled record path.
#' @param debug Logical flag that attaches low-level C++ round-trip diagnostics to the result.
#' @param time_tolerance Positive numeric resolution used consistently for
#'   event-time batches, scenario-tree cleanup, and matrix/tree comparisons.
#' @param ... Reserved compatibility arguments. V34 reads `debug_snapshots`,
#'   `debug_snapshot_scope`, `debug_snapshot_steps`, and
#'   `debug_snapshot_times` to materialize selected batch trajectory prefixes.
#' @return A list of V34 public PST trees plus optional TT and summary-statistics payloads.
#' @export
get_pst_trees_v34 <- function(
    phylo_charmap,
    root_state = NULL,
    tt_mats = TRUE,
    ss = TRUE,
    include_tt = NULL,
    include_ss = NULL,
    include_legacy_summary = FALSE,
    include_scenario_mats = TRUE,
    include_path_maps = TRUE,
    fast_labelled = FALSE,
    debug = FALSE,
    time_tolerance = 1e-6,
    ...) {
  total_started <- proc.time()[["elapsed"]]
  if (!is.numeric(time_tolerance) || length(time_tolerance) != 1L ||
      !is.finite(time_tolerance) || time_tolerance <= 0) {
    stop("time_tolerance must be one finite, positive numeric value", call. = FALSE)
  }
  time_tolerance <- as.numeric(time_tolerance)
  debug_options <- list(...)
  # V34 previously accepted a second `tolerance` argument for scenario-tree
  # cleanup. Silently accepting that name through `...` would make callers
  # believe they can still configure two different resolutions, so reject it
  # explicitly and direct every time-sensitive operation through the one
  # `time_tolerance` value above.
  if ("tolerance" %in% names(debug_options)) {
    stop(
      "V34 uses one time_tolerance argument; replace tolerance with time_tolerance",
      call. = FALSE
    )
  }
  # Read V34-only debug controls from `...` without changing the established
  # public wrapper signature that existing callers use.
  # Use the supplied flag only when callers named it; otherwise disable snapshots.
  debug_snapshots <- if ("debug_snapshots" %in% names(debug_options)) {
    isTRUE(debug_options$debug_snapshots)
  } else {
    FALSE
  }
  # Use the supplied scope only when callers named it; otherwise retain the batch default.
  debug_snapshot_scope <- if ("debug_snapshot_scope" %in% names(debug_options)) {
    debug_options$debug_snapshot_scope
  } else {
    "batch"
  }
  # Preserve named step selectors when present; otherwise leave future selection unset.
  debug_snapshot_steps <- if ("debug_snapshot_steps" %in% names(debug_options)) {
    debug_options$debug_snapshot_steps
  } else {
    NULL
  }
  # Preserve named time selectors when present; otherwise leave future selection unset.
  debug_snapshot_times <- if ("debug_snapshot_times" %in% names(debug_options)) {
    debug_options$debug_snapshot_times
  } else {
    NULL
  }
  # This private option exists only to exercise rollback deterministically in
  # the debug contract test; production callers should leave it unset.
  debug_fail_batch_id <- if (".debug_fail_batch_id" %in% names(debug_options)) {
    as.integer(debug_options$.debug_fail_batch_id)
  } else {
    NA_integer_
  }
  # Reject snapshot capture when normal debug output is disabled because the
  # snapshot implementation depends on the debug construction path.
  if (isTRUE(debug_snapshots) && !isTRUE(debug)) {
    stop("debug_snapshots=TRUE requires debug=TRUE", call. = FALSE)
  }
  # Keep the first V34 scaffold limited to the batch scope while preserving
  # the explicit API contract for future per-step and per-time snapshots.
  if (!identical(debug_snapshot_scope, "batch")) {
    stop("debug_snapshot_scope must be \"batch\"", call. = FALSE)
  }
  # Failure recovery is meaningful only when full debug snapshots were
  # requested; normal calls must continue to throw traversal errors.
  if (!is.na(debug_fail_batch_id) &&
      (!isTRUE(debug) || !isTRUE(debug_snapshots))) {
    stop(
      ".debug_fail_batch_id requires debug=TRUE and debug_snapshots=TRUE",
      call. = FALSE
    )
  }
  # Honor the legacy include flag when callers use the older TT option name.
  if (!is.null(include_tt)) {
    tt_mats <- isTRUE(include_tt)
  }
  # Honor the legacy include flag when callers use the older summary-statistics option name.
  if (!is.null(include_ss)) {
    ss <- isTRUE(include_ss)
  }
  # Warn on the accepted compatibility flag because V34 always uses the C++-first labelled path.
  if (isTRUE(fast_labelled)) {
    warning(
      "V34 is C++-first by default; fast_labelled is accepted for compatibility and has no separate record path.",
      call. = FALSE
    )
  }
  input_started <- proc.time()[["elapsed"]]
  payload <- pst_v34_encode_input_metadata(
    phylo_charmap,
    root_state = root_state,
    tt_mats = tt_mats,
    ss = ss,
    debug = debug,
    debug_snapshots = debug_snapshots,
    debug_snapshot_scope = debug_snapshot_scope,
    debug_snapshot_steps = debug_snapshot_steps,
    debug_snapshot_times = debug_snapshot_times,
    time_tolerance = time_tolerance
  )
  input_seconds <- proc.time()[["elapsed"]] - input_started
  traversal_started <- proc.time()[["elapsed"]]
  active_summary <- pst_v34_cpp_active_summary_tree_payload(
    phylo_charmap,
    payload$state_lookup,
    payload$root_anchor_state_ids,
    payload$root_policy$synthetic,
    recover_debug_failure = debug_snapshots,
    debug_fail_batch_id = debug_fail_batch_id,
    time_tolerance = time_tolerance
  )
  traversal_seconds <- proc.time()[["elapsed"]] - traversal_started
  requested_tt <- isTRUE(tt_mats)
  requested_ss <- isTRUE(ss)
  result <- pst_v34_build_cpp_traversal_output(
    phylo_charmap,
    payload,
    active_summary,
    tt_mats = requested_tt || debug_snapshots,
    ss = requested_ss || debug_snapshots
  )
  # Preserve the validated request and materialize selected batch snapshots only
  # on the explicit debug path. Normal traversal allocates no snapshot copies.
  attr(result, "pst_v34_debug_snapshot_request") <- payload$options$debug_snapshots
  debug_failure <- isTRUE(active_summary$metadata$debug_failed)
  if (isTRUE(debug_snapshots)) {
    result <- pst_v34_attach_debug_snapshots(
      result,
      steps = if (debug_failure) NULL else debug_snapshot_steps,
      times = if (debug_failure) NULL else debug_snapshot_times,
      time_tolerance = time_tolerance
    )
  }

  # A recovered traversal returns the complete trajectory at the last committed
  # checkpoint, never the partially mutated attempted batch.
  if (debug_failure) {
    snapshots <- attr(result, "pst_v34_debug_trajectory_steps")
    step_index <- attr(result, "pst_v34_debug_step_index")
    if (!length(snapshots)) {
      stop("V34 debug recovery produced no committed trajectory", call. = FALSE)
    }
    recovered <- snapshots[[length(snapshots)]]
    attr(recovered, "pst_v34_debug_trajectory_steps") <- snapshots
    attr(recovered, "pst_v34_debug_step_index") <- step_index
    attr(recovered, "pst_v34_debug_status") <- list(
      status = "failed",
      scope = "batch",
      saved_steps = step_index$step_id,
      attempted_batch_id = active_summary$metadata$debug_attempted_batch_id,
      attempted_time = active_summary$metadata$debug_attempted_time,
      message = active_summary$metadata$debug_failure_message
    )
    class(recovered) <- unique(c("pst_v34_debug_failure", class(recovered)))
    result <- recovered
  }

  # The final public object obeys caller inclusion flags even though snapshot
  # mode internally required both surfaces for its complete trajectory objects.
  if (!requested_tt) {
    result$tt <- NULL
  }
  if (!requested_ss) {
    result$ss <- NULL
  }
  if (!isTRUE(include_legacy_summary) && !is.null(result$ss)) {
    result$ss$summary <- NULL
  }

  # Scenario matrices are an optional public inspection surface. Their exact
  # edge/step coordinates are built only when the caller retains that surface;
  # TT and SS may still use traversal-owned matrices internally above.
  if (isTRUE(include_scenario_mats)) {
    result <- pst_v34_complete_scenario_mats(
      result,
      include_paths = isTRUE(include_path_maps),
      time_tolerance = time_tolerance
    )
  } else {
    result$scenario_mats <- NULL
    result$scenario_mats_active <- NULL
    # Phylogeny path maps remain useful without dense matrices, so attach them
    # independently when the caller requested path-map output.
    if (isTRUE(include_path_maps)) {
      result <- pst_v34_attach_phylo_path_maps(result)
    }
  }

  # Path maps are independently optional. Internal construction above has
  # already consumed scenario paths where required for TT/SS, so removing these
  # returned fields cannot change any computed tree or statistic.
  if (!isTRUE(include_path_maps)) {
    result$phylo$path_maps <- NULL
    result$scenario$path_maps <- NULL
    if (!is.null(result$scenario_mats)) {
      result$scenario_mats$scenario_paths <- NULL
      result$scenario_mats$phylo_paths <- NULL
    }
    if (!is.null(result$scenario_mats_active)) {
      result$scenario_mats_active$scenario_paths <- NULL
      result$scenario_mats_active$phylo_paths <- NULL
    }
  }

  # Attach the compact C++ input round-trip only when low-level debug is asked
  # for.
  if (isTRUE(debug)) {
    debug_payload <- pst_v34_encode_input(
      phylo_charmap,
      root_state = root_state,
      tt_mats = tt_mats,
      ss = ss,
      debug = debug,
      debug_snapshots = debug_snapshots,
      debug_snapshot_scope = debug_snapshot_scope,
      debug_snapshot_steps = debug_snapshot_steps,
      debug_snapshot_times = debug_snapshot_times,
      time_tolerance = time_tolerance
    )
    attr(result, "pst_v34_cpp_roundtrip") <- pst_v34_cpp_roundtrip_payload(debug_payload)
  }

  materialization_timings <- attr(result, "pst_v34_materialization_timings")
  cpp_timings <- active_summary$metadata$stage_timings
  attr(result, "pst_v34_stage_timings") <- c(
    list(
      input_validation = input_seconds,
      cpp_bridge_total = traversal_seconds
    ),
    cpp_timings,
    materialization_timings,
    list(total = proc.time()[["elapsed"]] - total_started)
  )
  # Preserve the one construction tolerance on the completed trajectory so
  # reference tests and later matrix/tree comparisons use the same resolution.
  attr(result, "pst_v34_time_tolerance") <- time_tolerance

  result
}

#' Compatibility alias for callers that request the V34 trajectory object by name.
#'
#' The same object returned by `get_pst_trees_v34()`.
#' @param phylo_charmap Input `phylo`/`simmap` object whose `maps` and `mapped.edge` fields describe mapped character histories.
#' @param root_state Optional explicit root-anchor state; when `NULL`, V34 infers or synthesizes one from root-descending edges.
#' @param include_tt Logical flag requesting through-time matrix construction.
#' @param include_ss Logical flag requesting empirical summary-statistics construction.
#' @param include_legacy_summary Logical flag exposing the non-authoritative
#'   legacy `ss$summary` compatibility payload. It is omitted by default.
#' @param include_scenario_mats Logical flag controlling whether completed
#'   scenario matrices are returned. Defaults to `TRUE`.
#' @param include_path_maps Logical flag controlling whether phylogeny/scenario
#'   path maps and dense path matrices are returned. Defaults to `TRUE`.
#' @param fast_labelled Compatibility flag accepted by V34; the C++-first path has no separate fast-labelled record path.
#' @param debug Logical flag that attaches low-level C++ round-trip diagnostics to the result.
#' @param time_tolerance Positive numeric resolution used consistently for
#'   event-time batches, scenario-tree cleanup, and matrix/tree comparisons.
#' @param run_tests Logical flag that runs the scenario-matrix, transition-size,
#'   scenario-edge-type, and SS `lindiff` reference tests on the completed
#'   trajectory. The resulting report is attached as `reference_tests`.
#' @param ... Reserved compatibility arguments, including V34 debug-snapshot
#'   request controls handled by `get_pst_trees_v34()`.
#' @return The same object returned by `get_pst_trees_v34()`.
#' @export
get_trajectory_obj_v34 <- function(
    phylo_charmap,
    root_state = NULL,
    include_tt = TRUE,
    include_ss = TRUE,
    include_legacy_summary = FALSE,
    include_scenario_mats = TRUE,
    include_path_maps = TRUE,
    fast_labelled = FALSE,
    debug = FALSE,
    time_tolerance = 1e-6,
    run_tests = FALSE,
    ...) {
  # Reference tests require the complete coordinate and path surface even when
  # the caller does not want to retain it. Build that temporary surface for the
  # tests, then remove only the fields excluded by the public request below.
  build_test_surface <- isTRUE(run_tests)
  trajectory_obj <- get_pst_trees_v34(
    phylo_charmap = phylo_charmap,
    root_state = root_state,
    tt_mats = isTRUE(include_tt),
    ss = isTRUE(include_ss),
    include_legacy_summary = isTRUE(include_legacy_summary),
    include_scenario_mats = isTRUE(include_scenario_mats) || build_test_surface,
    include_path_maps = isTRUE(include_path_maps) || build_test_surface,
    fast_labelled = fast_labelled,
    debug = debug,
    time_tolerance = time_tolerance,
    ...
  )

  # Run the global-label independent field-to-field reference checks only when
  # the caller explicitly requests them. The default construction path
  # therefore retains its existing output and avoids the additional
  # matrix-processing cost.
  if (isTRUE(run_tests)) {
    reports <- run_v34_reference_suites(
      phylo_charmap = phylo_charmap,
      trajectory_obj = trajectory_obj,
      root_state = root_state,
      time_tolerance = time_tolerance
    )
    trajectory_obj$reference_tests <- reports$new
    trajectory_obj$reference_tests_old <- reports$old
    trajectory_obj$reference_tests_new <- reports$new
    trajectory_obj$reference_test_agreement <- reports$agreement
    trajectory_obj$reference_tests_all_passed <- reports$all_passed
    trajectory_obj$reference_test_terminal_extensions <-
      reports$terminal_extensions
    trajectory_obj$reference_test_input_classification <-
      reports$input_classification
    trajectory_obj$reference_test_historical_input_modified <-
      reports$historical_input_modified
    trajectory_obj$reference_test_time_tolerance <-
      reports$time_tolerance
    trajectory_obj$reference_test_tip_time_keys <-
      reports$tip_time_keys
  }

  # Honor omission flags after optional validation. This changes only the
  # returned inspection fields; the already-computed PST trees, TT, SS, and
  # reference-test report remain untouched.
  if (!isTRUE(include_path_maps)) {
    trajectory_obj$phylo$path_maps <- NULL
    trajectory_obj$scenario$path_maps <- NULL
    if (!is.null(trajectory_obj$scenario_mats)) {
      trajectory_obj$scenario_mats$scenario_paths <- NULL
      trajectory_obj$scenario_mats$phylo_paths <- NULL
    }
    if (!is.null(trajectory_obj$scenario_mats_active)) {
      trajectory_obj$scenario_mats_active$scenario_paths <- NULL
      trajectory_obj$scenario_mats_active$phylo_paths <- NULL
    }
  }
  if (!isTRUE(include_scenario_mats)) {
    trajectory_obj$scenario_mats <- NULL
    trajectory_obj$scenario_mats_active <- NULL
  }

  trajectory_obj
}

#' Replace only `scenario_mats$scenario_edge_ids` with finalized edge row ids.
#'
#' This compatibility helper remains available for older saved V34 objects that
#' contain traversal ids but predate the complete coordinate matrices. New
#' construction uses `pst_v34_complete_scenario_mats()` below.
#' @param trajectory_obj Complete V34 trajectory object.
#' @param time_tolerance Positive V34 time-key resolution. By default this is
#'   read from the trajectory metadata and falls back to `1e-6` for older objects.
#' @return The trajectory object with only `scenario_mats$scenario_edge_ids` replaced.
#' @noRd
pst_v34_finalize_scenario_edge_ids <- function(
    trajectory_obj,
    time_tolerance = attr(
      trajectory_obj,
      "pst_v34_time_tolerance",
      exact = TRUE
    )) {
  # Older saved objects may predate the tolerance attribute. Use the current
  # V34 default rather than introducing a second reconstruction resolution.
  if (is.null(time_tolerance)) {
    time_tolerance <- 1e-6
  }
  if (!is.numeric(time_tolerance) || length(time_tolerance) != 1L ||
      !is.finite(time_tolerance) || time_tolerance <= 0) {
    stop("time_tolerance must be one finite, positive numeric value", call. = FALSE)
  }

  # Validate only the fields consumed by the pure C++ matrix builder. Missing
  # inputs stop before R assigns anything into the trajectory copy.
  if (is.null(trajectory_obj$scenario) ||
      is.null(trajectory_obj$phylo$tip.label) ||
      is.null(trajectory_obj$scenario_mats$time_vec) ||
      is.null(trajectory_obj$scenario_mats$scenario_edge_ids)) {
    stop("Trajectory lacks fields required to finalize scenario-edge ids", call. = FALSE)
  }

  # C++ returns one new matrix. This is the only assignment performed by the
  # helper, which makes the additive field-level footprint explicit.
  trajectory_obj$scenario_mats$scenario_edge_ids <-
    pst_v34_cpp_final_scenario_edge_ids_payload(
      scenario_tree = trajectory_obj$scenario,
      phylo_tip_labels = trajectory_obj$phylo$tip.label,
      time_vec = trajectory_obj$scenario_mats$time_vec,
      existing_ids = trajectory_obj$scenario_mats$scenario_edge_ids,
      time_tolerance = as.numeric(time_tolerance)
    )

  trajectory_obj
}

#' Attach root-to-state path maps to the returned phylogeny.
#'
#' V34 already uses the native pathifier for derived transition-tree work. This
#' helper applies the same root policy and escaping rules to the input phylogeny,
#' then copies only its edge-parallel path-labelled maps onto the public state
#' tree. Topology, state maps, lengths, labels, and every non-phylogeny field are
#' left unchanged.
#' @param trajectory_obj Complete V34 trajectory object.
#' @return Trajectory with `phylo$path_maps` attached.
#' @noRd
pst_v34_attach_phylo_path_maps <- function(trajectory_obj) {
  if (is.null(trajectory_obj$phylo) || is.null(trajectory_obj$root_policy)) {
    stop("Trajectory lacks phylogeny or root policy for path maps", call. = FALSE)
  }

  phylo_path_tree <- pst_v34_cpp_pathify_state_tree_payload(
    trajectory_obj$phylo,
    root_policy = trajectory_obj$root_policy
  )
  trajectory_obj$phylo$path_maps <- phylo_path_tree$maps
  trajectory_obj
}

#' Complete the saved scenario-matrix addressing surface.
#'
#' One native pass rebuilds public scenario-edge ids and adds exact step ids for
#' both scenario and phylogeny edges. When paths are requested, the same native
#' pass resolves the corresponding path labels from edge-parallel path maps.
#' Assignments are confined to the three coordinate matrices, two optional path
#' matrices, and the optional phylogeny path-map list.
#' @param trajectory_obj Complete V34 trajectory with traversal scenario mats.
#' @param include_paths Whether to attach/return path maps and dense path values.
#' @param time_tolerance Positive shared V34 integer-time-key resolution.
#' @return Trajectory with a finalized, step-addressable `scenario_mats` list.
#' @noRd
pst_v34_complete_scenario_mats <- function(
    trajectory_obj,
    include_paths = TRUE,
    time_tolerance = attr(
      trajectory_obj,
      "pst_v34_time_tolerance",
      exact = TRUE
    )) {
  completion_started <- proc.time()[["elapsed"]]
  if (is.null(time_tolerance)) {
    time_tolerance <- 1e-6
  }
  if (!is.numeric(time_tolerance) || length(time_tolerance) != 1L ||
      !is.finite(time_tolerance) || time_tolerance <= 0) {
    stop("time_tolerance must be one finite, positive numeric value", call. = FALSE)
  }
  if (is.null(trajectory_obj$scenario) ||
      is.null(trajectory_obj$phylo) ||
      is.null(trajectory_obj$scenario_mats$time_vec) ||
      is.null(trajectory_obj$scenario_mats$scenario_edge_ids) ||
      is.null(trajectory_obj$scenario_mats$phylo_edge_ids)) {
    stop("Trajectory lacks fields required for scenario-matrix coordinates", call. = FALSE)
  }

  # Paths are constructed before crossing the coordinate boundary because the
  # native step resolver reads edge-parallel path labels only when requested.
  phylo_path_started <- proc.time()[["elapsed"]]
  if (isTRUE(include_paths)) {
    trajectory_obj <- pst_v34_attach_phylo_path_maps(trajectory_obj)
  }
  phylo_path_seconds <- proc.time()[["elapsed"]] - phylo_path_started

  coordinate_started <- proc.time()[["elapsed"]]
  coordinates <- pst_v34_cpp_scenario_matrix_coordinates_payload(
    scenario_tree = trajectory_obj$scenario,
    phylo_tree = trajectory_obj$phylo,
    time_vec = trajectory_obj$scenario_mats$time_vec,
    existing_scenario_ids = trajectory_obj$scenario_mats$scenario_edge_ids,
    existing_phylo_ids = trajectory_obj$scenario_mats$phylo_edge_ids,
    include_paths = isTRUE(include_paths),
    time_tolerance = as.numeric(time_tolerance)
  )
  coordinate_seconds <- proc.time()[["elapsed"]] - coordinate_started

  # These are the only scenario-matrix assignments made by completion. Existing
  # states, times, events, phylogeny-edge ids, and their dimnames are preserved.
  assignment_started <- proc.time()[["elapsed"]]
  trajectory_obj$scenario_mats$scenario_edge_ids <-
    coordinates$scenario_edge_ids
  trajectory_obj$scenario_mats$scenario_edge_step_ids <-
    coordinates$scenario_edge_step_ids
  trajectory_obj$scenario_mats$phylo_edge_step_ids <-
    coordinates$phylo_edge_step_ids
  if (isTRUE(include_paths)) {
    trajectory_obj$scenario_mats$scenario_paths <- coordinates$scenario_paths
    trajectory_obj$scenario_mats$phylo_paths <- coordinates$phylo_paths
  }
  assignment_seconds <- proc.time()[["elapsed"]] - assignment_started

  active_mask_started <- proc.time()[["elapsed"]]
  trajectory_obj$scenario_mats_active <- pst_v34_mask_inactive_scenario_mats(
    scenario_mats = trajectory_obj$scenario_mats,
    terminal_events = trajectory_obj$scenario$terminal_events,
    time_tolerance = as.numeric(time_tolerance)
  )
  active_mask_seconds <- proc.time()[["elapsed"]] - active_mask_started

  # Add the previously unmeasured completion work to the same timing list used
  # by `get_pst_trees_v34()`. These labels make the full constructor total
  # reconcilable without changing any returned matrix, path, or tree value.
  materialization_timings <- attr(
    trajectory_obj,
    "pst_v34_materialization_timings",
    exact = TRUE
  )
  attr(trajectory_obj, "pst_v34_materialization_timings") <- c(
    materialization_timings,
    list(
      phylo_path_map_completion = phylo_path_seconds,
      scenario_matrix_coordinate_completion = coordinate_seconds,
      scenario_matrix_coordinate_assignment = assignment_seconds,
      scenario_matrix_active_mask = active_mask_seconds,
      scenario_matrix_completion =
        proc.time()[["elapsed"]] - completion_started
    )
  )

  trajectory_obj
}

#' Clone scenario matrices and mask cells strictly after each tip endpoint.
#'
#' @param scenario_mats Final endpoint-inclusive V34 scenario-matrix list.
#' @param terminal_events One canonical terminal-event row per phylogenetic tip.
#' @param time_tolerance Positive V34 integer-time-key resolution.
#' @return A same-shaped scenario-matrix list with post-terminal matrix cells `NA`.
#' @noRd
pst_v34_mask_inactive_scenario_mats <- function(
    scenario_mats,
    terminal_events,
    time_tolerance) {
  if (!is.list(scenario_mats) || !is.numeric(scenario_mats$time_vec) ||
      !is.data.frame(terminal_events)) {
    stop("V34 active scenario matrices require final matrices and terminal events", call. = FALSE)
  }
  matrix_names <- names(scenario_mats)[vapply(
    scenario_mats,
    is.matrix,
    logical(1)
  )]
  time_keys <- pst_v34_cpp_batch_time_keys_payload(
    scenario_mats$time_vec,
    time_tolerance
  )
  active <- scenario_mats
  for (matrix_name in matrix_names) {
    matrix_value <- scenario_mats[[matrix_name]]
    if (nrow(matrix_value) != nrow(terminal_events) ||
        ncol(matrix_value) != length(time_keys)) {
      stop("V34 scenario matrix is not aligned to terminal events and times", call. = FALSE)
    }
    for (event_row in seq_len(nrow(terminal_events))) {
      tip_id <- terminal_events$phylo_tip_id[[event_row]]
      terminal_time_keys <- attr(
        terminal_events,
        "pst_v34_terminal_time_keys",
        exact = TRUE
      )
      terminal_key <- if (!is.null(terminal_time_keys)) {
        terminal_time_keys[[event_row]]
      } else {
        pst_v34_cpp_batch_time_keys_payload(
          terminal_events$terminal_time[[event_row]],
          time_tolerance
        )[[1L]]
      }
      post_terminal <- time_keys > terminal_key
      matrix_value[tip_id, post_terminal] <- NA
    }
    active[[matrix_name]] <- matrix_value
  }
  active
}
