# V34 batch-level debug trajectory assembly.

#' Truncate one TT family to a committed batch prefix.
#'
#' @param family Named TT family list.
#' @param step_id One-based committed TT row to retain through.
#' @param total_rows Number of rows in the complete TT trajectory.
#' @return TT family with every row-parallel matrix truncated to `step_id`.
#' @noRd
pst_v34_debug_prefix_tt_family <- function(family, step_id, total_rows) {
  # Preserve non-list metadata values unchanged; TT counter families are lists.
  if (!is.list(family)) {
    return(family)
  }

  out <- family
  # Visit every family component. Matrices with one row per batch are truncated;
  # lookup tables and other non-row-parallel values remain unchanged.
  for (component in names(out)) {
    value <- out[[component]]
    if ((is.matrix(value) || is.data.frame(value)) && nrow(value) == total_rows) {
      out[[component]] <- value[seq_len(step_id), , drop = FALSE]
    }
  }
  out
}

#' Build a TT object containing committed rows through one batch.
#'
#' @param tt Complete V34 TT object.
#' @param step_id One-based committed row to retain through.
#' @param time_tolerance Shared V34 resolution for committed-time comparisons.
#' @return V34 TT object with the same field contract and prefix rows.
#' @noRd
pst_v34_debug_prefix_tt <- function(tt, step_id, time_tolerance) {
  total_rows <- length(tt$time)
  cutoff_time <- tt$time[[step_id]]
  out <- tt
  out$time <- tt$time[seq_len(step_id)]

  # Truncate every public TT counter family while preserving its columns and
  # class attributes.
  for (family in c("ltt", "stt", "tt", "ctt", "ldif", "puniq", "tlen")) {
    out[[family]] <- pst_v34_debug_prefix_tt_family(
      tt[[family]],
      step_id,
      total_rows
    )
  }

  metadata <- tt$metadata
  # Event buffers are diagnostic row records; retain only events committed by
  # this snapshot time. Lookup/support metadata keeps its stable public shape.
  for (field in c(
    "change_records",
    "lineage_differentiating_records",
    "path_unique_records"
  )) {
    records <- metadata[[field]]
    if (is.data.frame(records) && "time" %in% names(records)) {
      metadata[[field]] <- records[
        records$time <= cutoff_time + time_tolerance,
        ,
        drop = FALSE
      ]
    }
  }
  out$metadata <- metadata
  class(out) <- class(tt)
  out
}

#' Truncate scenario matrices to one committed batch prefix.
#'
#' @param scenario_mats Complete V34 scenario-matrix object.
#' @param step_id One-based committed row/column to retain through.
#' @return Scenario matrices with batch-parallel columns and vectors truncated.
#' @noRd
pst_v34_debug_prefix_scenario_mats <- function(scenario_mats, step_id) {
  out <- scenario_mats
  # Matrix columns correspond one-to-one with committed TT times.
  for (field in c(
    "scenarios",
    "edges",
    "phylo_edge_ids",
    "scenario_edge_ids"
  )) {
    value <- out[[field]]
    if (is.matrix(value)) {
      out[[field]] <- value[, seq_len(step_id), drop = FALSE]
    }
  }
  # Time, event code, and event-tip records are parallel to matrix columns.
  for (field in c("time_vec", "event_vec", "event_tips")) {
    out[[field]] <- out[[field]][seq_len(step_id)]
  }
  out
}

#' Resolve selected debug snapshot rows from optional step/time selectors.
#'
#' @param time_vec Complete committed TT time vector.
#' @param steps Optional one-based step ids.
#' @param times Optional batch times matched at the shared V34 resolution.
#' @param time_tolerance Shared V34 resolution for committed-time comparisons.
#' @return Sorted unique one-based snapshot row ids.
#' @noRd
pst_v34_debug_selected_steps <- function(
    time_vec,
    steps = NULL,
    times = NULL,
    time_tolerance) {
  selected <- integer()
  # Explicit step selectors must reference existing committed rows.
  if (!is.null(steps)) {
    steps <- as.integer(steps)
    if (anyNA(steps) || any(steps < 1L | steps > length(time_vec))) {
      stop("debug_snapshot_steps contains an invalid committed step", call. = FALSE)
    }
    selected <- c(selected, steps)
  }
  # Time selectors map to every committed row within the shared time tolerance.
  if (!is.null(times)) {
    times <- as.numeric(times)
    if (anyNA(times) || any(!is.finite(times))) {
      stop("debug_snapshot_times must contain finite numeric times", call. = FALSE)
    }
    for (requested_time in times) {
      matched <- which(abs(time_vec - requested_time) <= time_tolerance)
      if (!length(matched)) {
        stop(
          sprintf("debug_snapshot_times did not match batch time %.17g", requested_time),
          call. = FALSE
        )
      }
      selected <- c(selected, matched)
    }
  }
  # With no selectors, snapshot every committed TT batch.
  if (is.null(steps) && is.null(times)) {
    selected <- seq_along(time_vec)
  }
  sort(unique(selected))
}

#' Attach selected batch-level debug trajectory snapshots.
#'
#' @param result Complete internally materialized V34 trajectory.
#' @param steps Optional one-based step selectors.
#' @param times Optional committed-time selectors.
#' @param time_tolerance Shared V34 resolution for committed-time comparisons.
#' @return `result` with trajectory, index, and status attributes.
#' @noRd
pst_v34_attach_debug_snapshots <- function(
    result,
    steps = NULL,
    times = NULL,
    time_tolerance) {
  traversal_summary <- attr(
    result,
    "pst_v34_cpp_traversal_summary",
    exact = TRUE
  )
  debug_step_metadata <- traversal_summary$debug_step_metadata
  # Snapshot mode requires the compact C++ event identity aligned with every
  # committed TT row. Refuse to infer mixed biology from compatibility matrices.
  if (!is.list(debug_step_metadata) ||
      length(debug_step_metadata) != length(result$tt$time)) {
    stop("V34 traversal omitted aligned debug step metadata", call. = FALSE)
  }
  selected_steps <- pst_v34_debug_selected_steps(
    result$tt$time,
    steps = steps,
    times = times,
    time_tolerance = time_tolerance
  )
  snapshots <- vector("list", length(selected_steps))
  index_rows <- vector("list", length(selected_steps))

  # Materialize each selected committed prefix as a full top-level trajectory;
  # C++ clips trees below while TT and scenario matrices retain exact prefixes.
  for (snapshot_id in seq_along(selected_steps)) {
    step_id <- selected_steps[[snapshot_id]]
    snapshot <- result
    snapshot$tt <- pst_v34_debug_prefix_tt(
      result$tt,
      step_id,
      time_tolerance
    )
    snapshot$scenario_mats <- pst_v34_debug_prefix_scenario_mats(
      result$scenario_mats,
      step_id
    )
    cutoff_time <- result$tt$time[[step_id]]
    # C++ clips every public tree at the exact committed time and represents
    # active frontier branches as deterministic provisional tips.
    snapshot$phylo <- pst_v34_cpp_debug_clip_tree_payload(
      result$phylo,
      cutoff_time,
      time_tolerance
    )
    snapshot$scenario <- pst_v34_cpp_debug_clip_tree_payload(
      result$scenario,
      cutoff_time,
      time_tolerance
    )
    snapshot$trans <- pst_v34_cpp_debug_clip_tree_payload(
      result$trans,
      cutoff_time,
      time_tolerance
    )
    step_metadata <- debug_step_metadata[[step_id]]
    required_metadata <- c(
      "event_types", "causal_phases", "event_ids", "phylogeny_edges",
      "scenario_ids", "terminal_ids", "batch_size"
    )
    # Every selected row must expose the complete C++ diagnostic contract;
    # missing fields would make failure localization partial or misleading.
    if (!is.list(step_metadata) ||
        !all(required_metadata %in% names(step_metadata))) {
      stop("V34 debug step metadata is incomplete", call. = FALSE)
    }
    event_types <- as.character(step_metadata$event_types)
    causal_phases <- as.character(step_metadata$causal_phases)
    event_ids <- as.integer(step_metadata$event_ids)
    event_phylo_edges <- as.integer(step_metadata$phylogeny_edges)
    event_scenario_ids <- as.integer(step_metadata$scenario_ids)
    terminal_ids <- as.integer(step_metadata$terminal_ids)
    batch_size <- as.integer(step_metadata$batch_size)[[1L]]
    # Scalar labels keep the index easy to scan in RStudio while list columns
    # below preserve every event and phase in a simultaneous transaction.
    event_type <- paste(event_types, collapse = "+")
    phase <- paste(unique(causal_phases), collapse = "+")

    # Capture metrics from the three clipped C++ trees as part of snapshot
    # materialization. Prefix SS then consumes only these compact parts and the
    # committed TT prefix, never the returned snapshot trees themselves.
    # Match the compatibility summary's state domain to labels represented by
    # the clipped phylogeny, excluding synthetic root-only lookup entries.
    snapshot_state_levels <- sort(unique(colnames(snapshot$phylo$mapped.edge)))
    snapshot$phylo <- pst_v34_cpp_state_only_tree_payload(
      snapshot$phylo,
      summary_suffix = "P",
      state_levels = snapshot_state_levels,
      include_edge_boundaries = FALSE,
      source_labels_are_state_only = TRUE
    )
    snapshot$scenario <- pst_v34_cpp_state_only_tree_payload(
      snapshot$scenario,
      summary_suffix = "S",
      state_levels = snapshot_state_levels,
      include_edge_boundaries = TRUE
    )
    snapshot$trans <- pst_v34_cpp_state_only_tree_payload(
      snapshot$trans,
      summary_suffix = "T",
      state_levels = snapshot_state_levels,
      include_edge_boundaries = TRUE
    )
    snapshot_tree_parts <- list(
      P = attr(snapshot$phylo, "pst_v34_tree_summary_parts"),
      S = attr(snapshot$scenario, "pst_v34_tree_summary_parts"),
      T = attr(snapshot$trans, "pst_v34_tree_summary_parts")
    )
    # Snapshot SS remains populated under the V30-compatible field contract.
    # Clipped metric parts and prefix TT ensure branch and event totals stop at
    # the same committed transaction boundary.
    snapshot$ss <- pst_v34_build_summary_stats_canonical(
      tt = snapshot$tt,
      tree_parts = snapshot_tree_parts,
      path_size_maps = snapshot$trans$path_size_maps
    )
    attr(snapshot, "pst_v34_tree_summary_parts") <- snapshot_tree_parts
    attr(snapshot, "pst_v34_debug_step_id") <- step_id
    attr(snapshot, "pst_v34_debug_time") <- cutoff_time
    attr(snapshot, "pst_v34_debug_event_type") <- event_types
    attr(snapshot, "pst_v34_debug_event_ids") <- event_ids
    attr(snapshot, "pst_v34_debug_batch_id") <- step_id
    attr(snapshot, "pst_v34_debug_phase") <- causal_phases
    snapshots[[snapshot_id]] <- snapshot

    index_rows[[snapshot_id]] <- data.frame(
      step_id = step_id,
      time = cutoff_time,
      event_type = event_type,
      batch_id = step_id,
      phase = phase,
      batch_size = batch_size,
      stringsAsFactors = FALSE
    )
    index_rows[[snapshot_id]]$event_types <- I(list(event_types))
    index_rows[[snapshot_id]]$event_ids <- I(list(event_ids))
    index_rows[[snapshot_id]]$causal_phases <- I(list(causal_phases))
    index_rows[[snapshot_id]]$phylogeny_edges <- I(list(event_phylo_edges))
    index_rows[[snapshot_id]]$scenario_ids <- I(list(event_scenario_ids))
    index_rows[[snapshot_id]]$terminal_ids <- I(list(terminal_ids))
  }

  step_index <- do.call(rbind, index_rows)
  rownames(step_index) <- NULL
  attr(result, "pst_v34_debug_trajectory_steps") <- snapshots
  attr(result, "pst_v34_debug_step_index") <- step_index
  attr(result, "pst_v34_debug_status") <- list(
    status = "complete",
    scope = "batch",
    saved_steps = selected_steps
  )
  result
}
