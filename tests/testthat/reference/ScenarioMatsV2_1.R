# ScenarioMats V2.1
#
# This is the readable, pure-R ScenarioMats implementation. Its central rule is that parallel
# history is reconstructed separately for each tip pair beginning where those two lineages
# separate. An inherited path from before separation can therefore never be mistaken for an
# independently repeated path.

.smv21_tip_index <- function(intervals, tip) {
    if (is.character(tip) && length(tip) == 1L) {
        tip <- match(tip, intervals$tip_labels)
    }
    if (!is.numeric(tip) || length(tip) != 1L || !is.finite(tip) || tip != floor(tip) || tip < 1L ||
        tip > length(intervals$tip_labels)) {
        stop("tip must identify one row in the prepared intervals", call. = FALSE)
    }
    as.integer(tip)
}

.smv21_view_name <- function(maximum_transition_count) {
    if (maximum_transition_count == 0L) {
        return("state_only")
    }
    paste0("using_up_to_", maximum_transition_count, "_", if (maximum_transition_count == 1L)
        "transition" else "transitions")
}

.smv21_transition_label <- function(transition_count) {
    paste0(transition_count, if (transition_count == 1L)
        "Transition" else "Transitions")
}

.smv21_depth_partition_label <- function(transition_count, maximum_transition_count) {
    if (transition_count == maximum_transition_count) {
        return(paste0(transition_count, "OrMoreTransitions"))
    }
    .smv21_transition_label(transition_count)
}

.smv21_format_state_path <- function(states) {
    escaped_states <- gsub("%", "%25", as.character(states), fixed = TRUE)
    escaped_states <- gsub("|", "%7C", escaped_states, fixed = TRUE)
    paste0("|", paste(escaped_states, collapse = "|"))
}

.smv21_project_state_path <- function(states, maximum_transition_count) {
    if (!length(states))
        return(character())
    if (maximum_transition_count == 0L)
        return(tail(states, 1L))
    tail(states, min(length(states), maximum_transition_count + 1L))
}

.smv21_project_existing_path <- function(path, maximum_transition_count) {
    if (is.na(path) || !nzchar(path))
        return(NA_character_)
    if (!startsWith(path, "|")) {
        stop("active phylo_paths entries must begin with '|'", call. = FALSE)
    }
    components <- strsplit(substring(path, 2L), "|", fixed = TRUE)[[1L]]
    if (maximum_transition_count == 0L) {
        component <- tail(components, 1L)
        component <- gsub("%7C", "|", component, fixed = TRUE)
        return(gsub("%25", "%", component, fixed = TRUE))
    }
    paste0("|", paste(tail(components, min(length(components), maximum_transition_count + 1L)), collapse = "|"))
}

.smv21_tree_times <- function(tree) {
    tip_count <- length(tree$tip.label)
    node_count <- tip_count + tree$Nnode
    root_node <- tip_count + 1L
    node_times <- rep(NA_real_, node_count)
    node_times[[root_node]] <- 0
    unresolved_edges <- seq_len(nrow(tree$edge))

    # Resolve edges only after their parent time is known.
    while (length(unresolved_edges)) {
        resolved_now <- unresolved_edges[!is.na(node_times[tree$edge[unresolved_edges, 1L]])]
        if (!length(resolved_now)) {
            stop("phylo edges cannot be traversed outward from the root", call. = FALSE)
        }
        for (edge_index in resolved_now) {
            parent <- tree$edge[edge_index, 1L]
            child <- tree$edge[edge_index, 2L]
            node_times[[child]] <- node_times[[parent]] + tree$edge.length[[edge_index]]
        }
        unresolved_edges <- setdiff(unresolved_edges, resolved_now)
    }
    node_times[seq_len(tip_count)]
}

.smv21_validate_depths <- function(maximum_transition_counts) {
    if (!is.numeric(maximum_transition_counts) || !length(maximum_transition_counts) || any(!is.finite(maximum_transition_counts)) ||
        any(maximum_transition_counts < 0) || any(maximum_transition_counts != floor(maximum_transition_counts)) ||
        any(maximum_transition_counts >= .Machine$integer.max) || anyDuplicated(maximum_transition_counts) ||
        !0L %in% maximum_transition_counts) {
        stop(paste("maximum_transition_counts must contain unique non-negative whole", "numbers including 0"),
            call. = FALSE)
    }
    as.integer(maximum_transition_counts)
}

.smv21_resolve_time_tolerance <- function(trajectory_obj, time_tolerance) {
    if (is.null(time_tolerance))
        time_tolerance <- trajectory_obj$time_tolerance
    if (is.null(time_tolerance)) {
        time_tolerance <- trajectory_obj$scenario_mats$time_tolerance
    }
    if (is.null(time_tolerance)) {
        time_tolerance <- attr(trajectory_obj, "scenario_mats_time_tolerance", exact = TRUE)
    }
    if (is.null(time_tolerance)) {
        legacy_names <- paste0("pst_v", c(35L, 34L, 33L, 32L), "_time_tolerance")
        for (legacy_name in legacy_names) {
            time_tolerance <- attr(trajectory_obj, legacy_name, exact = TRUE)
            if (!is.null(time_tolerance))
                break
        }
    }
    if (!is.numeric(time_tolerance) || length(time_tolerance) != 1L || !is.finite(time_tolerance) ||
        time_tolerance <= 0) {
        stop(paste("time_tolerance must be one finite positive number, supplied directly", "or stored on the trajectory object"),
            call. = FALSE)
    }
    as.numeric(time_tolerance)
}

# Validate the trajectory's rectangular scenario matrices and convert their point-in-time columns
# into half-open intervals with explicit durations.
# maximum_transition_counts = c(0,1,2,3)
prepare_scenario_intervals_v2_1 <- function(trajectory_obj, maximum_transition_counts = 0L, time_tolerance = NULL) {
    if (!is.list(trajectory_obj) || is.null(trajectory_obj$phylo) || is.null(trajectory_obj$scenario_mats) ||
        is.null(trajectory_obj$root_policy)) {
        stop("trajectory_obj must contain phylo, scenario_mats, and root_policy", call. = FALSE)
    }
    maximum_transition_counts <- .smv21_validate_depths(maximum_transition_counts)
    time_tolerance <- .smv21_resolve_time_tolerance(trajectory_obj, time_tolerance)

    tree <- trajectory_obj$phylo
    scenario_mats <- trajectory_obj$scenario_mats
    tip_labels <- as.character(tree$tip.label)
    if (length(tip_labels) < 2L || anyNA(tip_labels) || any(!nzchar(tip_labels)) || anyDuplicated(tip_labels)) {
        stop("phylo must contain at least two uniquely labelled tips", call. = FALSE)
    }
    if (is.null(tree$edge.length) || length(tree$edge.length) != nrow(tree$edge) || any(!is.finite(tree$edge.length)) ||
        any(tree$edge.length < 0)) {
        stop("phylo must contain finite non-negative edge lengths", call. = FALSE)
    }

    required_fields <- c("scenarios", "phylo_paths", "phylo_edge_ids", "scenario_edge_ids", "scenario_edge_step_ids")
    matrices <- scenario_mats[required_fields]
    if (any(vapply(matrices, is.null, logical(1))) || any(!vapply(matrices, is.matrix, logical(1)))) {
        stop("scenario_mats is missing required matrix fields", call. = FALSE)
    }
    matrix_dimensions <- dim(matrices[[1L]])
    if (matrix_dimensions[[1L]] != length(tip_labels) || matrix_dimensions[[2L]] < 2L || any(!vapply(matrices,
        function(value) identical(dim(value), matrix_dimensions), logical(1)))) {
        stop("scenario matrices must have identical tip-by-time dimensions", call. = FALSE)
    }
    for (field_name in required_fields) {
        row_labels <- rownames(matrices[[field_name]])
        if (!is.null(row_labels) && !identical(row_labels, tip_labels)) {
            stop("scenario matrix rows must exactly match phylo tip labels", call. = FALSE)
        }
        rownames(matrices[[field_name]]) <- tip_labels
    }

    time_vector <- scenario_mats$time_vec
    if (!is.numeric(time_vector) || length(time_vector) != matrix_dimensions[[2L]] || any(!is.finite(time_vector)) ||
        any(diff(time_vector) <= 0)) {
        stop("scenario_mats$time_vec must be finite, column-aligned, and increasing", call. = FALSE)
    }
    time_vector <- as.numeric(time_vector)
    tip_end_times <- .smv21_tree_times(tree)
    names(tip_end_times) <- tip_labels
    tip_end_columns <- integer(length(tip_labels))

    # Every tip must end on exactly one existing time column.
    for (tip_index in seq_along(tip_labels)) {
        matching_columns <- which(abs(time_vector - tip_end_times[[tip_index]]) <= time_tolerance)
        if (length(matching_columns) != 1L) {
            stop("every tip termination must match exactly one time_vec column", call. = FALSE)
        }
        tip_end_columns[[tip_index]] <- matching_columns
    }
    names(tip_end_columns) <- tip_labels

    interval_count <- length(time_vector) - 1L
    interval_columns <- seq_len(interval_count)
    interval_start <- head(time_vector, -1L)
    interval_end <- tail(time_vector, -1L)
    interval_duration <- interval_end - interval_start
    active <- matrix(FALSE, nrow = length(tip_labels), ncol = interval_count, dimnames = list(tip_labels,
        NULL))

    # A lineage is active through every interval ending at or before its tip.
    for (tip_index in seq_along(tip_labels)) {
        last_active_column <- tip_end_columns[[tip_index]] - 1L
        if (last_active_column > 0L) {
            active[tip_index, seq_len(last_active_column)] <- TRUE
        }
    }

    states <- matrix(as.character(matrices$scenarios[, interval_columns, drop = FALSE]), nrow = length(tip_labels),
        dimnames = list(tip_labels, colnames(matrices$scenarios)[interval_columns]))
    paths <- matrix(as.character(matrices$phylo_paths[, interval_columns, drop = FALSE]), nrow = length(tip_labels),
        dimnames = list(tip_labels, colnames(matrices$phylo_paths)[interval_columns]))
    if (any(active & (is.na(states) | states == "")) || any(active & (is.na(paths) | paths == ""))) {
        stop("active intervals require nonmissing states and paths", call. = FALSE)
    }
    for (field_name in c("phylo_edge_ids", "scenario_edge_ids", "scenario_edge_step_ids")) {
        values <- matrices[[field_name]][, interval_columns, drop = FALSE]
        if (!is.numeric(values) || any(active & (is.na(values) | values < 1 | values != floor(values)))) {
            stop("active edge and step ids must be positive whole numbers", call. = FALSE)
        }
    }

    available_transition_count <- matrix(NA_integer_, nrow = nrow(paths), ncol = ncol(paths), dimnames = dimnames(paths))
    for (cell_index in seq_along(paths)) {
        if (!is.na(paths[[cell_index]]) && nzchar(paths[[cell_index]])) {
            components <- strsplit(substring(paths[[cell_index]], 2L), "|", fixed = TRUE)[[1L]]
            available_transition_count[[cell_index]] <- length(components) - 1L
        }
    }

    paths_by_depth <- list()
    for (maximum_transition_count in maximum_transition_counts) {
        projected_paths <- matrix(NA_character_, nrow = nrow(paths), ncol = ncol(paths), dimnames = dimnames(paths))
        for (cell_index in seq_along(paths)) {
            projected_paths[[cell_index]] <- .smv21_project_existing_path(paths[[cell_index]], maximum_transition_count)
        }
        paths_by_depth[[.smv21_view_name(maximum_transition_count)]] <- projected_paths
    }

    list(tip_labels = tip_labels, time = data.frame(interval_id = seq_len(interval_count), source_column = interval_columns,
        start = interval_start, end = interval_end, duration = interval_duration), states = states, paths = paths,
        paths_using_up_to_transition_count = paths_by_depth, available_transition_count = available_transition_count,
        phylo_edge_ids = matrices$phylo_edge_ids[, interval_columns, drop = FALSE], scenario_edge_ids = matrices$scenario_edge_ids[,
            interval_columns, drop = FALSE], scenario_edge_step_ids = matrices$scenario_edge_step_ids[,
            interval_columns, drop = FALSE], active = active, duration = interval_duration, tip_end_times = tip_end_times,
        tip_end_columns = tip_end_columns, root_anchor = trajectory_obj$root_policy$anchor, root_policy = trajectory_obj$root_policy,
        maximum_transition_counts = maximum_transition_counts, state_levels = sort(unique(states[active])),
        time_tolerance = time_tolerance)
}

# Build the readable interval ledger for one pair. The post-separation paths in this ledger contain
# raw biological state vectors as list columns and formatted labels beside them for direct
# inspection.
build_pair_history_v2_1 <- function(intervals, first_tip, second_tip) {
    first_tip_index <- .smv21_tip_index(intervals, first_tip)
    second_tip_index <- .smv21_tip_index(intervals, second_tip)
    if (first_tip_index == second_tip_index) {
        stop("the two tips must be different", call. = FALSE)
    }

    first_lineage_active <- intervals$active[first_tip_index, ]
    second_lineage_active <- intervals$active[second_tip_index, ]
    comparable <- first_lineage_active & second_lineage_active
    comparable_columns <- which(comparable)
    total_comparable_time <- sum(intervals$duration[comparable_columns])
    if (!length(comparable_columns) || total_comparable_time <= intervals$time_tolerance) {
        stop("pair has zero total comparable time", call. = FALSE)
    }
    history_columns <- which(first_lineage_active | second_lineage_active)
    first_total_history_time <- sum(intervals$duration[first_lineage_active])
    second_total_history_time <- sum(intervals$duration[second_lineage_active])

    same_phylogeny_edge <- intervals$phylo_edge_ids[first_tip_index, comparable_columns] == intervals$phylo_edge_ids[second_tip_index,
        comparable_columns]
    shared_positions <- which(same_phylogeny_edge)
    if (length(shared_positions)) {
        last_shared_position <- max(shared_positions)
        if (last_shared_position < length(same_phylogeny_edge) && any(same_phylogeny_edge[seq.int(last_shared_position +
            1L, length(same_phylogeny_edge))])) {
            stop("a tip pair cannot return to a shared phylogeny edge", call. = FALSE)
        }
    }

    root_is_synthetic <- isTRUE(intervals$root_policy$synthetic) || identical(intervals$root_policy$provenance,
        "synthetic-root")
    last_shared_position <- if (length(shared_positions)) {
        max(shared_positions)
    } else {
        0L
    }
    pair_mrca_state <- if (last_shared_position > 0L) {
        shared_column <- comparable_columns[[last_shared_position]]
        as.character(intervals$states[first_tip_index, shared_column])
    } else if (!root_is_synthetic) {
        as.character(intervals$root_anchor)
    } else {
        NA_character_
    }
    if (!root_is_synthetic && (length(pair_mrca_state) != 1L || is.na(pair_mrca_state) || !nzchar(pair_mrca_state))) {
        stop("a biological MRCA requires one nonmissing state", call. = FALSE)
    }

    first_post_separation_states <- character()
    second_post_separation_states <- character()
    records <- vector("list", length(history_columns))
    first_state_histories <- vector("list", length(history_columns))
    second_state_histories <- vector("list", length(history_columns))

    # Walk through the union of both histories. Rows after the shorter tip ends
    # remain available to asynchronous calculations, but have no synchronous
    # class because the two lineages are no longer simultaneously active.
    for (record_index in seq_along(history_columns)) {
        column <- history_columns[[record_index]]
        first_is_active <- first_lineage_active[[column]]
        second_is_active <- second_lineage_active[[column]]
        both_are_active <- first_is_active && second_is_active
        first_state <- if (first_is_active) as.character(intervals$states[first_tip_index, column]) else NA_character_
        second_state <- if (second_is_active) as.character(intervals$states[second_tip_index, column]) else NA_character_
        first_phylogeny_edge <- if (first_is_active) intervals$phylo_edge_ids[first_tip_index, column] else NA_integer_
        second_phylogeny_edge <- if (second_is_active) intervals$phylo_edge_ids[second_tip_index, column] else NA_integer_
        interval_is_shared <- both_are_active && first_phylogeny_edge == second_phylogeny_edge

        requested_depth_used <- NA_integer_
        actual_transition_count_used <- NA_integer_
        matching_path <- NA_character_
        phase <- "Asynchronous_only"
        synchronous_class <- NA_character_
        asynchronous_class_i <- "Inactive"
        asynchronous_class_j <- "Inactive"
        first_path_label <- NA_character_
        second_path_label <- NA_character_
        first_available_count <- NA_integer_
        second_available_count <- NA_integer_

        if (interval_is_shared) {
            phase <- "Shared"
            synchronous_class <- "Shared"
            asynchronous_class_i <- "Shared"
            asynchronous_class_j <- "Shared"
            first_path_label <- NA_character_
            second_path_label <- NA_character_
            first_available_count <- NA_integer_
            second_available_count <- NA_integer_
        } else {
            # Initialize and extend each active lineage separately. This allows
            # the longer lineage to retain changes made after the shorter tip
            # has ended.
            if (first_is_active) {
                if (!length(first_post_separation_states)) {
                  if (root_is_synthetic && last_shared_position == 0L) {
                    first_post_separation_states <- first_state
                  } else {
                    first_post_separation_states <- pair_mrca_state
                  }
                }
                if (tail(first_post_separation_states, 1L) != first_state) {
                  first_post_separation_states <- c(first_post_separation_states, first_state)
                }
                first_available_count <- length(first_post_separation_states) - 1L
                first_path_label <- .smv21_format_state_path(first_post_separation_states)
                first_lineage_is_conserved <- !root_is_synthetic && first_available_count == 0L &&
                  identical(first_state, pair_mrca_state)
                asynchronous_class_i <- if (first_lineage_is_conserved) "Conserved" else "Independent"
            } else {
                first_lineage_is_conserved <- FALSE
            }

            if (second_is_active) {
                if (!length(second_post_separation_states)) {
                  if (root_is_synthetic && last_shared_position == 0L) {
                    second_post_separation_states <- second_state
                  } else {
                    second_post_separation_states <- pair_mrca_state
                  }
                }
                if (tail(second_post_separation_states, 1L) != second_state) {
                  second_post_separation_states <- c(second_post_separation_states, second_state)
                }
                second_available_count <- length(second_post_separation_states) - 1L
                second_path_label <- .smv21_format_state_path(second_post_separation_states)
                second_lineage_is_conserved <- !root_is_synthetic && second_available_count == 0L &&
                  identical(second_state, pair_mrca_state)
                asynchronous_class_j <- if (second_lineage_is_conserved) "Conserved" else "Independent"
            } else {
                second_lineage_is_conserved <- FALSE
            }

            # Only rows where both lineages exist receive a synchronous class.
            if (both_are_active) {
                interval_is_conserved <- first_lineage_is_conserved && second_lineage_is_conserved
                if (interval_is_conserved) {
                  phase <- "Conserved"
                  synchronous_class <- "Conserved"
                } else {
                  phase <- "Independent"
                  if (!identical(first_state, second_state)) {
                    synchronous_class <- "Divergent"
                  } else if (xor(first_lineage_is_conserved, second_lineage_is_conserved)) {
                    phase <- "Conserved_Homoplasy"
                    synchronous_class <- "Conserved_Homoplasy"
                  } else {
                    deepest_actual_match <- 0L
                    positive_depths <- intervals$maximum_transition_counts[intervals$maximum_transition_counts >
                      0L]

                    # Compare only suffixes accumulated after this pair separated.
                    for (requested_depth in positive_depths) {
                      first_projection <- .smv21_project_state_path(first_post_separation_states, requested_depth)
                      second_projection <- .smv21_project_state_path(second_post_separation_states, requested_depth)
                      actual_match <- min(requested_depth, first_available_count, second_available_count)
                      if (actual_match > 0L && identical(first_projection, second_projection) && actual_match >
                        deepest_actual_match) {
                        deepest_actual_match <- actual_match
                        requested_depth_used <- requested_depth
                        matching_path <- .smv21_format_state_path(first_projection)
                      }
                    }

                    if (deepest_actual_match > 0L) {
                      synchronous_class <- "Parallel"
                      actual_transition_count_used <- deepest_actual_match
                    } else {
                      synchronous_class <- "Convergent"
                      requested_depth_used <- NA_integer_
                      matching_path <- NA_character_
                    }
                  }
                }
            }
        }

        first_state_histories[[record_index]] <- first_post_separation_states
        second_state_histories[[record_index]] <- second_post_separation_states
        records[[record_index]] <- data.frame(interval_id = intervals$time$interval_id[[column]], start = intervals$time$start[[column]],
            end = intervals$time$end[[column]], duration = intervals$duration[[column]], state_i = first_state,
            state_j = second_state, first_lineage_active = first_is_active, second_lineage_active = second_is_active,
            path_i = if (first_is_active) as.character(intervals$paths[first_tip_index, column]) else NA_character_,
            path_j = if (second_is_active) as.character(intervals$paths[second_tip_index, column]) else NA_character_,
            post_separation_path_i = first_path_label,
            post_separation_path_j = second_path_label, phylo_edge_i = first_phylogeny_edge, phylo_edge_j = second_phylogeny_edge,
            scenario_edge_i = if (first_is_active) intervals$scenario_edge_ids[first_tip_index, column] else NA_integer_,
            scenario_edge_j = if (second_is_active) intervals$scenario_edge_ids[second_tip_index, column] else NA_integer_,
            scenario_step_i = if (first_is_active) intervals$scenario_edge_step_ids[first_tip_index, column] else NA_integer_,
            scenario_step_j = if (second_is_active) intervals$scenario_edge_step_ids[second_tip_index, column] else NA_integer_,
            pair_mrca_state = pair_mrca_state,
            root_is_synthetic = root_is_synthetic, phase = phase, sync_class = synchronous_class,
            asynchronous_class_i = asynchronous_class_i, asynchronous_class_j = asynchronous_class_j,
            available_transition_count_i = first_available_count,
            available_transition_count_j = second_available_count, requested_maximum_transition_count = requested_depth_used,
            actual_transition_count_used = actual_transition_count_used, matching_path = matching_path)
    }

    interval_records <- do.call(rbind, records)
    interval_records$post_separation_states_i <- I(first_state_histories)
    interval_records$post_separation_states_j <- I(second_state_histories)
    list(row_i = first_tip_index, row_j = second_tip_index, tip_i = intervals$tip_labels[[first_tip_index]],
        tip_j = intervals$tip_labels[[second_tip_index]], pair_mrca_state = pair_mrca_state, root_is_synthetic = root_is_synthetic,
        total_comparable_time = total_comparable_time, first_total_history_time = first_total_history_time,
        second_total_history_time = second_total_history_time, intervals = interval_records)
}

.smv21_component_definitions <- function(states, maximum_transition_counts) {
    synchronous_names <- character()
    metadata_rows <- list()
    metadata_row <- 0L

    add_metadata <- function(matrix_name, state, phase, kind, requested_depth, actual_depth, view) {
        metadata_row <<- metadata_row + 1L
        metadata_rows[[metadata_row]] <<- data.frame(matrix_name = matrix_name, state = state, path = NA_character_,
            phase = phase, kind = kind, requested_maximum_transition_count = requested_depth, actual_transition_count_used = actual_depth,
            view = view)
    }

    for (state in states) {
        synchronous_names <- c(synchronous_names, paste0(state, "_Shared"), paste0(state, "_Conserved"),
            paste0(state, "_Conserved_Homoplasy"), paste0(state, "_Convergent"))
        add_metadata(paste0(state, "_Shared"), state, "Shared", "state", NA_integer_, NA_integer_, "sync")
        add_metadata(paste0(state, "_Conserved"), state, "Conserved", "state", NA_integer_, NA_integer_,
            "sync")
        add_metadata(paste0(state, "_Conserved_Homoplasy"), state, "Conserved_Homoplasy", "state",
            NA_integer_, NA_integer_, "sync")
        add_metadata(paste0(state, "_Convergent"), state, "Independent", "state", NA_integer_, 0L, "sync")
    }
    largest_depth <- max(maximum_transition_counts)
    if (largest_depth > 0L) {
        for (actual_depth in seq_len(largest_depth)) {
            for (state in states) {
                matrix_name <- paste0(state, "_Parallel_", .smv21_transition_label(actual_depth))
                synchronous_names <- c(synchronous_names, matrix_name)
                add_metadata(matrix_name, state, "Independent", "parallel", NA_integer_, actual_depth,
                  "sync")
            }
        }
    }
    synchronous_names <- c(synchronous_names, "Divergent")
    add_metadata("Divergent", NA_character_, "Independent", "residual", NA_integer_, NA_integer_, "sync")

    asynchronous_names <- list()
    for (requested_depth in maximum_transition_counts) {
        view_name <- .smv21_view_name(requested_depth)
        view_names <- character()
        for (state in states) {
            view_names <- c(view_names, paste0(state, "_Shared"), paste0(state, "_Conserved"),
                paste0(state, "_Conserved_Homoplasy"))
            add_metadata(paste0(state, "_Shared"), state, "Shared", "state", requested_depth, NA_integer_,
                view_name)
            add_metadata(paste0(state, "_Conserved"), state, "Conserved", "state", requested_depth, NA_integer_,
                view_name)
            add_metadata(paste0(state, "_Conserved_Homoplasy"), state, "Conserved_Homoplasy", "state",
                requested_depth, NA_integer_, view_name)
            if (requested_depth == 0L) {
                view_names <- c(view_names, paste0(state, "_Convergent"))
                add_metadata(paste0(state, "_Convergent"), state, "Independent", "convergent", 0L, 0L,
                  view_name)
            } else {
                for (actual_depth in seq_len(requested_depth)) {
                  parallel_name <- paste0(state, "_Parallel_", .smv21_depth_partition_label(actual_depth,
                    requested_depth))
                  view_names <- c(view_names, parallel_name)
                  add_metadata(parallel_name, state, "Independent", "parallel", requested_depth, actual_depth,
                    view_name)
                }
                view_names <- c(view_names, paste0(state, "_Convergent"))
                add_metadata(paste0(state, "_Convergent"), state, "Independent", "convergent", requested_depth,
                  0L, view_name)
            }
        }
        view_names <- c(view_names, "Divergent")
        add_metadata("Divergent", NA_character_, "Independent", "residual", requested_depth, NA_integer_,
            view_name)
        asynchronous_names[[view_name]] <- view_names
    }

    list(sync = synchronous_names, async = asynchronous_names, metadata = do.call(rbind, metadata_rows))
}

.smv21_similarity <- function(first_duration, second_duration, first_total_history_time,
    second_total_history_time, metric) {
    if (metric == "bhattacharyya") {
        sqrt(first_duration * second_duration)/sqrt(first_total_history_time * second_total_history_time)
    } else {
        min(first_duration, second_duration)/min(first_total_history_time, second_total_history_time)
    }
}

.smv21_guard_values <- function(values, tolerance) {
    if (any(!is.finite(values)) || any(values < -tolerance) || any(values > 1 + tolerance)) {
        stop("similarity component fell outside [0, 1]", call. = FALSE)
    }
    values[values < 0] <- 0
    total <- sum(values)
    if (abs(total - 1) > tolerance) {
        stop("similarity components do not sum to one", call. = FALSE)
    }
    if (total != 1)
        values <- values/total
    values
}

# Convert one pair's readable history ledger into synchronous interval mass and asynchronous
# duration-profile overlap.
calculate_pair_similarity_v2_1 <- function(intervals, first_tip, second_tip, async_metric = c("bhattacharyya",
    "minimum")) {
    if (missing(async_metric))
        async_metric <- "bhattacharyya"
    if (!is.character(async_metric) || length(async_metric) != 1L || is.na(async_metric) || !async_metric %in%
        c("bhattacharyya", "minimum")) {
        stop("async_metric must be either 'bhattacharyya' or 'minimum'", call. = FALSE)
    }

    pair_history <- build_pair_history_v2_1(intervals, first_tip, second_tip)
    records <- pair_history$intervals
    total_comparable_time <- pair_history$total_comparable_time
    first_total_history_time <- pair_history$first_total_history_time
    second_total_history_time <- pair_history$second_total_history_time
    component_definitions <- .smv21_component_definitions(intervals$state_levels, intervals$maximum_transition_counts)
    synchronous_values <- setNames(numeric(length(component_definitions$sync)), component_definitions$sync)
    synchronous_parallel_rows <- list()
    synchronous_parallel_row <- 0L

    # Every interval contributes to exactly one synchronous component.
    for (record_index in seq_len(nrow(records))) {
        record <- records[record_index, ]
        if (is.na(record$sync_class))
            next
        component_name <- switch(record$sync_class, Shared = paste0(record$state_i, "_Shared"), Conserved = paste0(record$state_i,
            "_Conserved"), Conserved_Homoplasy = paste0(record$state_i, "_Conserved_Homoplasy"),
            Convergent = paste0(record$state_i, "_Convergent"), Parallel = paste0(record$state_i,
            "_Parallel_", .smv21_transition_label(record$actual_transition_count_used)), Divergent = "Divergent")
        synchronous_values[[component_name]] <- synchronous_values[[component_name]] + record$duration/total_comparable_time

        if (record$sync_class == "Parallel") {
            synchronous_parallel_row <- synchronous_parallel_row + 1L
            synchronous_parallel_rows[[synchronous_parallel_row]] <- data.frame(requested_maximum_transition_count = record$requested_maximum_transition_count,
                actual_transition_count_used = record$actual_transition_count_used, path = record$matching_path,
                terminal_state = record$state_i, duration = record$duration)
        }
    }
    synchronous_values <- .smv21_guard_values(synchronous_values, intervals$time_tolerance)

    shared_first <- shared_second <- setNames(numeric(length(intervals$state_levels)), intervals$state_levels)
    conserved_first <- conserved_second <- shared_first
    independent_first <- independent_second <- shared_first

    # Build the state-duration profiles used by every asynchronous depth.
    for (record_index in seq_len(nrow(records))) {
        record <- records[record_index, ]
        if (record$asynchronous_class_i == "Shared") {
            shared_first[[record$state_i]] <- shared_first[[record$state_i]] + record$duration
        } else if (record$asynchronous_class_i == "Conserved") {
            conserved_first[[record$state_i]] <- conserved_first[[record$state_i]] + record$duration
        } else if (record$asynchronous_class_i == "Independent") {
            independent_first[[record$state_i]] <- independent_first[[record$state_i]] + record$duration
        }

        if (record$asynchronous_class_j == "Shared") {
            shared_second[[record$state_j]] <- shared_second[[record$state_j]] + record$duration
        } else if (record$asynchronous_class_j == "Conserved") {
            conserved_second[[record$state_j]] <- conserved_second[[record$state_j]] + record$duration
        } else if (record$asynchronous_class_j == "Independent") {
            independent_second[[record$state_j]] <- independent_second[[record$state_j]] + record$duration
        }
    }

    shared_similarity <- conserved_similarity <- independent_similarity <- post_separation_state_similarity <-
        conserved_homoplasy_similarity <- setNames(numeric(length(intervals$state_levels)),
        intervals$state_levels)
    for (state in intervals$state_levels) {
        shared_similarity[[state]] <- .smv21_similarity(shared_first[[state]], shared_second[[state]],
            first_total_history_time, second_total_history_time, async_metric)
        conserved_similarity[[state]] <- .smv21_similarity(conserved_first[[state]], conserved_second[[state]],
            first_total_history_time, second_total_history_time, async_metric)
        independent_similarity[[state]] <- .smv21_similarity(independent_first[[state]], independent_second[[state]],
            first_total_history_time, second_total_history_time, async_metric)
        post_separation_state_similarity[[state]] <- .smv21_similarity(
            conserved_first[[state]] + independent_first[[state]],
            conserved_second[[state]] + independent_second[[state]],
            first_total_history_time, second_total_history_time, async_metric)
        conserved_homoplasy_similarity[[state]] <- post_separation_state_similarity[[state]] -
            conserved_similarity[[state]] - independent_similarity[[state]]
        if (conserved_homoplasy_similarity[[state]] < -intervals$time_tolerance) {
            stop("conserved-homoplasy similarity became negative", call. = FALSE)
        }
        conserved_homoplasy_similarity[[state]] <- max(0, conserved_homoplasy_similarity[[state]])
    }
    asynchronous_divergent <- max(0, 1 - sum(shared_similarity) - sum(post_separation_state_similarity))

    asynchronous_values <- list()
    asynchronous_independent_rows <- list()
    asynchronous_parallel_rows <- list()
    independent_row <- 0L
    parallel_row <- 0L

    # This function already represents one tip pair. Each pass builds that
    # pair's two path-duration profiles at one requested depth.
    for (requested_depth in intervals$maximum_transition_counts) {
        view_name <- .smv21_view_name(requested_depth)
        view_values <- setNames(numeric(length(component_definitions$async[[view_name]])), component_definitions$async[[view_name]])

        # Shared and conserved similarity do not depend on path depth, so copy
        # the previously calculated state values into every depth view.
        for (state in intervals$state_levels) {
            view_values[[paste0(state, "_Shared")]] <- shared_similarity[[state]]
            view_values[[paste0(state, "_Conserved")]] <- conserved_similarity[[state]]
            view_values[[paste0(state, "_Conserved_Homoplasy")]] <- conserved_homoplasy_similarity[[state]]
        }

        # At depth zero there is no path comparison. All same-state similarity
        # that is neither shared nor conserved is provisionally convergent.
        if (requested_depth == 0L) {
            for (state in intervals$state_levels) {
                view_values[[paste0(state, "_Convergent")]] <- independent_similarity[[state]]
            }
        } else {
            # These named vectors are duration profiles: each path label points
            # to the total time that one tip spent following that path.
            first_path_durations <- numeric()
            second_path_durations <- numeric()
            path_states <- character()
            path_transition_counts <- integer()

            # Each pass handles one time interval for this tip pair. Shared and
            # conserved intervals are already accounted for above; the remaining
            # intervals are projected to this depth and added to each tip's
            # path-duration profile. Similarity is not calculated in this loop.
            for (record_index in seq_len(nrow(records))) {
                record <- records[record_index, ]
                first_lineage_is_independent <- record$asynchronous_class_i == "Independent"
                second_lineage_is_independent <- record$asynchronous_class_j == "Independent"
                if (!first_lineage_is_independent && !second_lineage_is_independent)
                  next

                if (first_lineage_is_independent) {
                  first_projection <- .smv21_project_state_path(records$post_separation_states_i[[record_index]],
                    requested_depth)
                  first_path <- .smv21_format_state_path(first_projection)
                  first_count <- length(first_projection) - 1L

                  if (!first_path %in% names(first_path_durations)) {
                    first_path_durations[[first_path]] <- 0
                  }
                  first_path_durations[[first_path]] <- first_path_durations[[first_path]] + record$duration
                  path_states[[first_path]] <- record$state_i
                  path_transition_counts[[first_path]] <- first_count

                  # Expose this unaggregated lineage contribution for complete-mode auditing.
                  independent_row <- independent_row + 1L
                  asynchronous_independent_rows[[independent_row]] <- data.frame(lineage = "i", tip = pair_history$tip_i,
                    requested_maximum_transition_count = requested_depth, available_transition_count = record$available_transition_count_i,
                    actual_transition_count_used = first_count, path = first_path, terminal_state = record$state_i,
                    duration = record$duration, proportion = record$duration/first_total_history_time)
                }

                if (second_lineage_is_independent) {
                  second_projection <- .smv21_project_state_path(records$post_separation_states_j[[record_index]],
                    requested_depth)
                  second_path <- .smv21_format_state_path(second_projection)
                  second_count <- length(second_projection) - 1L

                  if (!second_path %in% names(second_path_durations)) {
                    second_path_durations[[second_path]] <- 0
                  }
                  second_path_durations[[second_path]] <- second_path_durations[[second_path]] + record$duration
                  path_states[[second_path]] <- record$state_j
                  path_transition_counts[[second_path]] <- second_count

                  # Expose this unaggregated lineage contribution for complete-mode auditing.
                  independent_row <- independent_row + 1L
                  asynchronous_independent_rows[[independent_row]] <- data.frame(lineage = "j", tip = pair_history$tip_j,
                    requested_maximum_transition_count = requested_depth, available_transition_count = record$available_transition_count_j,
                    actual_transition_count_used = second_count, path = second_path, terminal_state = record$state_j,
                    duration = record$duration, proportion = record$duration/second_total_history_time)
                }
            }

            parallel_by_state <- setNames(numeric(length(intervals$state_levels)), intervals$state_levels)

            # Only now are the two tips compared. A path can contribute parallel
            # similarity only when its label occurs in both duration profiles.
            common_paths <- intersect(names(first_path_durations), names(second_path_durations))

            # Matching zero-transition paths remain convergent.
            for (path in common_paths) {
                transition_count <- path_transition_counts[[path]]
                if (transition_count == 0L)
                  next
                similarity <- .smv21_similarity(first_path_durations[[path]], second_path_durations[[path]],
                  first_total_history_time, second_total_history_time, async_metric)
                if (similarity == 0)
                  next
                terminal_state <- path_states[[path]]
                component_name <- paste0(terminal_state, "_Parallel_", .smv21_depth_partition_label(transition_count,
                  requested_depth))
                view_values[[component_name]] <- view_values[[component_name]] + similarity
                parallel_by_state[[terminal_state]] <- parallel_by_state[[terminal_state]] + similarity
                parallel_row <- parallel_row + 1L
                asynchronous_parallel_rows[[parallel_row]] <- data.frame(requested_maximum_transition_count = requested_depth,
                  actual_transition_count_used = transition_count, path = path, terminal_state = terminal_state,
                  similarity = similarity, metric = async_metric)
            }

            for (state in intervals$state_levels) {
                remaining_convergence <- independent_similarity[[state]] - parallel_by_state[[state]]
                if (remaining_convergence < -intervals$time_tolerance) {
                  stop("parallel path similarity exceeds same-state similarity", call. = FALSE)
                }
                view_values[[paste0(state, "_Convergent")]] <- max(0, remaining_convergence)
            }
        }
        view_values[["Divergent"]] <- asynchronous_divergent
        asynchronous_values[[view_name]] <- .smv21_guard_values(view_values, intervals$time_tolerance)
    }

    empty_sync_parallel <- data.frame(requested_maximum_transition_count = integer(), actual_transition_count_used = integer(),
        path = character(), terminal_state = character(), duration = numeric())
    empty_independent <- data.frame(lineage = character(), tip = character(), requested_maximum_transition_count = integer(),
        available_transition_count = integer(), actual_transition_count_used = integer(), path = character(),
        terminal_state = character(), duration = numeric(), proportion = numeric())
    empty_async_parallel <- data.frame(requested_maximum_transition_count = integer(), actual_transition_count_used = integer(),
        path = character(), terminal_state = character(), similarity = numeric(), metric = character())

    list(total_comparable_time = total_comparable_time, first_total_history_time = first_total_history_time,
        second_total_history_time = second_total_history_time, sync = synchronous_values, async = asynchronous_values,
        details = list(interval_classification = records, sync_parallel_time_by_path = if (length(synchronous_parallel_rows)) {
            do.call(rbind, synchronous_parallel_rows)
        } else {
            empty_sync_parallel
        }, async_independent_time_by_path = if (length(asynchronous_independent_rows)) {
            do.call(rbind, asynchronous_independent_rows)
        } else {
            empty_independent
        }, async_parallel_similarity_by_path = if (length(asynchronous_parallel_rows)) {
            do.call(rbind, asynchronous_parallel_rows)
        } else {
            empty_async_parallel
        }))
}

.smv21_summarize_family <- function(matrices, available_similarity_by_pair = NULL) {
    if (!length(matrices)) {
        return(list(totals = numeric(), tree_wide_proportions = numeric(),
            means = numeric(), total_across_matrices = 0))
    }
    upper_triangle <- upper.tri(matrices[[1L]])
    if (!is.null(available_similarity_by_pair)) {
        if (!identical(dim(available_similarity_by_pair), dim(matrices[[1L]])) ||
            any(!is.finite(available_similarity_by_pair)) || any(available_similarity_by_pair < 0)) {
            stop("available-similarity weights must match the component matrices", call. = FALSE)
        }
        available_similarity_total <- sum(available_similarity_by_pair[upper_triangle])
        if (available_similarity_total <= 0) {
            stop("available-similarity weights must have a positive upper-triangle sum", call. = FALSE)
        }
        totals <- vapply(matrices, function(matrix_value) sum(matrix_value[upper_triangle] *
            available_similarity_by_pair[upper_triangle]), numeric(1))
        total_across_matrices <- sum(totals)
        tree_wide_proportions <- if (total_across_matrices > 0) {
            totals/total_across_matrices
        } else {
            totals
        }
        return(list(totals = totals,
            tree_wide_proportions = tree_wide_proportions,
            means = tree_wide_proportions,
            total_across_matrices = total_across_matrices,
            available_similarity_total = available_similarity_total))
    }
    totals <- vapply(matrices, function(matrix_value) sum(matrix_value[upper_triangle]), numeric(1))
    total_across_matrices <- sum(totals)
    tree_wide_proportions <- if (total_across_matrices > 0) {
        totals/total_across_matrices
    } else {
        totals
    }
    list(totals = totals,
        tree_wide_proportions = tree_wide_proportions,
        means = tree_wide_proportions,
        total_across_matrices = total_across_matrices)
}

.smv21_add_depth_reports <- function(result) {
    metadata <- result$summaries$components
    asynchronous_summary_weights <- if (isTRUE(result$summaries$weighting$available_similarity_weighting_applied)) {
        result$summaries$weighting$asynchronous_available_similarity_by_pair
    } else {
        NULL
    }
    matrix_reports <- list()
    summary_reports <- list()

    # Present the state-specific matrices and their totals in the same additional report used by
    # V3.
    for (view_name in names(result$matrices$async)) {
        rows <- metadata[metadata$view == view_name & metadata$kind %in% c("convergent", "parallel"),
            , drop = FALSE]
        requested_depth <- unique(rows$requested_maximum_transition_count)
        requested_depth <- requested_depth[!is.na(requested_depth)]
        if (length(requested_depth) != 1L || requested_depth == 0L)
            next

        state_matrices <- result$matrices$async[[view_name]][rows$matrix_name]
        total_matrices <- list()
        for (row_index in seq_len(nrow(rows))) {
            total_name <- if (rows$kind[[row_index]] == "convergent") {
                "Convergent"
            } else if (rows$actual_transition_count_used[[row_index]] == requested_depth) {
                paste0("Parallel_", requested_depth, "OrMoreTransitions")
            } else {
                paste0("Parallel_", .smv21_transition_label(rows$actual_transition_count_used[[row_index]]))
            }
            if (is.null(total_matrices[[total_name]])) {
                total_matrices[[total_name]] <- state_matrices[[row_index]]
            } else {
                total_matrices[[total_name]] <- total_matrices[[total_name]] + state_matrices[[row_index]]
            }
        }
        matrix_reports[[view_name]] <- list(state = state_matrices, total = total_matrices)
        summary_reports[[view_name]] <- list(state = .smv21_summarize_family(state_matrices,
            asynchronous_summary_weights), total = .smv21_summarize_family(total_matrices,
            asynchronous_summary_weights))
    }
    result$matrices$async_convergence_by_transition_depth <- matrix_reports
    result$summaries$async_convergence_by_transition_depth <- summary_reports
    result
}

# Calculate every unordered tip pair and assemble the V3-compatible public matrix and summary
# structure.
#async_metric = "bhattacharyya"
# maximum_transition_counts = c(0,1,2,3)
# return_mode = "complete"
# time_tolerance = NULL
trajectory_similarity_v2_1 <- function(trajectory_obj, maximum_transition_counts = 0L, async_metric = c("bhattacharyya",
    "minimum"), return_mode = c("matrices_and_summaries", "complete"), time_tolerance = NULL) {
    if (missing(trajectory_obj)) {
        stop("argument \"trajectory_obj\" is missing", call. = FALSE)
    }
    maximum_transition_counts <- .smv21_validate_depths(maximum_transition_counts)
    if (missing(async_metric))
        async_metric <- "bhattacharyya"
    if (!is.character(async_metric) || length(async_metric) != 1L || is.na(async_metric) || !async_metric %in%
        c("bhattacharyya", "minimum")) {
        stop("async_metric must be either 'bhattacharyya' or 'minimum'", call. = FALSE)
    }
    if (missing(return_mode))
        return_mode <- "matrices_and_summaries"
    if (!is.character(return_mode) || length(return_mode) != 1L || is.na(return_mode) || !return_mode %in%
        c("matrices_and_summaries", "complete")) {
        stop("return_mode must be either 'matrices_and_summaries' or 'complete'", call. = FALSE)
    }
    time_tolerance <- .smv21_resolve_time_tolerance(trajectory_obj, time_tolerance)
    intervals <- prepare_scenario_intervals_v2_1(trajectory_obj, maximum_transition_counts, time_tolerance)
    component_definitions <- .smv21_component_definitions(intervals$state_levels, maximum_transition_counts)
    tip_count <- length(intervals$tip_labels)
    new_tip_matrix <- function() {
        matrix(0, nrow = tip_count, ncol = tip_count, dimnames = list(intervals$tip_labels, intervals$tip_labels))
    }
    synchronous_matrices <- setNames(lapply(component_definitions$sync, function(name) new_tip_matrix()),
        component_definitions$sync)
    asynchronous_matrices <- list()
    for (view_name in names(component_definitions$async)) {
        asynchronous_matrices[[view_name]] <- setNames(lapply(component_definitions$async[[view_name]],
            function(name) new_tip_matrix()), component_definitions$async[[view_name]])
    }
    total_comparable_time_by_pair <- new_tip_matrix()
    asynchronous_available_similarity_by_pair <- new_tip_matrix()
    tip_history_times <- as.numeric(intervals$active %*% intervals$duration)
    names(tip_history_times) <- intervals$tip_labels
    tree_is_ultrametric <- max(tip_history_times) - min(tip_history_times) <= time_tolerance
    pair_detail_rows <- list(interval_classification = list(), sync_parallel_time_by_path = list(), async_independent_time_by_path = list(),
        async_parallel_similarity_by_path = list())
    pair_id <- 0L

    # Calculate each unordered pair once and copy its values symmetrically.
    #first_tip_index = 1
    #second_tip_index = 3
    for (first_tip_index in seq_len(tip_count - 1L)) {
        for (second_tip_index in seq.int(first_tip_index + 1L, tip_count)) {
            pair_id <- pair_id + 1L
            pair_result <- calculate_pair_similarity_v2_1(intervals, first_tip = first_tip_index,second_tip =  second_tip_index,
                async_metric)
            for (component_name in names(pair_result$sync)) {
                value <- pair_result$sync[[component_name]]
                synchronous_matrices[[component_name]][first_tip_index, second_tip_index] <- value
                synchronous_matrices[[component_name]][second_tip_index, first_tip_index] <- value
            }
            for (view_name in names(pair_result$async)) {
                for (component_name in names(pair_result$async[[view_name]])) {
                  value <- pair_result$async[[view_name]][[component_name]]
                  asynchronous_matrices[[view_name]][[component_name]][first_tip_index, second_tip_index] <- value
                  asynchronous_matrices[[view_name]][[component_name]][second_tip_index, first_tip_index] <- value
                }
            }
            total_comparable_time_by_pair[first_tip_index, second_tip_index] <- pair_result$total_comparable_time
            total_comparable_time_by_pair[second_tip_index, first_tip_index] <- pair_result$total_comparable_time
            available_asynchronous_similarity <- if (async_metric == "bhattacharyya") {
                sqrt(pair_result$first_total_history_time * pair_result$second_total_history_time)
            } else {
                min(pair_result$first_total_history_time, pair_result$second_total_history_time)
            }
            asynchronous_available_similarity_by_pair[first_tip_index, second_tip_index] <- available_asynchronous_similarity
            asynchronous_available_similarity_by_pair[second_tip_index, first_tip_index] <- available_asynchronous_similarity

            # Complete mode retains every readable pair ledger and path calculation.
            if (return_mode == "complete") {
                for (detail_name in names(pair_detail_rows)) {
                  detail <- pair_result$details[[detail_name]]
                  if (!nrow(detail))
                    next
                  detail$row_i <- first_tip_index
                  detail$row_j <- second_tip_index
                  detail$tip_i <- intervals$tip_labels[[first_tip_index]]
                  detail$tip_j <- intervals$tip_labels[[second_tip_index]]
                  pair_detail_rows[[detail_name]][[pair_id]] <- detail
                }
            }
        }
    }

    matrices <- list(sync = synchronous_matrices, async = asynchronous_matrices)
    synchronous_summary_weights <- if (tree_is_ultrametric) NULL else total_comparable_time_by_pair
    asynchronous_summary_weights <- if (tree_is_ultrametric) NULL else asynchronous_available_similarity_by_pair
    summaries <- list(sync = .smv21_summarize_family(synchronous_matrices, synchronous_summary_weights),
        async = lapply(asynchronous_matrices, .smv21_summarize_family,
            available_similarity_by_pair = asynchronous_summary_weights),
        components = component_definitions$metadata,
        weighting = list(tree_is_ultrametric = tree_is_ultrametric,
            available_similarity_weighting_applied = !tree_is_ultrametric,
            tip_history_times = tip_history_times,
            synchronous_available_similarity_by_pair = total_comparable_time_by_pair,
            asynchronous_available_similarity_by_pair = asynchronous_available_similarity_by_pair))

    if (return_mode == "matrices_and_summaries") {
        return(.smv21_add_depth_reports(list(matrices = matrices, summaries = summaries)))
    }

    bound_details <- list()
    for (detail_name in names(pair_detail_rows)) {
        rows <- pair_detail_rows[[detail_name]]
        if (length(rows)) {
            bound_details[[detail_name]] <- do.call(rbind, rows)
        } else {
            bound_details[[detail_name]] <- data.frame()
        }
    }
    result <- list(intervals = intervals, matrices = matrices, total_comparable_time_by_pair = total_comparable_time_by_pair,
        available_similarity_by_pair = list(synchronous = total_comparable_time_by_pair,
            asynchronous = asynchronous_available_similarity_by_pair),
        pair_details = bound_details, summaries = summaries, components = component_definitions$metadata,
        metadata = list(maximum_transition_counts = maximum_transition_counts, async_metric = async_metric,
            tree_is_ultrametric = tree_is_ultrametric,
            available_similarity_weighting_applied = !tree_is_ultrametric,
            normalization = list(synchronous = "shared_active_time",
                asynchronous_bhattacharyya = "sqrt(first_history_time * second_history_time)",
                asynchronous_minimum = "min(first_history_time, second_history_time)"),
            time_tolerance = time_tolerance))
    .smv21_add_depth_reports(result)
}

