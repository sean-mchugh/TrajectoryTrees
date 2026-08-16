# Package-native ScenarioMats similarity interface.

.scenario_mats_v3_add_similarity_vectors <- function(
    result, state_levels = NULL, tolerance = 1e-8) {
  metadata <- result$summaries$components
  required <- c(
    "matrix_name", "state", "phase", "kind",
    "requested_maximum_transition_count",
    "actual_transition_count_used", "view"
  )
  if (!all(required %in% names(metadata))) {
    stop("component metadata cannot build similarity vectors", call. = FALSE)
  }

  component_group <- rep(NA_character_, nrow(metadata))
  component_group[metadata$phase == "Shared"] <- "shared"
  component_group[metadata$phase == "Conserved"] <- "conserved"
  component_group[
    metadata$phase == "Conserved_Homoplasy"
  ] <- "conserved_homoplasy"
  component_group[metadata$kind == "parallel"] <- "parallel"
  component_group[metadata$kind == "convergent"] <- "convergent"
  component_group[
    metadata$phase == "Independent" &
      metadata$kind == "state" &
      metadata$actual_transition_count_used == 0
  ] <- "convergent"
  component_group[metadata$kind == "residual"] <- "divergent"
  if (anyNA(component_group)) {
    stop("a similarity-vector component is unclassified", call. = FALSE)
  }
  metadata$component_group <- component_group

  observed_states <- unique(metadata$state[!is.na(metadata$state)])
  states <- if (is.null(state_levels)) {
    observed_states
  } else {
    unique(as.character(state_levels))
  }
  if (!length(states) || anyNA(states) || any(!nzchar(states)) ||
      !all(observed_states %in% states)) {
    stop(
      "state_levels must contain every observed nonempty state label",
      call. = FALSE
    )
  }
  positive_depths <- sort(unique(
    metadata$requested_maximum_transition_count[
      !is.na(metadata$requested_maximum_transition_count) &
        metadata$requested_maximum_transition_count > 0
    ]
  ))
  maximum_depth <- if (length(positive_depths)) {
    max(positive_depths)
  } else {
    0L
  }
  if (maximum_depth > 0L &&
      !identical(as.integer(positive_depths), seq_len(maximum_depth))) {
    result$summaries$similarity_vectors <- list(
      available = FALSE,
      reason = paste(
        "combined vectors require contiguous positive depths from 1",
        "through the maximum requested depth"
      ),
      requested_positive_depths = as.integer(positive_depths)
    )
    return(result)
  }

  view_summary <- function(view) {
    if (view == "sync") {
      result$summaries$sync
    } else {
      result$summaries$async[[view]]
    }
  }
  component_value <- function(
      view, group, state = NULL,
      requested_depth = NULL, actual_depth = NULL) {
    rows <- metadata$view == view &
      metadata$component_group == group
    if (!is.null(state)) {
      rows <- rows & !is.na(metadata$state) & metadata$state == state
    }
    if (!is.null(requested_depth)) {
      rows <- rows &
        !is.na(metadata$requested_maximum_transition_count) &
        metadata$requested_maximum_transition_count == requested_depth
    }
    if (!is.null(actual_depth)) {
      rows <- rows & !is.na(metadata$actual_transition_count_used) &
        metadata$actual_transition_count_used == actual_depth
    }
    summary <- view_summary(view)
    proportions <- summary$tree_wide_proportions
    if (is.null(proportions)) proportions <- summary$means
    sum(proportions[metadata$matrix_name[rows]])
  }
  depth_view <- function(depth) {
    if (depth == 1L) {
      "using_up_to_1_transition"
    } else {
      paste0("using_up_to_", depth, "_transitions")
    }
  }
  parallel_names <- if (maximum_depth) {
    paste0("parallel_", seq_len(maximum_depth))
  } else {
    character()
  }
  async_parallel_names <- parallel_names
  if (maximum_depth) {
    async_parallel_names[[maximum_depth]] <- paste0(
      async_parallel_names[[maximum_depth]], "_plus"
    )
  }

  synchronous_state <- do.call(rbind, lapply(states, function(state) {
    parallel <- if (maximum_depth) {
      vapply(seq_len(maximum_depth), function(depth) {
        component_value(
          "sync", "parallel", state = state, actual_depth = depth
        )
      }, numeric(1))
    } else {
      numeric()
    }
    names(parallel) <- parallel_names
    c(
      shared = component_value("sync", "shared", state = state),
      conserved = component_value("sync", "conserved", state = state),
      conserved_homoplasy = component_value(
        "sync", "conserved_homoplasy", state = state
      ),
      convergent = component_value(
        "sync", "convergent", state = state
      ),
      parallel
    )
  }))
  rownames(synchronous_state) <- states

  asynchronous_state <- do.call(rbind, lapply(states, function(state) {
    state_only_convergence <- component_value(
      "state_only", "convergent", state = state
    )
    nested <- if (maximum_depth) {
      vapply(seq_len(maximum_depth), function(depth) {
        component_value(
          depth_view(depth), "parallel", state = state,
          requested_depth = depth, actual_depth = depth
        )
      }, numeric(1))
    } else {
      numeric()
    }
    if (maximum_depth &&
        (nested[[1L]] > state_only_convergence + tolerance ||
          (maximum_depth > 1L && any(diff(nested) > tolerance)))) {
      stop(
        "asynchronous terminal parallel bins are not nested for state ",
        state,
        call. = FALSE
      )
    }
    if (maximum_depth) {
      nested[[1L]] <- min(state_only_convergence, nested[[1L]])
      if (maximum_depth > 1L) {
        for (depth in 2:maximum_depth) {
          nested[[depth]] <- min(nested[[depth - 1L]], nested[[depth]])
        }
      }
      parallel <- numeric(maximum_depth)
      if (maximum_depth > 1L) {
        parallel[seq_len(maximum_depth - 1L)] <-
          nested[seq_len(maximum_depth - 1L)] -
          nested[seq.int(2L, maximum_depth)]
      }
      parallel[[maximum_depth]] <- nested[[maximum_depth]]
      parallel[abs(parallel) <= tolerance] <- 0
      names(parallel) <- async_parallel_names
      convergent <- state_only_convergence - nested[[1L]]
    } else {
      parallel <- numeric()
      convergent <- state_only_convergence
    }
    if (convergent < -tolerance || any(parallel < -tolerance)) {
      stop("combined asynchronous vector is negative", call. = FALSE)
    }
    c(
      shared = component_value("state_only", "shared", state = state),
      conserved = component_value(
        "state_only", "conserved", state = state
      ),
      conserved_homoplasy = component_value(
        "state_only", "conserved_homoplasy", state = state
      ),
      convergent = max(0, convergent),
      parallel
    )
  }))
  rownames(asynchronous_state) <- states

  cumulative_state <- function(values) {
    parallel_columns <- grepl("^parallel_", colnames(values))
    parallel <- if (any(parallel_columns)) {
      rowSums(values[, parallel_columns, drop = FALSE])
    } else {
      rep(0, nrow(values))
    }
    cbind(
      all_parallel = parallel,
      convergent_plus_parallel = values[, "convergent"] + parallel,
      conserved_homoplasy_plus_convergent_plus_parallel =
        values[, "conserved_homoplasy"] +
        values[, "convergent"] + parallel
    )
  }
  all_states <- function(state_values, divergent, view) {
    values <- c(colSums(state_values), divergent = divergent)
    total <- sum(values)
    if (abs(total - 1) > tolerance) {
      condition <- structure(
        list(
          message = paste0(
            "the ", view, " all-states similarity vector closes to ",
            format(total, digits = 17), " rather than one"
          ),
          call = NULL,
          view = view,
          values = values,
          total = total,
          tolerance = tolerance
        ),
        class = c(
          "trajectorytrees_similarity_closure_error",
          "error",
          "condition"
        )
      )
      stop(condition)
    }
    values
  }
  cumulative_all_states <- function(values) {
    parallel <- sum(values[grepl("^parallel_", names(values))])
    c(
      all_parallel = parallel,
      convergent_plus_parallel = values[["convergent"]] + parallel,
      conserved_homoplasy_plus_convergent_plus_parallel =
        values[["conserved_homoplasy"]] +
        values[["convergent"]] + parallel
    )
  }

  synchronous_all <- all_states(
    synchronous_state,
    component_value("sync", "divergent"),
    "sync"
  )
  asynchronous_all <- all_states(
    asynchronous_state,
    component_value("state_only", "divergent"),
    "async_combined"
  )
  result$summaries$similarity_vectors <- list(
    available = TRUE,
    maximum_depth = maximum_depth,
    normalization = "global_tree_wide",
    sync = list(
      state = list(
        mutually_exclusive = synchronous_state,
        cumulative = cumulative_state(synchronous_state)
      ),
      all_states = list(
        mutually_exclusive = synchronous_all,
        cumulative = cumulative_all_states(synchronous_all)
      )
    ),
    async_combined = list(
      state = list(
        mutually_exclusive = asynchronous_state,
        cumulative = cumulative_state(asynchronous_state)
      ),
      all_states = list(
        mutually_exclusive = asynchronous_all,
        cumulative = cumulative_all_states(asynchronous_all)
      )
    )
  )
  result
}

.scenario_mats_v3_add_depth_reports <- function(result) {
  if (is.null(result$matrices)) return(result)
  metadata <- result$summaries$components
  asynchronous_weights <- if (
      isTRUE(result$summaries$weighting$
        available_similarity_weighting_applied)) {
    result$summaries$weighting$asynchronous_available_similarity_by_pair
  } else {
    NULL
  }
  matrix_reports <- list()
  summary_reports <- list()
  summarize_family <- function(matrices, weights = NULL) {
    upper_triangle <- upper.tri(matrices[[1L]])
    totals <- vapply(
      matrices,
      function(matrix) {
        if (is.null(weights)) {
          sum(matrix[upper_triangle])
        } else {
          sum(matrix[upper_triangle] * weights[upper_triangle])
        }
      },
      numeric(1)
    )
    total <- sum(totals)
    proportions <- if (total > 0) totals / total else totals
    summary <- list(
      totals = totals,
      tree_wide_proportions = proportions,
      means = proportions,
      total_across_matrices = total
    )
    if (!is.null(weights)) {
      summary$available_similarity_total <- sum(weights[upper_triangle])
    }
    summary
  }

  for (view_name in names(result$matrices$async)) {
    rows <- metadata[
      metadata$view == view_name &
        metadata$kind %in% c("convergent", "parallel"),
      ,
      drop = FALSE
    ]
    requested_depth <- unique(rows$requested_maximum_transition_count)
    requested_depth <- requested_depth[!is.na(requested_depth)]
    if (length(requested_depth) != 1L || requested_depth == 0L) next

    state_matrices <- result$matrices$async[[view_name]][rows$matrix_name]
    total_matrices <- list()
    for (row_index in seq_len(nrow(rows))) {
      total_name <- if (rows$kind[[row_index]] == "convergent") {
        "Convergent"
      } else if (
        rows$actual_transition_count_used[[row_index]] == requested_depth
      ) {
        paste0("Parallel_", requested_depth, "OrMoreTransitions")
      } else {
        transition_count <- rows$actual_transition_count_used[[row_index]]
        paste0(
          "Parallel_", transition_count,
          if (transition_count == 1L) "Transition" else "Transitions"
        )
      }
      if (is.null(total_matrices[[total_name]])) {
        total_matrices[[total_name]] <- state_matrices[[row_index]]
      } else {
        total_matrices[[total_name]] <-
          total_matrices[[total_name]] + state_matrices[[row_index]]
      }
    }
    matrix_reports[[view_name]] <- list(
      state = state_matrices,
      total = total_matrices
    )
    summary_reports[[view_name]] <- list(
      state = summarize_family(state_matrices, asynchronous_weights),
      total = summarize_family(total_matrices, asynchronous_weights)
    )
  }
  result$matrices$async_convergence_by_transition_depth <- matrix_reports
  result$summaries$async_convergence_by_transition_depth <- summary_reports
  result
}

#' Calculate similarity among evolutionary trajectories
#'
#' @param trajectory_obj A trajectory object containing `phylo`,
#'   `scenario_mats`, and `root_policy`. No PST version is required.
#' @param maximum_transition_counts Unique nonnegative whole-number path depths,
#'   including zero for the state-only view.
#' @param async_metric Either `"bhattacharyya"` or `"minimum"`.
#' @param return_mode Return matrices and summaries by default, summaries only
#'   with `"summaries_only"`, or all calculation details with `"complete"`.
#' @param time_tolerance Positive numeric resolution used for time comparisons.
#'   When omitted, a version-neutral tolerance stored on the trajectory is used;
#'   legacy PST attributes remain readable for saved-object compatibility.
#' @return A ScenarioMats similarity result.
#' @export
trajectory_similarity <- function(
    trajectory_obj,
    maximum_transition_counts = 0L,
    async_metric = c("bhattacharyya", "minimum"),
    return_mode = c("matrices_and_summaries", "summaries_only", "complete"),
    time_tolerance = NULL) {
  if (missing(trajectory_obj)) {
    stop("argument \"trajectory_obj\" is missing", call. = FALSE)
  }
  if (!is.list(trajectory_obj) ||
      is.null(trajectory_obj$phylo) ||
      is.null(trajectory_obj$scenario_mats) ||
      is.null(trajectory_obj$root_policy)) {
    stop(
      "trajectory_obj must contain phylo, scenario_mats, and root_policy",
      call. = FALSE
    )
  }
  if (!is.numeric(maximum_transition_counts) ||
      !length(maximum_transition_counts) ||
      any(!is.finite(maximum_transition_counts)) ||
      any(maximum_transition_counts < 0) ||
      any(maximum_transition_counts != floor(maximum_transition_counts))) {
    stop(
      "maximum_transition_counts must contain non-negative whole numbers",
      call. = FALSE
    )
  }
  if (any(maximum_transition_counts >= .Machine$integer.max)) {
    stop(
      paste(
        "maximum_transition_counts must be representable R integers smaller",
        "than .Machine$integer.max"
      ),
      call. = FALSE
    )
  }
  if (anyDuplicated(maximum_transition_counts) ||
      !0L %in% maximum_transition_counts) {
    stop(
      paste(
        "maximum_transition_counts must be unique and include 0 for the",
        "state-only comparison"
      ),
      call. = FALSE
    )
  }

  if (missing(async_metric)) {
    async_metric <- "bhattacharyya"
  } else if (!is.character(async_metric) || length(async_metric) != 1L ||
             is.na(async_metric) ||
             !async_metric %in% c("bhattacharyya", "minimum")) {
    stop(
      "async_metric must be either 'bhattacharyya' or 'minimum'",
      call. = FALSE
    )
  }
  if (missing(return_mode)) {
    return_mode <- "matrices_and_summaries"
  } else if (!is.character(return_mode) || length(return_mode) != 1L ||
             is.na(return_mode) ||
             !return_mode %in% c("matrices_and_summaries", "summaries_only", "complete")) {
    stop(
      paste(
        "return_mode must be either 'matrices_and_summaries',",
        "'summaries_only', or 'complete'"
      ),
      call. = FALSE
    )
  }
  maximum_transition_counts <- as.integer(maximum_transition_counts)

  if (is.null(time_tolerance)) {
    time_tolerance <- trajectory_obj$time_tolerance
  }
  if (is.null(time_tolerance)) {
    time_tolerance <- trajectory_obj$scenario_mats$time_tolerance
  }
  if (is.null(time_tolerance)) {
    time_tolerance <- attr(
      trajectory_obj, "scenario_mats_time_tolerance", exact = TRUE
    )
  }
  # Saved trajectories created before the version-neutral contract used a
  # version-labelled attribute. Keep those objects readable without exposing
  # their version labels to C++ or requiring them from new callers.
  if (is.null(time_tolerance)) {
    legacy_names <- paste0("pst_v", c(35L, 34L, 33L, 32L), "_time_tolerance")
    for (legacy_name in legacy_names) {
      time_tolerance <- attr(trajectory_obj, legacy_name, exact = TRUE)
      if (!is.null(time_tolerance)) break
    }
  }
  if (!is.numeric(time_tolerance) || length(time_tolerance) != 1L ||
      !is.finite(time_tolerance) || time_tolerance <= 0) {
    stop(
      paste(
        "time_tolerance must be one finite positive number, supplied directly",
        "or stored on the trajectory object"
      ),
      call. = FALSE
    )
  }

  scenario_mats_v3_load_cpp()
  result <- scenario_mats_v3_calculate_cpp(
    trajectory_obj,
    maximum_transition_counts,
    async_metric,
    identical(return_mode, "complete"),
    !identical(return_mode, "summaries_only"),
    as.numeric(time_tolerance)
  )
  result <- .scenario_mats_v3_add_depth_reports(result)
  .scenario_mats_v3_add_similarity_vectors(
    result,
    state_levels = colnames(trajectory_obj$phylo$mapped.edge)
  )
}
