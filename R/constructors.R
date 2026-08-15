#' Construct trajectory trees from stochastic character maps
#'
#' `make_trajectory_objects()` is the main construction entry point. It accepts
#' one `simmap`, a list of `simmap` objects, or a list of such lists. Names and
#' nesting are retained. The production implementation is private so public
#' code does not depend on an implementation-version suffix.
#'
#' @param simmaps A `simmap`, a nonempty list of `simmap` objects, or a
#'   nonempty list of nonempty `simmap` lists.
#' @param root_state Optional explicit root state.
#' @param include_tt,include_ss Whether to construct through-time and summary
#'   statistic surfaces.
#' @param include_scenario_mats Whether to retain dense scenario-coordinate
#'   matrices.
#' @param include_path_maps Whether to retain path-labelled maps.
#' @param time_tolerance Positive event-time resolution. `NULL` chooses the
#'   smaller of `1e-10` and one tenth of the shortest positive map duration.
#' @param keep_internal_attributes Whether to retain private V35 construction,
#'   timing, provenance, and diagnostic attributes. Defaults to `FALSE`.
#' @param fields Optional top-level fields to retain after construction.
#' @param ... Advanced construction controls.
#'
#' @return One `trajectory_tree`, a named list of them, or a shape-preserving
#'   named list of lists.
#' @export
make_trajectory_objects <- function(
    simmaps,
    root_state = NULL,
    include_tt = TRUE,
    include_ss = TRUE,
    include_scenario_mats = TRUE,
    include_path_maps = TRUE,
    time_tolerance = NULL,
    keep_internal_attributes = FALSE,
    fields = NULL,
    ...) {
  normalized <- trajectory_normalize_simmaps(simmaps)
  built <- lapply(normalized$sets, function(simmap_set) {
    lapply(simmap_set, function(simmap) {
      tolerance <- trajectory_resolve_tolerance(simmap, time_tolerance)
      result <- get_trajectory_obj_v35(
        phylo_charmap = simmap,
        root_state = root_state,
        include_tt = include_tt,
        include_ss = include_ss,
        include_scenario_mats = include_scenario_mats,
        include_path_maps = include_path_maps,
        time_tolerance = tolerance,
        run_tests = FALSE,
        keep_internal_attributes = keep_internal_attributes,
        ...
      )
      class(result) <- unique(c("trajectory_tree", class(result)))
      if (!is.null(fields)) {
        result <- select_trajectory_fields(result, include = fields)
      }
      result
    })
  })
  trajectory_restore_shape(built, normalized)
}

#' Select fields from a trajectory tree
#'
#' @param trajectory A `trajectory_tree`.
#' @param include Character vector of top-level fields to retain.
#' @param ss_include Optional names within `ss` to retain.
#' @param tt_include Optional names within `tt` to retain. When TT is retained,
#'   its shared `time` coordinate is always kept.
#' @param keep_internal_attributes Whether to retain private construction
#'   attributes. The default keeps only names and class.
#'
#' @return A compact `trajectory_tree`.
#' @export
select_trajectory_fields <- function(
    trajectory,
    include = names(trajectory),
    ss_include = NULL,
    tt_include = NULL,
    keep_internal_attributes = FALSE) {
  if (!inherits(trajectory, "trajectory_tree") || !is.list(trajectory)) {
    stop("trajectory must be a trajectory_tree", call. = FALSE)
  }
  if (!is.character(include) || anyNA(include)) {
    stop("include must be a character vector without missing values", call. = FALSE)
  }
  unknown <- setdiff(include, names(trajectory))
  if (length(unknown)) {
    stop("unknown trajectory field(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  original_attributes <- attributes(trajectory)
  result <- trajectory[include]
  if ("ss" %in% names(result) && !is.null(ss_include)) {
    result$ss <- trajectory_select_subfields(result$ss, ss_include, "ss_include")
  }
  if ("tt" %in% names(result) && !is.null(tt_include)) {
    tt_include <- unique(c("time", tt_include))
    result$tt <- trajectory_select_subfields(result$tt, tt_include, "tt_include")
  }
  if (isTRUE(keep_internal_attributes)) {
    original_attributes$names <- names(result)
    attributes(result) <- original_attributes
  } else {
    class(result) <- unique(c("trajectory_tree", class(result)))
  }
  result
}

trajectory_select_subfields <- function(value, include, argument) {
  if (!is.list(value)) {
    stop(argument, " cannot be used because the surface is not a list", call. = FALSE)
  }
  if (!is.character(include) || anyNA(include)) {
    stop(argument, " must be a character vector without missing values", call. = FALSE)
  }
  unknown <- setdiff(include, names(value))
  if (length(unknown)) {
    stop("unknown ", argument, " field(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  value[include]
}

trajectory_is_simmap <- function(value) {
  inherits(value, "simmap") && inherits(value, "phylo") &&
    is.matrix(value$edge) && is.list(value$maps)
}

trajectory_is_simmap_set <- function(value) {
  is.list(value) && length(value) > 0L &&
    all(vapply(value, trajectory_is_simmap, logical(1)))
}

trajectory_normalize_simmaps <- function(simmaps) {
  if (trajectory_is_simmap(simmaps)) {
    sets <- list(list(simmaps))
    shape <- "single"
  } else if (trajectory_is_simmap_set(simmaps)) {
    sets <- list(simmaps)
    shape <- "flat"
  } else if (is.list(simmaps) && length(simmaps) > 0L &&
             all(vapply(simmaps, trajectory_is_simmap_set, logical(1)))) {
    sets <- unname(simmaps)
    shape <- "nested"
  } else {
    stop(
      "simmaps must be a simmap, a nonempty simmap list, or a uniformly nested nonempty list of simmap lists",
      call. = FALSE
    )
  }
  list(
    shape = shape,
    sets = sets,
    outer_names = if (identical(shape, "nested")) names(simmaps) else NULL,
    inner_names = lapply(sets, names)
  )
}

trajectory_restore_shape <- function(sets, normalized) {
  for (index in seq_along(sets)) {
    names(sets[[index]]) <- normalized$inner_names[[index]]
  }
  if (identical(normalized$shape, "single")) {
    return(sets[[1L]][[1L]])
  }
  if (identical(normalized$shape, "flat")) {
    return(sets[[1L]])
  }
  names(sets) <- normalized$outer_names
  sets
}

trajectory_resolve_tolerance <- function(simmap, tolerance) {
  if (!is.null(tolerance)) {
    if (!is.numeric(tolerance) || length(tolerance) != 1L ||
        !is.finite(tolerance) || tolerance <= 0) {
      stop("time_tolerance must be NULL or one finite positive number", call. = FALSE)
    }
    return(as.numeric(tolerance))
  }
  durations <- unlist(simmap$maps, use.names = FALSE)
  positive <- durations[is.finite(durations) & durations > 0]
  if (!length(positive)) {
    return(1e-10)
  }
  min(1e-10, min(positive) / 10)
}

#' @export
print.trajectory_tree <- function(x, ...) {
  tip_count <- if (inherits(x$phylo, "phylo")) length(x$phylo$tip.label) else NA_integer_
  states <- if (!is.null(x$phylo$mapped.edge)) colnames(x$phylo$mapped.edge) else character()
  cat("<trajectory_tree>\n")
  cat("  tips:", tip_count, "\n")
  cat("  states:", if (length(states)) paste(states, collapse = ", ") else "not available", "\n")
  cat("  surfaces:", paste(names(x), collapse = ", "), "\n")
  invisible(x)
}
