#' Extract summary-statistic surfaces
#'
#' @param x A `trajectory_tree`, `trajectory_distribution`, or a list of either.
#' @param component Summary component to return. `"all"` returns the complete
#'   summary surface.
#'
#' @return A summary list, one requested component, or a shape-preserving list.
#' @export
trajectory_summary <- function(
    x,
    component = c(
      "all", "original_metrics", "tree_wide", "by_state",
      "by_state_norm", "by_path", "by_event"
    )) {
  component <- match.arg(component)
  if (inherits(x, "trajectory_tree")) {
    surface <- x$ss
  } else if (inherits(x, "trajectory_distribution")) {
    surface <- x$SS
  } else if (is.list(x)) {
    return(lapply(x, trajectory_summary, component = component))
  } else {
    stop("x must be a trajectory tree, distribution, or list", call. = FALSE)
  }
  if (is.null(surface)) {
    stop("x has no summary-statistic surface", call. = FALSE)
  }
  if (identical(component, "all")) {
    return(surface)
  }
  if (!component %in% names(surface)) {
    stop("summary component is unavailable: ", component, call. = FALSE)
  }
  surface[[component]]
}

#' Extract through-time surfaces
#'
#' @param x A `trajectory_tree`, `trajectory_distribution`, or list of either.
#' @param metric Optional metric such as `"ltt"`, `"stt"`, `"tt"`, `"ctt"`,
#'   `"ldif"`, `"puniq"`, or `"tlen"`.
#' @param component Optional component such as `"tot"`, `"state"`, `"event"`,
#'   `"arr"`, `"leave"`, or `"trans"`.
#'
#' @return The requested through-time surface.
#' @export
trajectory_through_time <- function(x, metric = NULL, component = NULL) {
  if (inherits(x, "trajectory_tree")) {
    surface <- x$tt
  } else if (inherits(x, "trajectory_distribution")) {
    surface <- x$TT
  } else if (is.list(x)) {
    return(lapply(
      x, trajectory_through_time,
      metric = metric, component = component
    ))
  } else {
    stop("x must be a trajectory tree, distribution, or list", call. = FALSE)
  }
  if (is.null(surface)) {
    stop("x has no through-time surface", call. = FALSE)
  }
  if (is.null(metric)) {
    if (!is.null(component)) {
      stop("component requires a metric", call. = FALSE)
    }
    return(surface)
  }
  if (!metric %in% names(surface)) {
    stop("through-time metric is unavailable: ", metric, call. = FALSE)
  }
  result <- surface[[metric]]
  if (is.null(component)) {
    return(result)
  }
  if (!is.list(result) || !component %in% names(result)) {
    stop("through-time component is unavailable: ", component, call. = FALSE)
  }
  result[[component]]
}

#' Build posterior trajectory distributions
#'
#' A flat list of trajectory trees produces one distribution. A named list of
#' trajectory-tree lists produces a named list of distributions.
#'
#' @param trajectory_objects Nonempty trajectory-tree list or list of such
#'   lists.
#' @param probs Lower and upper empirical quantile probabilities.
#' @param tolerance Nonnegative alignment tolerance.
#' @param include Any of `"P"`, `"T"`, `"SS"`, and `"TT"`.
#' @param backend Exact reduction backend, `"cpp"` or `"R"`.
#' @param block_bytes Positive C++ working-buffer target in bytes.
#'
#' @return A `trajectory_distribution` or a shape-preserving named list.
#' @export
make_trajectory_distributions <- function(
    trajectory_objects,
    probs = c(0.025, 0.975),
    tolerance = 1e-8,
    include = c("P", "T", "SS", "TT"),
    backend = "cpp",
    block_bytes = 16 * 1024^2) {
  build <- function(samples) {
    if (!is.list(samples) || !length(samples) ||
        !all(vapply(samples, inherits, logical(1), "trajectory_tree"))) {
      stop("each distribution requires a nonempty trajectory_tree list", call. = FALSE)
    }
    result <- get_pst_distribution(
      samples = samples,
      probs = probs,
      tolerance = tolerance,
      include = include,
      backend = backend,
      block_bytes = block_bytes
    )
    class(result) <- unique(c(
      "trajectory_distribution",
      setdiff(class(result), "pst_distribution")
    ))
    result
  }

  if (is.list(trajectory_objects) && length(trajectory_objects) &&
      all(vapply(trajectory_objects, inherits, logical(1), "trajectory_tree"))) {
    return(build(trajectory_objects))
  }
  if (is.list(trajectory_objects) && length(trajectory_objects) &&
      all(vapply(trajectory_objects, function(value) {
        is.list(value) && length(value) &&
          all(vapply(value, inherits, logical(1), "trajectory_tree"))
      }, logical(1)))) {
    return(lapply(trajectory_objects, build))
  }
  stop(
    "trajectory_objects must be a trajectory_tree list or a list of trajectory_tree lists",
    call. = FALSE
  )
}

#' Compare two trajectory distributions
#'
#' Matching numeric leaves are compared recursively, producing one compact
#' diagnostic row per surface. Structural leaves that occur in only one object
#' are also reported.
#'
#' @param reference,candidate `trajectory_distribution` objects.
#' @param surfaces Any of `"P"`, `"T"`, `"SS"`, and `"TT"`.
#' @param tolerance Absolute numeric comparison tolerance.
#'
#' @return A data frame with paths, dimensions, maximum absolute differences,
#'   and equality flags.
#' @export
compare_trajectory_distributions <- function(
    reference,
    candidate,
    surfaces = c("P", "T", "SS", "TT"),
    tolerance = 1e-8) {
  if (!inherits(reference, "trajectory_distribution") ||
      !inherits(candidate, "trajectory_distribution")) {
    stop("reference and candidate must be trajectory_distribution objects", call. = FALSE)
  }
  surfaces <- unique(toupper(surfaces))
  if (!length(surfaces) || any(!surfaces %in% c("P", "T", "SS", "TT"))) {
    stop("surfaces may contain only P, T, SS, and TT", call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be one nonnegative finite number", call. = FALSE)
  }
  rows <- unlist(lapply(surfaces, function(surface) {
    trajectory_compare_nodes(
      reference[[surface]], candidate[[surface]], surface, tolerance
    )
  }), recursive = FALSE)
  if (!length(rows)) {
    return(data.frame(
      path = character(), reference_length = integer(),
      candidate_length = integer(), max_abs_diff = numeric(),
      equal = logical()
    ))
  }
  do.call(rbind, rows)
}

trajectory_compare_nodes <- function(reference, candidate, path, tolerance) {
  if (is.numeric(reference) || is.numeric(candidate)) {
    same_shape <- identical(dim(reference), dim(candidate)) &&
      length(reference) == length(candidate)
    difference <- if (same_shape) {
      suppressWarnings(max(abs(reference - candidate), na.rm = TRUE))
    } else {
      Inf
    }
    if (!is.finite(difference) && same_shape &&
        all(is.na(reference) == is.na(candidate))) {
      difference <- 0
    }
    return(list(data.frame(
      path = path,
      reference_length = length(reference),
      candidate_length = length(candidate),
      max_abs_diff = difference,
      equal = same_shape &&
        identical(is.na(reference), is.na(candidate)) &&
        difference <= tolerance,
      stringsAsFactors = FALSE
    )))
  }
  if (is.list(reference) || is.list(candidate)) {
    fields <- union(names(reference), names(candidate))
    if (!length(fields)) {
      return(list())
    }
    return(unlist(lapply(fields, function(field) {
      trajectory_compare_nodes(
        reference[[field]], candidate[[field]],
        paste(path, field, sep = "/"), tolerance
      )
    }), recursive = FALSE))
  }
  list()
}

#' @export
print.trajectory_distribution <- function(x, ...) {
  cat("<trajectory_distribution>\n")
  cat("  samples:", x$metadata$n_samples, "\n")
  cat(
    "  surfaces:",
    paste(names(x)[vapply(x, Negate(is.null), logical(1))], collapse = ", "),
    "\n"
  )
  invisible(x)
}
