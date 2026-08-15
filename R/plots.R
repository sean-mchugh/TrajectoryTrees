#' Plot one tree from a trajectory object
#'
#' @param x A `trajectory_tree` or a `phylo`/`simmap`.
#' @param tree Which trajectory tree to draw.
#' @param colors Optional named state-color vector.
#' @param ... Plot arguments passed to `phytools::plotSimmap()` when phytools is
#'   installed, otherwise to [ape::plot.phylo()].
#'
#' @return The plotted tree, invisibly.
#' @export
plot_trajectory_tree <- function(
    x,
    tree = c("phylo", "scenario", "trans"),
    colors = NULL,
    ...) {
  tree <- match.arg(tree)
  value <- if (inherits(x, "trajectory_tree")) x[[tree]] else x
  if (!inherits(value, "phylo")) {
    stop("the selected value is not a phylogenetic tree", call. = FALSE)
  }
  if (inherits(value, "simmap") && requireNamespace("phytools", quietly = TRUE)) {
    phytools::plotSimmap(value, colors = colors, ...)
  } else {
    ape::plot.phylo(value, ...)
  }
  invisible(value)
}

#' Plot a through-time trajectory surface
#'
#' @param x A `trajectory_tree`.
#' @param metric Through-time metric.
#' @param component Metric component.
#' @param series Optional series names to retain.
#' @param title Optional plot title.
#'
#' @return A `ggplot` object.
#' @export
plot_trajectory_through_time <- function(
    x,
    metric = "ltt",
    component = "tot",
    series = NULL,
    title = NULL) {
  if (!inherits(x, "trajectory_tree")) {
    stop("x must be a trajectory_tree", call. = FALSE)
  }
  values <- trajectory_through_time(x, metric, component)
  data <- trajectory_matrix_long(values, x$tt$time)
  if (!is.null(series)) {
    data <- data[data$series %in% series, , drop = FALSE]
  }
  ggplot2::ggplot(
    data,
    ggplot2::aes(x = time, y = value, color = series, group = series)
  ) +
    ggplot2::geom_step(direction = "hv", linewidth = 0.7) +
    ggplot2::labs(
      x = "Time from root", y = metric, color = component, title = title
    ) +
    ggplot2::theme_bw()
}

#' Plot an interval surface from a trajectory distribution
#'
#' @param x A `trajectory_distribution`.
#' @param surface One of `"P"`, `"T"`, `"SS"`, or `"TT"`.
#' @param path Optional character vector naming nested fields beneath the
#'   selected surface. With `NULL`, the first plottable interval node is used.
#' @param title Optional plot title.
#'
#' @return A `ggplot` object.
#' @export
plot_trajectory_distribution <- function(
    x,
    surface = c("SS", "TT", "T", "P"),
    path = NULL,
    title = NULL) {
  if (!inherits(x, "trajectory_distribution")) {
    stop("x must be a trajectory_distribution", call. = FALSE)
  }
  surface <- match.arg(surface)
  node <- x[[surface]]
  if (is.null(node)) {
    stop("distribution surface is unavailable: ", surface, call. = FALSE)
  }
  if (!is.null(path)) {
    for (field in path) {
      if (!is.list(node) || !field %in% names(node)) {
        stop("distribution path is unavailable: ", paste(path, collapse = "/"), call. = FALSE)
      }
      node <- node[[field]]
    }
  } else {
    found <- trajectory_find_interval_node(node)
    if (is.null(found)) {
      stop("surface has no plottable lower/mean/upper interval node", call. = FALSE)
    }
    node <- found$value
    path <- found$path
  }
  data <- trajectory_interval_long(node)
  ggplot2::ggplot(
    data,
    ggplot2::aes(x = index, y = mean, color = series, group = series)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper, fill = series),
      alpha = 0.18, color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::labs(
      x = "Index", y = "Posterior value", color = "Series", fill = "Series",
      title = title %||% paste(c(surface, path), collapse = " / ")
    ) +
    ggplot2::theme_bw()
}

#' Plot two trajectory trees as a tanglegram
#'
#' @param x A `trajectory_tree`.
#' @param left,right Trees to place on the left and right.
#' @param association Optional two-column character matrix of associated tip
#'   labels. When omitted, exact and fixed-substring label matches are used.
#' @param ... Arguments passed to `phytools::plot.cophylo()`.
#'
#' @return A `cophylo` object, invisibly.
#' @export
plot_trajectory_tanglegram <- function(
    x,
    left = c("phylo", "scenario", "trans"),
    right = c("scenario", "trans", "phylo"),
    association = NULL,
    ...) {
  if (!inherits(x, "trajectory_tree")) {
    stop("x must be a trajectory_tree", call. = FALSE)
  }
  left <- match.arg(left)
  right <- match.arg(right)
  if (identical(left, right)) {
    stop("left and right must select different trees", call. = FALSE)
  }
  if (!requireNamespace("phytools", quietly = TRUE)) {
    stop("plot_trajectory_tanglegram() requires the suggested phytools package", call. = FALSE)
  }
  left_tree <- x[[left]]
  right_tree <- x[[right]]
  if (!inherits(left_tree, "phylo") || !inherits(right_tree, "phylo")) {
    stop("both selected trajectory fields must be phylogenetic trees", call. = FALSE)
  }
  if (is.null(association)) {
    association <- trajectory_tip_association(left_tree, right_tree)
  }
  if (!is.matrix(association) || ncol(association) != 2L || !nrow(association)) {
    stop("association must be a nonempty two-column matrix", call. = FALSE)
  }
  object <- phytools::cophylo(
    left_tree, right_tree, assoc = association, rotate = TRUE
  )
  plot(object, ...)
  invisible(object)
}

trajectory_matrix_long <- function(values, time) {
  if (is.null(dim(values))) {
    values <- matrix(values, ncol = 1L)
  }
  if (!is.matrix(values) || nrow(values) != length(time)) {
    stop("through-time values are not aligned to tt$time", call. = FALSE)
  }
  labels <- colnames(values)
  if (is.null(labels)) labels <- paste0("series_", seq_len(ncol(values)))
  keep <- !grepl("(^time$|^times$|time_vec$)", labels)
  values <- values[, keep, drop = FALSE]
  labels <- labels[keep]
  data.frame(
    time = rep(as.numeric(time), times = ncol(values)),
    series = rep(labels, each = nrow(values)),
    value = as.numeric(values),
    stringsAsFactors = FALSE
  )
}

trajectory_find_interval_node <- function(node, path = character()) {
  if (is.list(node) && all(c("lower", "mean", "upper") %in% names(node)) &&
      is.numeric(node$lower) && is.numeric(node$mean) && is.numeric(node$upper)) {
    return(list(value = node, path = path))
  }
  if (!is.list(node)) return(NULL)
  for (field in names(node)) {
    found <- trajectory_find_interval_node(node[[field]], c(path, field))
    if (!is.null(found)) return(found)
  }
  NULL
}

trajectory_interval_long <- function(node) {
  if (!all(c("lower", "mean", "upper") %in% names(node))) {
    stop("selected distribution node lacks lower/mean/upper values", call. = FALSE)
  }
  template <- node$mean
  if (!identical(dim(node$lower), dim(template)) ||
      !identical(dim(node$upper), dim(template))) {
    stop("distribution interval arrays have different dimensions", call. = FALSE)
  }
  if (is.null(dim(template))) {
    labels <- names(template)
    if (is.null(labels)) labels <- paste0("value_", seq_along(template))
    return(data.frame(
      index = seq_along(template), series = labels,
      lower = as.numeric(node$lower), mean = as.numeric(template),
      upper = as.numeric(node$upper), stringsAsFactors = FALSE
    ))
  }
  matrix_mean <- matrix(template, nrow = dim(template)[[1L]])
  matrix_lower <- matrix(node$lower, nrow = dim(template)[[1L]])
  matrix_upper <- matrix(node$upper, nrow = dim(template)[[1L]])
  labels <- if (length(dim(template)) == 2L && !is.null(colnames(template))) {
    colnames(template)
  } else {
    paste0("series_", seq_len(ncol(matrix_mean)))
  }
  data.frame(
    index = rep(seq_len(nrow(matrix_mean)), times = ncol(matrix_mean)),
    series = rep(labels, each = nrow(matrix_mean)),
    lower = as.numeric(matrix_lower),
    mean = as.numeric(matrix_mean),
    upper = as.numeric(matrix_upper),
    stringsAsFactors = FALSE
  )
}

trajectory_tip_association <- function(left, right) {
  pairs <- lapply(left$tip.label, function(left_label) {
    matches <- right$tip.label[
      vapply(right$tip.label, function(right_label) {
        identical(left_label, right_label) ||
          grepl(left_label, right_label, fixed = TRUE) ||
          grepl(right_label, left_label, fixed = TRUE)
      }, logical(1))
    ]
    if (!length(matches)) return(NULL)
    cbind(left_label, matches)
  })
  pairs <- Filter(Negate(is.null), pairs)
  if (!length(pairs)) {
    stop("no tip associations could be inferred; supply association explicitly", call. = FALSE)
  }
  do.call(rbind, pairs)
}

`%||%` <- function(left, right) {
  if (is.null(left)) right else left
}
