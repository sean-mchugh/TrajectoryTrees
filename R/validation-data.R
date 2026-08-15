#' Load the bundled 626-tree validation corpus
#'
#' The corpus contains the 471 established ultrametric cases and 155
#' non-ultrametric cases used by the release validation. Filters are applied to
#' the manifest before any tree files are read.
#'
#' @param population `"all"`, `"ultrametric"`, or `"nonultrametric"`.
#' @param case_ids Optional exact case identifiers to include.
#' @param include_label,exclude_label Optional regular expressions applied to
#'   case identifiers.
#' @param min_tips,max_tips Optional inclusive tip-count bounds.
#' @param exclude_tips Optional vector of exact tip counts to exclude.
#' @param collection Optional collection names such as `"adversarial"` or
#'   `"generated"`.
#' @param start One-based offset after filtering.
#' @param limit Maximum number of cases to return; `Inf` keeps all.
#' @param load Whether to read the selected RDS files. `FALSE` returns the
#'   filtered inventory.
#'
#' @return A named list of `simmap` trees, or a
#'   `trajectory_tree_inventory` data frame when `load = FALSE`.
#' @export
load_validation_trees <- function(
    population = c("all", "ultrametric", "nonultrametric"),
    case_ids = NULL,
    include_label = NULL,
    exclude_label = NULL,
    min_tips = NULL,
    max_tips = NULL,
    exclude_tips = NULL,
    collection = NULL,
    start = 1L,
    limit = Inf,
    load = TRUE) {
  population <- match.arg(population)
  data_dir <- system.file(
    "extdata", "validation",
    package = "TrajectoryTrees",
    mustWork = TRUE
  )
  inventory <- utils::read.csv(
    file.path(data_dir, "corpus_manifest.csv"),
    stringsAsFactors = FALSE
  )
  if (identical(population, "ultrametric")) {
    inventory <- inventory[inventory$ultrametric, , drop = FALSE]
  } else if (identical(population, "nonultrametric")) {
    inventory <- inventory[!inventory$ultrametric, , drop = FALSE]
  }
  if (!is.null(case_ids)) {
    inventory <- inventory[inventory$case_id %in% case_ids, , drop = FALSE]
  }
  if (!is.null(include_label)) {
    inventory <- inventory[
      grepl(include_label, inventory$case_id, perl = TRUE),
      , drop = FALSE
    ]
  }
  if (!is.null(exclude_label)) {
    inventory <- inventory[
      !grepl(exclude_label, inventory$case_id, perl = TRUE),
      , drop = FALSE
    ]
  }
  if (!is.null(min_tips)) {
    inventory <- inventory[inventory$n_tip >= min_tips, , drop = FALSE]
  }
  if (!is.null(max_tips)) {
    inventory <- inventory[inventory$n_tip <= max_tips, , drop = FALSE]
  }
  if (!is.null(exclude_tips)) {
    inventory <- inventory[!inventory$n_tip %in% exclude_tips, , drop = FALSE]
  }
  if (!is.null(collection)) {
    inventory <- inventory[inventory$collection %in% collection, , drop = FALSE]
  }
  if (!is.numeric(start) || length(start) != 1L || !is.finite(start) ||
      start < 1 || start != as.integer(start)) {
    stop("start must be one positive integer", call. = FALSE)
  }
  if (!is.numeric(limit) || length(limit) != 1L || is.na(limit) ||
      limit < 0 || (!is.infinite(limit) && limit != as.integer(limit))) {
    stop("limit must be a nonnegative integer or Inf", call. = FALSE)
  }
  if (!is.logical(load) || length(load) != 1L || is.na(load)) {
    stop("load must be TRUE or FALSE", call. = FALSE)
  }

  if (nrow(inventory) < start || identical(limit, 0L)) {
    inventory <- inventory[0, , drop = FALSE]
  } else {
    last <- if (is.infinite(limit)) {
      nrow(inventory)
    } else {
      min(nrow(inventory), as.integer(start + limit - 1L))
    }
    inventory <- inventory[seq.int(as.integer(start), last), , drop = FALSE]
  }
  rownames(inventory) <- NULL
  if (!isTRUE(load)) {
    class(inventory) <- unique(c("trajectory_tree_inventory", class(inventory)))
    return(inventory)
  }
  paths <- file.path(data_dir, inventory$tree_path)
  trees <- lapply(paths, readRDS)
  names(trees) <- inventory$case_id
  attr(trees, "inventory") <- inventory
  trees
}

#' Load bundled Anolis stochastic-map ensembles
#'
#' @param set `"empirical"` for 500 maps or `"simulated"` for 1,000 maps.
#' @param start One-based starting index.
#' @param limit Maximum maps to load; `Inf` loads the remainder.
#'
#' @return A named list of `simmap` objects.
#' @export
load_anolis_trees <- function(
    set = c("empirical", "simulated"),
    start = 1L,
    limit = Inf) {
  set <- match.arg(set)
  filename <- switch(
    set,
    empirical = "anolis_empirical_500.rds",
    simulated = "anolis_simulated_1000.rds"
  )
  path <- system.file(
    "extdata", "anolis", filename,
    package = "TrajectoryTrees",
    mustWork = TRUE
  )
  trees <- readRDS(path)
  if (!is.numeric(start) || length(start) != 1L || !is.finite(start) ||
      start < 1 || start != as.integer(start)) {
    stop("start must be one positive integer", call. = FALSE)
  }
  if (!is.numeric(limit) || length(limit) != 1L || is.na(limit) ||
      limit < 0 || (!is.infinite(limit) && limit != as.integer(limit))) {
    stop("limit must be a nonnegative integer or Inf", call. = FALSE)
  }
  if (start > length(trees) || identical(limit, 0L)) {
    return(trees[0])
  }
  last <- if (is.infinite(limit)) {
    length(trees)
  } else {
    min(length(trees), as.integer(start + limit - 1L))
  }
  trees[seq.int(as.integer(start), last)]
}
