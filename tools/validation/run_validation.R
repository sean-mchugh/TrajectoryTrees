#!/usr/bin/env Rscript

# Standalone release validation for the public TrajectoryTrees package.
#
# Normal package users never load V34. This runner loads it only into a private
# environment when --compare-v34=true is requested.

validation_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (!length(file_arg)) {
    stop("run_validation.R must be executed with Rscript", call. = FALSE)
  }
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
}

validation_parse_args <- function(args) {
  defaults <- list(
    scopes = "corpus,known-true,anolis-empirical,anolis-simulated",
    population = "all",
    case_ids = "",
    include_label = "",
    exclude_label = "",
    min_tips = "",
    max_tips = "",
    exclude_tips = "",
    start = "1",
    limit = "Inf",
    compare_v34 = "false",
    tolerance = "1e-8",
    output = ""
  )
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) {
      stop("arguments must use --name=value: ", arg, call. = FALSE)
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", pieces[[1L]], fixed = TRUE)
    if (!key %in% names(defaults)) {
      stop("unknown argument: --", pieces[[1L]], call. = FALSE)
    }
    defaults[[key]] <- paste(pieces[-1L], collapse = "=")
  }
  defaults
}

validation_bool <- function(value, name) {
  normalized <- tolower(value)
  if (!normalized %in% c("true", "false")) {
    stop("--", name, " must be true or false", call. = FALSE)
  }
  identical(normalized, "true")
}

validation_optional_number <- function(value) {
  if (!nzchar(value)) return(NULL)
  number <- suppressWarnings(as.numeric(value))
  if (length(number) != 1L || is.na(number)) {
    stop("invalid numeric filter: ", value, call. = FALSE)
  }
  number
}

validation_csv <- function(value, mode = "character") {
  if (!nzchar(value)) return(NULL)
  result <- strsplit(value, ",", fixed = TRUE)[[1L]]
  if (identical(mode, "numeric")) {
    result <- suppressWarnings(as.numeric(result))
    if (anyNA(result)) stop("invalid numeric list: ", value, call. = FALSE)
  }
  result
}

validation_filter_inventory <- function(inventory, options) {
  if (!identical(options$population, "all")) {
    if (!options$population %in% c("ultrametric", "nonultrametric")) {
      stop("--population must be all, ultrametric, or nonultrametric", call. = FALSE)
    }
    keep_ultrametric <- identical(options$population, "ultrametric")
    inventory <- inventory[inventory$ultrametric == keep_ultrametric, , drop = FALSE]
  }
  case_ids <- validation_csv(options$case_ids)
  if (!is.null(case_ids)) {
    inventory <- inventory[inventory$case_id %in% case_ids, , drop = FALSE]
  }
  if (nzchar(options$include_label)) {
    inventory <- inventory[
      grepl(options$include_label, inventory$case_id, perl = TRUE),
      , drop = FALSE
    ]
  }
  if (nzchar(options$exclude_label)) {
    inventory <- inventory[
      !grepl(options$exclude_label, inventory$case_id, perl = TRUE),
      , drop = FALSE
    ]
  }
  min_tips <- validation_optional_number(options$min_tips)
  max_tips <- validation_optional_number(options$max_tips)
  excluded <- validation_csv(options$exclude_tips, "numeric")
  if (!is.null(min_tips)) {
    inventory <- inventory[inventory$n_tip >= min_tips, , drop = FALSE]
  }
  if (!is.null(max_tips)) {
    inventory <- inventory[inventory$n_tip <= max_tips, , drop = FALSE]
  }
  if (!is.null(excluded)) {
    inventory <- inventory[!inventory$n_tip %in% excluded, , drop = FALSE]
  }
  start <- as.integer(options$start)
  limit <- suppressWarnings(as.numeric(options$limit))
  if (is.na(start) || start < 1L || is.na(limit) || limit < 0) {
    stop("--start and --limit must be nonnegative valid indices", call. = FALSE)
  }
  if (!nrow(inventory) || start > nrow(inventory) || identical(limit, 0)) {
    return(inventory[0, , drop = FALSE])
  }
  last <- if (is.infinite(limit)) {
    nrow(inventory)
  } else {
    min(nrow(inventory), start + as.integer(limit) - 1L)
  }
  inventory[seq.int(start, last), , drop = FALSE]
}

validation_load_v34 <- function(validation_root) {
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("hidden V34 comparison requires Rcpp", call. = FALSE)
  }
  v34_root <- file.path(validation_root, "reference", "v34")
  source_dir <- file.path(v34_root, "src", "versions", "v34")
  environment <- new.env(parent = globalenv())
  environment$pst_project_root <- local({
    root <- v34_root
    function() root
  })
  environment$pst_observed_states <- function(simmap) {
    unique(unlist(lapply(simmap$maps, names), use.names = FALSE))
  }
  environment$pst_root_edge_indices <- function(simmap) {
    root <- length(simmap$tip.label) + 1L
    which(simmap$edge[, 1L] == root)
  }
  for (path in sort(list.files(source_dir, pattern = "[.]R$", full.names = TRUE))) {
    sys.source(path, envir = environment)
  }
  environment
}

validation_compare <- function(reference, candidate, tolerance) {
  comparison <- all.equal(
    reference, candidate,
    tolerance = tolerance,
    check.attributes = FALSE
  )
  if (isTRUE(comparison)) "" else paste(comparison, collapse = " | ")
}

validation_project_legacy_candidate <- function(candidate) {
  projected <- candidate
  projected$scenario$terminal_events <- NULL
  projected$scenario$edge_type <- NULL
  projected$scenario$path_maps <- NULL
  projected$phylo$path_maps <- NULL
  projected$scenario_mats_active <- NULL
  projected$scenario_mats$scenario_edge_step_ids <- NULL
  projected$scenario_mats$phylo_edge_step_ids <- NULL
  projected$scenario_mats$scenario_paths <- NULL
  projected$scenario_mats$phylo_paths <- NULL
  projected$tt$cltt <- NULL
  projected$tt$cstt <- NULL
  projected$tt$event_family <- NULL
  projected$tt$metadata <- NULL
  projected$ss$summary <- NULL
  if (is.list(projected$ss$metadata)) {
    projected$ss$metadata$construction_source <- NULL
  }
  projected
}

validation_project_legacy_reference <- function(reference) {
  projected <- reference
  projected$tt$event_family <- NULL
  projected$tt$metadata <- NULL
  projected$ss$summary <- NULL
  projected
}

validation_construct <- function(tree) {
  TrajectoryTrees::make_trajectory_objects(
    tree,
    include_tt = TRUE,
    include_ss = TRUE,
    include_scenario_mats = TRUE,
    include_path_maps = TRUE,
    fields = NULL,
    include_legacy_summary = TRUE
  )
}

validation_construct_v34 <- function(environment, tree, tolerance) {
  environment$get_trajectory_obj_v34(
    tree,
    include_tt = TRUE,
    include_ss = TRUE,
    include_legacy_summary = TRUE,
    include_scenario_mats = TRUE,
    include_path_maps = TRUE,
    time_tolerance = min(
      1e-10,
      min(unlist(tree$maps, use.names = FALSE)[
        unlist(tree$maps, use.names = FALSE) > 0
      ]) / 10
    ),
    run_tests = FALSE
  )
}

validation_one_case <- function(
    case_id,
    scope,
    tree,
    tolerance,
    v34_environment = NULL,
    known_reference = NULL,
    reference_flags = NULL) {
  before <- digest::digest(tree, algo = "xxhash64")
  started <- proc.time()[["elapsed"]]
  error_message <- ""
  comparison <- ""
  candidate <- tryCatch(
    validation_construct(tree),
    error = function(error) {
      error_message <<- conditionMessage(error)
      NULL
    }
  )
  modified <- !identical(before, digest::digest(tree, algo = "xxhash64"))
  if (!is.null(candidate) && !is.null(v34_environment)) {
    baseline <- tryCatch(
      validation_construct_v34(v34_environment, tree, tolerance),
      error = function(error) {
        error_message <<- paste0("V34: ", conditionMessage(error))
        NULL
      }
    )
    if (!is.null(baseline)) {
      comparison <- validation_compare(baseline, unclass(candidate), tolerance)
    }
  }
  if (!is.null(candidate) && !is.null(known_reference)) {
    fields <- character()
    if (isTRUE(reference_flags$check_pst)) {
      fields <- c(fields, "phylo", "scenario", "trans", "root_policy")
    }
    if (isTRUE(reference_flags$check_ss) ||
        isTRUE(reference_flags$check_ss_original_metric)) {
      fields <- c(fields, "ss")
    }
    if (isTRUE(reference_flags$check_tt)) {
      fields <- c(fields, "tt")
    }
    fields <- intersect(unique(fields), intersect(
      names(known_reference), names(candidate)
    ))
    known_reference <- validation_project_legacy_reference(known_reference)
    candidate <- validation_project_legacy_candidate(unclass(candidate))
    comparison <- validation_compare(
      known_reference[fields], candidate[fields], tolerance
    )
  }
  data.frame(
    scope = scope,
    case_id = case_id,
    n_tip = length(tree$tip.label),
    constructed = !is.null(candidate),
    input_modified = modified,
    comparison_passed = !nzchar(comparison),
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    error = error_message,
    comparison = comparison,
    stringsAsFactors = FALSE
  )
}

script_path <- validation_script_path()
validation_root <- dirname(script_path)
package_root <- normalizePath(file.path(validation_root, "..", ".."), mustWork = TRUE)
options <- validation_parse_args(commandArgs(trailingOnly = TRUE))
scopes <- validation_csv(options$scopes)
allowed_scopes <- c(
  "corpus", "known-true", "anolis-empirical", "anolis-simulated"
)
if (!length(scopes) || any(!scopes %in% allowed_scopes)) {
  stop("--scopes may contain only: ", paste(allowed_scopes, collapse = ", "), call. = FALSE)
}
tolerance <- as.numeric(options$tolerance)
if (!is.finite(tolerance) || tolerance < 0) {
  stop("--tolerance must be one nonnegative finite number", call. = FALSE)
}
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("validation requires the pkgload package", call. = FALSE)
}
pkgload::load_all(package_root, quiet = TRUE)
compare_v34 <- validation_bool(options$compare_v34, "compare-v34")
v34_environment <- if (compare_v34) {
  validation_load_v34(validation_root)
} else {
  NULL
}
rows <- list()

if ("corpus" %in% scopes) {
  inventory_path <- file.path(
    package_root, "inst", "extdata", "validation", "corpus_manifest.csv"
  )
  inventory <- validation_filter_inventory(
    utils::read.csv(inventory_path, stringsAsFactors = FALSE),
    options
  )
  for (index in seq_len(nrow(inventory))) {
    tree <- readRDS(file.path(
      dirname(inventory_path), inventory$tree_path[[index]]
    ))
    rows[[length(rows) + 1L]] <- validation_one_case(
      inventory$case_id[[index]], "corpus", tree, tolerance, v34_environment
    )
  }
}

if ("known-true" %in% scopes) {
  true_root <- file.path(validation_root, "reference", "true")
  manifest <- utils::read.csv(
    file.path(true_root, "manifest_complete.csv"),
    stringsAsFactors = FALSE
  )
  manifest$n_tip <- vapply(seq_len(nrow(manifest)), function(index) {
    reference <- readRDS(file.path(true_root, manifest$reference_path[[index]]))
    length(reference$phylo$tip.label)
  }, integer(1))
  manifest$ultrametric <- vapply(seq_len(nrow(manifest)), function(index) {
    reference <- readRDS(file.path(true_root, manifest$reference_path[[index]]))
    isTRUE(ape::is.ultrametric(reference$phylo))
  }, logical(1))
  manifest <- validation_filter_inventory(manifest, options)
  for (index in seq_len(nrow(manifest))) {
    reference <- readRDS(file.path(true_root, manifest$reference_path[[index]]))
    rows[[length(rows) + 1L]] <- validation_one_case(
      manifest$case_id[[index]], "known-true", reference$phylo, tolerance,
      known_reference = reference,
      reference_flags = manifest[index, , drop = FALSE]
    )
  }
}

for (scope in intersect(scopes, c("anolis-empirical", "anolis-simulated"))) {
  set <- if (identical(scope, "anolis-empirical")) "empirical" else "simulated"
  trees <- TrajectoryTrees::load_anolis_trees(set)
  inventory <- data.frame(
    case_id = names(trees),
    n_tip = vapply(trees, function(tree) length(tree$tip.label), integer(1)),
    ultrametric = vapply(trees, ape::is.ultrametric, logical(1)),
    stringsAsFactors = FALSE
  )
  inventory <- validation_filter_inventory(inventory, options)
  for (index in seq_len(nrow(inventory))) {
    case_id <- inventory$case_id[[index]]
    rows[[length(rows) + 1L]] <- validation_one_case(
      case_id, scope, trees[[case_id]], tolerance, v34_environment
    )
  }
}

results <- if (length(rows)) do.call(rbind, rows) else data.frame()
passed <- nrow(results) > 0L &&
  all(results$constructed) &&
  !any(results$input_modified) &&
  all(results$comparison_passed) &&
  !any(nzchar(results$error))
summary <- data.frame(
  status = if (passed) "pass" else "fail",
  cases = nrow(results),
  construction_failures = if (nrow(results)) sum(!results$constructed) else 0L,
  input_mutations = if (nrow(results)) sum(results$input_modified) else 0L,
  comparison_failures = if (nrow(results)) sum(!results$comparison_passed) else 0L,
  compare_v34 = compare_v34,
  tolerance = tolerance,
  stringsAsFactors = FALSE
)
if (nzchar(options$output)) {
  output <- normalizePath(options$output, mustWork = FALSE)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(results, file.path(output, "case_results.csv"), row.names = FALSE)
  utils::write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
}
print(summary)
if (!passed) quit(status = 1L)
