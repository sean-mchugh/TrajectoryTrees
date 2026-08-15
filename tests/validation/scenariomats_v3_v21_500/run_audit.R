#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

package_root <- normalizePath(".", mustWork = TRUE)
contract_directory <- file.path(
  package_root,
  "tests", "validation", "scenariomats_v3_v21_500"
)
saved_tree_manifest <- file.path(contract_directory, "selected_inventory.csv")
tree_directory <- file.path(contract_directory, "trees")

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments)) {
  stop(
    "run_audit.R takes no trailing arguments; run it from the TrajectoryTrees package root",
    call. = FALSE
  )
}
if (!file.exists(file.path(package_root, "DESCRIPTION")) ||
    !dir.exists(contract_directory)) {
  stop(
    "Run this command from the TrajectoryTrees package root: ",
    "Rscript tests/validation/scenariomats_v3_v21_500/run_audit.R",
    call. = FALSE
  )
}

case_limit <- Sys.getenv(
  "TRAJECTORYTREES_SIMILARITY_AUDIT_CASE_LIMIT",
  unset = ""
)
case_limit <- if (nzchar(case_limit)) as.integer(case_limit) else Inf
if (length(case_limit) != 1L || is.na(case_limit) || case_limit < 1) {
  stop(
    "TRAJECTORYTREES_SIMILARITY_AUDIT_CASE_LIMIT must be one positive integer",
    call. = FALSE
  )
}

numeric_tolerance <- 1e-10
maximum_transition_counts <- 0:3
metrics <- c("bhattacharyya", "minimum")
return_modes <- c("complete", "matrices_and_summaries")

suppressPackageStartupMessages(library(TrajectoryTrees))
if (utils::packageVersion("TrajectoryTrees") < "0.2.0") {
  stop("TrajectoryTrees 0.2.0 or newer is required", call. = FALSE)
}
source(file.path(
  package_root,
  "tests", "testthat", "reference", "ScenarioMatsV2_1.R"
))

discrepancies <- list()
checks <- list()

common_v2_1_surface <- function(result) {
  result$summaries$similarity_vectors <- NULL
  if (!is.null(result$pair_details$interval_classification)) {
    rownames(result$pair_details$interval_classification) <- NULL
  }
  result
}

one_line <- function(value) {
  if (length(value) == 0L) return("<length 0>")
  if (length(value) > 1L) {
    return(paste(vapply(value, one_line, character(1)), collapse = " | "))
  }
  if (is.na(value)) return("<NA>")
  as.character(value)
}

add_discrepancy <- function(
    case_id,
    metric,
    return_mode,
    field,
    component = "",
    tip_pair = "",
    reference_value = "",
    candidate_value = "",
    difference = NA_real_,
    detail = "") {
  discrepancies[[length(discrepancies) + 1L]] <<- data.frame(
    case_id = case_id,
    tip_pair = tip_pair,
    metric = metric,
    return_mode = return_mode,
    field = field,
    component = component,
    v2_1_value = one_line(reference_value),
    v3_value = one_line(candidate_value),
    difference = difference,
    absolute_difference = abs(difference),
    detail = detail,
    stringsAsFactors = FALSE
  )
}

path_label <- function(path) {
  if (!length(path)) return("<root>")
  paste(path, collapse = "$")
}

compare_exact_structure <- function(
    reference,
    candidate,
    case_id,
    metric,
    return_mode,
    path = character(),
    component = "",
    tip_pair = "") {
  field <- path_label(path)
  reference_class <- class(reference)
  candidate_class <- class(candidate)
  if (!identical(reference_class, candidate_class)) {
    add_discrepancy(
      case_id, metric, return_mode, field, component, tip_pair,
      paste(reference_class, collapse = "/"),
      paste(candidate_class, collapse = "/"),
      detail = "class mismatch"
    )
  }

  if (is.data.frame(reference)) {
    if (!is.data.frame(candidate)) return(invisible(NULL))
    if (!identical(names(reference), names(candidate))) {
      add_discrepancy(
        case_id, metric, return_mode, field, component, tip_pair,
        names(reference), names(candidate), detail = "column names or order mismatch"
      )
    }
    if (nrow(reference) != nrow(candidate)) {
      add_discrepancy(
        case_id, metric, return_mode, field, component, tip_pair,
        nrow(reference), nrow(candidate),
        difference = nrow(candidate) - nrow(reference),
        detail = "row-count mismatch"
      )
    }
    common_columns <- intersect(names(reference), names(candidate))
    common_rows <- seq_len(min(nrow(reference), nrow(candidate)))
    for (column_name in common_columns) {
      for (row_index in common_rows) {
        row_tip_pair <- tip_pair
        if (all(c("tip_i", "tip_j") %in% names(reference))) {
          row_tip_pair <- paste(
            reference$tip_i[[row_index]],
            reference$tip_j[[row_index]],
            sep = " :: "
          )
        }
        compare_exact_structure(
          reference[[column_name]][row_index],
          candidate[[column_name]][row_index],
          case_id,
          metric,
          return_mode,
          c(path, paste0(column_name, "[", row_index, "]")),
          component,
          row_tip_pair
        )
      }
    }
    return(invisible(NULL))
  }

  if (is.matrix(reference) || is.array(reference)) {
    if (!is.matrix(candidate) && !is.array(candidate)) return(invisible(NULL))
    if (!identical(dim(reference), dim(candidate))) {
      add_discrepancy(
        case_id, metric, return_mode, field, component, tip_pair,
        dim(reference), dim(candidate), detail = "dimension mismatch"
      )
      return(invisible(NULL))
    }
    if (!identical(dimnames(reference), dimnames(candidate))) {
      add_discrepancy(
        case_id, metric, return_mode, field, component, tip_pair,
        unlist(dimnames(reference), use.names = FALSE),
        unlist(dimnames(candidate), use.names = FALSE),
        detail = "dimension labels mismatch"
      )
    }
    for (index in seq_along(reference)) {
      location <- arrayInd(index, dim(reference))
      matrix_tip_pair <- tip_pair
      if (length(dim(reference)) == 2L && !is.null(rownames(reference)) &&
          !is.null(colnames(reference))) {
        matrix_tip_pair <- paste(
          rownames(reference)[[location[[1L]]]],
          colnames(reference)[[location[[2L]]]],
          sep = " :: "
        )
      }
      compare_exact_structure(
        reference[[index]],
        candidate[[index]],
        case_id,
        metric,
        return_mode,
        c(path, paste0("[", paste(location, collapse = ","), "]")),
        component,
        matrix_tip_pair
      )
    }
    return(invisible(NULL))
  }

  if (is.list(reference)) {
    if (!is.list(candidate)) return(invisible(NULL))
    if (!identical(names(reference), names(candidate))) {
      add_discrepancy(
        case_id, metric, return_mode, field, component, tip_pair,
        names(reference), names(candidate), detail = "names or order mismatch"
      )
    }
    reference_names <- names(reference)
    candidate_names <- names(candidate)
    named_lists <- !is.null(reference_names) &&
      !is.null(candidate_names) &&
      all(nzchar(reference_names)) &&
      all(nzchar(candidate_names)) &&
      !anyDuplicated(reference_names) &&
      !anyDuplicated(candidate_names)
    if (named_lists) {
      child_names <- reference_names[reference_names %in% candidate_names]
      reference_indices <- match(child_names, reference_names)
      candidate_indices <- match(child_names, candidate_names)
    } else {
      common_length <- min(length(reference), length(candidate))
      child_names <- paste0("[[", seq_len(common_length), "]]")
      reference_indices <- seq_len(common_length)
      candidate_indices <- seq_len(common_length)
    }
    for (child_index in seq_along(child_names)) {
      child_name <- child_names[[child_index]]
      child_component <- component
      parent_is_matrix_family <- length(path) >= 1L && (
        tail(path, 1L) %in% c("sync", "state_only") ||
          grepl("^using_up_to_", tail(path, 1L))
      )
      if (parent_is_matrix_family) {
        child_component <- child_name
      }
      compare_exact_structure(
        reference[[reference_indices[[child_index]]]],
        candidate[[candidate_indices[[child_index]]]],
        case_id,
        metric,
        return_mode,
        c(path, child_name),
        child_component,
        tip_pair
      )
    }
    if (length(reference) != length(candidate)) {
      add_discrepancy(
        case_id, metric, return_mode, field, component, tip_pair,
        length(reference), length(candidate),
        difference = length(candidate) - length(reference),
        detail = "length mismatch"
      )
    }
    return(invisible(NULL))
  }

  if (is.numeric(reference) && is.numeric(candidate)) {
    if (length(reference) != length(candidate)) {
      add_discrepancy(
        case_id, metric, return_mode, field, component, tip_pair,
        length(reference), length(candidate),
        difference = length(candidate) - length(reference),
        detail = "numeric length mismatch"
      )
      return(invisible(NULL))
    }
    if (!identical(is.na(reference), is.na(candidate))) {
      add_discrepancy(
        case_id, metric, return_mode, field, component, tip_pair,
        reference, candidate, detail = "missingness mismatch"
      )
      return(invisible(NULL))
    }
    comparable <- !is.na(reference)
    for (index in which(comparable)) {
      difference <- candidate[[index]] - reference[[index]]
      if (!is.finite(reference[[index]]) || !is.finite(candidate[[index]]) ||
          abs(difference) > numeric_tolerance) {
        add_discrepancy(
          case_id, metric, return_mode, field, component, tip_pair,
          reference[[index]], candidate[[index]], difference,
          detail = "numeric mismatch"
        )
      }
    }
    if (!identical(names(reference), names(candidate))) {
      add_discrepancy(
        case_id, metric, return_mode, field, component, tip_pair,
        names(reference), names(candidate), detail = "numeric names mismatch"
      )
    }
    return(invisible(NULL))
  }

  if (!identical(reference, candidate)) {
    add_discrepancy(
      case_id, metric, return_mode, field, component, tip_pair,
      reference, candidate, detail = "exact-value mismatch"
    )
  }
  invisible(NULL)
}

add_check <- function(case_id, metric, return_mode, check_name, passed, detail = "") {
  checks[[length(checks) + 1L]] <<- data.frame(
    case_id = case_id,
    metric = metric,
    return_mode = return_mode,
    check = check_name,
    status = if (isTRUE(passed)) "pass" else "fail",
    detail = detail,
    stringsAsFactors = FALSE
  )
}

check_matrix_family <- function(
    family,
    tip_labels,
    case_id,
    metric,
    return_mode,
    family_name) {
  matrices_are_valid <- length(family) > 0L && all(vapply(
    family,
    function(value) {
      is.matrix(value) &&
        identical(rownames(value), tip_labels) &&
        identical(colnames(value), tip_labels) &&
        all(is.finite(value)) &&
        all(value >= -numeric_tolerance) &&
        isTRUE(all.equal(value, t(value), tolerance = numeric_tolerance)) &&
        all(abs(diag(value)) <= numeric_tolerance)
    },
    logical(1)
  ))
  add_check(
    case_id, metric, return_mode,
    paste0(family_name, "_matrix_invariants"),
    matrices_are_valid
  )
  if (!matrices_are_valid) return(invisible(NULL))
  closure <- Reduce(`+`, family)
  add_check(
    case_id, metric, return_mode,
    paste0(family_name, "_pairwise_closure"),
    all(abs(closure[upper.tri(closure)] - 1) <= numeric_tolerance)
  )
  invisible(NULL)
}

check_summary_against_matrices <- function(
    family,
    summary,
    available_similarity_by_pair,
    case_id,
    metric,
    return_mode,
    family_name) {
  upper_triangle <- upper.tri(family[[1L]])
  expected_totals <- vapply(
    family,
    function(component_matrix) {
      values <- component_matrix[upper_triangle]
      if (!is.null(available_similarity_by_pair)) {
        values <- values * available_similarity_by_pair[upper_triangle]
      }
      sum(values)
    },
    numeric(1)
  )
  expected_total <- sum(expected_totals)
  expected_proportions <- expected_totals / expected_total
  values_match <- function(actual, expected) {
    is.numeric(actual) &&
      length(actual) == length(expected) &&
      identical(names(actual), names(expected)) &&
      all(is.finite(actual)) &&
      all(abs(actual - expected) <= numeric_tolerance)
  }

  add_check(
    case_id, metric, return_mode,
    paste0(family_name, "_summary_totals_from_matrices"),
    values_match(summary$totals, expected_totals)
  )
  add_check(
    case_id, metric, return_mode,
    paste0(family_name, "_summary_proportions_from_matrices"),
    values_match(summary$tree_wide_proportions, expected_proportions)
  )
  add_check(
    case_id, metric, return_mode,
    paste0(family_name, "_summary_means_from_matrices"),
    values_match(summary$means, expected_proportions)
  )
  add_check(
    case_id, metric, return_mode,
    paste0(family_name, "_summary_total_across_matrices"),
    is.numeric(summary$total_across_matrices) &&
      length(summary$total_across_matrices) == 1L &&
      is.finite(summary$total_across_matrices) &&
      abs(summary$total_across_matrices - expected_total) <= numeric_tolerance
  )
  if (!is.null(available_similarity_by_pair)) {
    expected_available <- sum(available_similarity_by_pair[upper_triangle])
    add_check(
      case_id, metric, return_mode,
      paste0(family_name, "_summary_available_similarity_total"),
      is.numeric(summary$available_similarity_total) &&
        length(summary$available_similarity_total) == 1L &&
        is.finite(summary$available_similarity_total) &&
        abs(summary$available_similarity_total - expected_available) <=
          numeric_tolerance
    )
  }
  invisible(NULL)
}

check_result <- function(result, trajectory, case_id, metric, return_mode) {
  tip_labels <- as.character(trajectory$phylo$tip.label)
  weighting <- result$summaries$weighting
  weighted <- isTRUE(weighting$available_similarity_weighting_applied)
  sync_weights <- if (weighted) {
    weighting$synchronous_available_similarity_by_pair
  } else {
    NULL
  }
  async_weights <- if (weighted) {
    weighting$asynchronous_available_similarity_by_pair
  } else {
    NULL
  }
  check_matrix_family(
    result$matrices$sync,
    tip_labels,
    case_id,
    metric,
    return_mode,
    "sync"
  )
  check_summary_against_matrices(
    result$matrices$sync,
    result$summaries$sync,
    sync_weights,
    case_id,
    metric,
    return_mode,
    "sync"
  )
  for (view_name in names(result$matrices$async)) {
    check_matrix_family(
      result$matrices$async[[view_name]],
      tip_labels,
      case_id,
      metric,
      return_mode,
      paste0("async_", view_name)
    )
    check_summary_against_matrices(
      result$matrices$async[[view_name]],
      result$summaries$async[[view_name]],
      async_weights,
      case_id,
      metric,
      return_mode,
      paste0("async_", view_name)
    )
  }
  summary_families <- c(
    list(sync = result$summaries$sync),
    result$summaries$async
  )
  for (summary_name in names(summary_families)) {
    proportions <- summary_families[[summary_name]]$tree_wide_proportions
    add_check(
      case_id, metric, return_mode,
      paste0("summary_closure_", summary_name),
      is.numeric(proportions) &&
        all(is.finite(proportions)) &&
        abs(sum(proportions) - 1) <= numeric_tolerance
    )
  }
  if (!is.null(weighting)) {
    expected_weighting <- !isTRUE(weighting$tree_is_ultrametric)
    add_check(
      case_id, metric, return_mode,
      "available_similarity_weighting",
      identical(
        isTRUE(weighting$available_similarity_weighting_applied),
        expected_weighting
      )
    )
  } else {
    add_check(
      case_id, metric, return_mode,
      "available_similarity_weighting",
      FALSE,
      "summaries$weighting is missing"
    )
  }
}

load_saved_tree_cases <- function() {
  if (!file.exists(saved_tree_manifest)) {
    stop("Missing promoted 500-case manifest: ", saved_tree_manifest, call. = FALSE)
  }
  inventory <- utils::read.csv(saved_tree_manifest, stringsAsFactors = FALSE)
  if (nrow(inventory) != 500L) {
    stop("the saved-tree equality corpus must contain 500 rows", call. = FALSE)
  }
  if (anyNA(inventory$tip_count) ||
      any(inventory$tip_count < 3L | inventory$tip_count > 10L)) {
    stop("every saved-tree audit case must contain 3 to 10 tips", call. = FALSE)
  }
  object_files <- file.path(tree_directory, basename(inventory$tree_path))
  if (anyDuplicated(object_files) || any(!file.exists(object_files))) {
    stop(
      "the package-local 500-tree fixture corpus is incomplete or duplicated",
      call. = FALSE
    )
  }
  inventory <- head(inventory, case_limit)
  lapply(seq_len(nrow(inventory)), function(index) {
    selected <- inventory[index, , drop = FALSE]
    list(
      case_id = selected$saved_tree_id,
      object_file = file.path(tree_directory, basename(selected$tree_path)),
      corpus = "saved_tree"
    )
  })
}

trajectory_construction_tolerance <- function(phylo) {
  durations <- suppressWarnings(as.numeric(unlist(
    phylo$maps,
    recursive = TRUE,
    use.names = FALSE
  )))
  positive <- durations[is.finite(durations) & durations > 0]
  if (!length(positive)) {
    stop("saved-tree input has no finite positive map duration", call. = FALSE)
  }
  min(1e-10, min(positive) / 10)
}

cases <- load_saved_tree_cases()

for (case_index in seq_along(cases)) {
  case <- cases[[case_index]]
  message(sprintf("[%d/%d] %s", case_index, length(cases), case$case_id))
  saved_object <- readRDS(case$object_file)
  phylo <- if (inherits(saved_object, "phylo")) {
    saved_object
  } else if (is.list(saved_object) && inherits(saved_object$phylo, "phylo")) {
    saved_object$phylo
  } else {
    stop("saved-tree object does not contain a phylo: ", case$case_id)
  }
  if (length(phylo$tip.label) < 3L || length(phylo$tip.label) > 10L) {
    stop("saved-tree phylo is outside the 3-to-10-tip contract: ", case$case_id)
  }
  tolerance <- trajectory_construction_tolerance(phylo)
  trajectory <- TrajectoryTrees::make_trajectory_objects(
    phylo,
    include_tt = FALSE,
    include_ss = FALSE,
    include_scenario_mats = TRUE,
    include_path_maps = TRUE,
    time_tolerance = tolerance
  )
  for (metric in metrics) {
    results <- list()
    for (return_mode in return_modes) {
      reference <- trajectory_similarity_v2_1(
        trajectory,
        maximum_transition_counts = maximum_transition_counts,
        async_metric = metric,
        return_mode = return_mode
      )
      candidate <- TrajectoryTrees::trajectory_similarity(
        trajectory,
        maximum_transition_counts = maximum_transition_counts,
        async_metric = metric,
        return_mode = return_mode
      )
      results[[return_mode]] <- candidate
      compare_exact_structure(
        common_v2_1_surface(reference),
        common_v2_1_surface(candidate),
        case$case_id,
        metric,
        return_mode
      )
      check_result(
        candidate,
        trajectory,
        case$case_id,
        metric,
        return_mode
      )
    }
    compare_exact_structure(
      results$complete[c("matrices", "summaries")],
      results$matrices_and_summaries,
      case$case_id,
      metric,
      "complete_vs_compact"
    )
  }
}

discrepancy_table <- if (length(discrepancies)) {
  do.call(rbind, discrepancies)
} else {
  data.frame(
    case_id = character(),
    tip_pair = character(),
    metric = character(),
    return_mode = character(),
    field = character(),
    component = character(),
    v2_1_value = character(),
    v3_value = character(),
    difference = numeric(),
    absolute_difference = numeric(),
    detail = character(),
    stringsAsFactors = FALSE
  )
}
check_table <- do.call(rbind, checks)
result_root <- file.path(contract_directory, "results")
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_directory <- file.path(
  result_root,
  sprintf(
    "audit_%03d_of_500_%s_pid%d",
    length(cases),
    run_stamp,
    Sys.getpid()
  )
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  discrepancy_table,
  file.path(output_directory, "discrepancies.csv"),
  row.names = FALSE
)
utils::write.csv(
  check_table,
  file.path(output_directory, "internal_checks.csv"),
  row.names = FALSE
)

failed_checks <- check_table[check_table$status == "fail", , drop = FALSE]
audit_status <- if (nrow(discrepancy_table) || nrow(failed_checks)) {
  "fail"
} else {
  "pass"
}
utils::write.csv(
  data.frame(
    status = audit_status,
    histories = length(cases),
    corpus_size = 500L,
    package_version = as.character(utils::packageVersion("TrajectoryTrees")),
    candidate = "TrajectoryTrees::trajectory_similarity",
    oracle = "tests/testthat/reference/ScenarioMatsV2_1.R",
    maximum_transition_counts = paste(maximum_transition_counts, collapse = ","),
    async_metrics = paste(metrics, collapse = ","),
    return_modes = paste(return_modes, collapse = ","),
    numeric_tolerance = numeric_tolerance,
    run_started_from = package_root,
    stringsAsFactors = FALSE
  ),
  file.path(output_directory, "audit_metadata.csv"),
  row.names = FALSE
)
if (nrow(discrepancy_table) || nrow(failed_checks)) {
  stop(
    sprintf(
      paste(
        "ScenarioMats V3 equality failed with %d discrepancies and",
        "%d failed internal checks. Results: %s"
      ),
      nrow(discrepancy_table),
      nrow(failed_checks),
      output_directory
    ),
    call. = FALSE
  )
}

message(sprintf(
  "TrajectoryTrees similarity equality passed for %d histories. Results: %s",
  length(cases),
  output_directory
))
