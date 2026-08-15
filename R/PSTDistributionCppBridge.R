# R bridge for the exact downstream distribution reducer registered in the
# TrajectoryTrees package DLL.

#' Convert a byte target to a safe positive C++ element batch size.
#'
#' @param block_bytes Positive requested byte target.
#' @param sample_count Positive posterior sample count.
#' @return Positive integer no larger than `.Machine$integer.max`.
#' @noRd
pst_distribution_cpp_elements_per_block <- function(block_bytes, sample_count) {
  raw_elements <- max(1, floor(block_bytes / (8 * sample_count)))
  as.integer(min(raw_elements, .Machine$integer.max))
}

#' Verify that the exact distribution native routines are loaded.
#'
#' @return Invisibly returns `TRUE` after the exported symbol is available.
#' @noRd
pst_distribution_load_cpp <- function() {
  if (!is.loaded("_TrajectoryTrees_pst_distribution_cpp_reduce_exact")) {
    stop(
      "TrajectoryTrees native distribution routines are not loaded; reinstall the package",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Reduce equal-shaped numeric arrays through bounded exact C++ blocks.
#'
#' @param arrays Numeric arrays, one per posterior sample.
#' @param probs Lower and upper type-7 probabilities.
#' @param block_bytes Positive working-buffer target in bytes.
#' @param path Diagnostic surface path.
#' @return Lower, mean, median, upper, and observation-count arrays.
#' @noRd
pst_distribution_reduce_arrays_cpp <- function(
    arrays,
    probs = c(0.025, 0.975),
    block_bytes = 16 * 1024^2,
    path = "distribution") {
  if (!length(arrays)) {
    stop(path, " requires at least one posterior sample", call. = FALSE)
  }
  template <- arrays[[1L]]
  if (any(vapply(arrays, function(value) !identical(dim(value), dim(template)), logical(1)))) {
    stop(path, " array dimensions do not match", call. = FALSE)
  }
  if (length(block_bytes) != 1L || !is.finite(block_bytes) || block_bytes <= 0) {
    stop("block_bytes must be one positive finite value", call. = FALSE)
  }
  pst_distribution_load_cpp()
  sample_count <- length(arrays)
  element_count <- length(template)
  elements_per_block <- pst_distribution_cpp_elements_per_block(
    block_bytes, sample_count
  )
  flattened <- lapply(arrays, as.numeric)
  fields <- c("lower", "mean", "median", "upper")
  output <- setNames(lapply(fields, function(field) numeric(element_count)), fields)
  output$n_observed <- integer(element_count)
  # Fill and reduce one bounded samples-by-elements block at a time.
  for (block_start in seq.int(1L, element_count, by = elements_per_block)) {
    block_end <- min(element_count, block_start + elements_per_block - 1L)
    indices <- block_start:block_end
    values <- matrix(NA_real_, nrow = sample_count, ncol = length(indices))
    # Copy each sample's contiguous element slice into the reusable block shape.
    for (sample_index in seq_len(sample_count)) {
      values[sample_index, ] <- flattened[[sample_index]][indices]
    }
    reduced <- tryCatch(
      pst_distribution_cpp_reduce_exact(values, probs),
      error = function(error) {
        stop(path, " block ", block_start, "-", block_end, ": ",
             conditionMessage(error), call. = FALSE)
      }
    )
    for (field in fields) {
      output[[field]][indices] <- reduced[[field]]
    }
    output$n_observed[indices] <- reduced$n_observed
  }
  # Restore the canonical dimensions and labels after all blocks are reduced.
  for (field in names(output)) {
    if (is.null(dim(template))) {
      names(output[[field]]) <- names(template)
    } else {
      output[[field]] <- array(output[[field]], dim = dim(template),
                               dimnames = dimnames(template))
    }
  }
  storage.mode(output$n_observed) <- "integer"
  output
}

#' Fuse TT step alignment and exact reduction in bounded C++ blocks.
#'
#' @param matrices Numeric TT matrices, one per sample.
#' @param source_times Per-sample time vectors.
#' @param canonical_time Canonical union time vector.
#' @param columns Canonical semantic columns.
#' @param probs Lower and upper type-7 probabilities.
#' @param tolerance Time comparison tolerance.
#' @param block_bytes Working-buffer target in bytes.
#' @param path Diagnostic surface path.
#' @return Four canonical matrices plus `n_observed`.
#' @noRd
pst_distribution_align_reduce_tt_cpp <- function(
    matrices,
    source_times,
    canonical_time,
    columns,
    probs,
    tolerance,
    block_bytes,
    path,
    alignment_context = NULL) {
  pst_distribution_load_cpp()
  sample_count <- length(matrices)
  ordered_matrices <- vector("list", sample_count)
  row_maps <- vector("list", sample_count)
  column_maps <- vector("list", sample_count)
  # Precompute each sample's sorted matrix and canonical step-row/column maps once.
  for (sample_index in seq_len(sample_count)) {
    order_index <- if (is.null(alignment_context)) {
      order(source_times[[sample_index]])
    } else {
      alignment_context$order_indices[[sample_index]]
    }
    # Preserve the original matrix SEXP so repeated trajectory objects can be
    # collapsed to weighted representatives inside C++.
    ordered_matrices[[sample_index]] <- matrices[[sample_index]]
    row_maps[[sample_index]] <- if (is.null(alignment_context)) {
      sorted_time <- source_times[[sample_index]][order_index]
      interval <- findInterval(canonical_time + tolerance, sorted_time)
      if (any(interval < 1L)) stop(path, " precedes a sample time range", call. = FALSE)
      order_index[interval] - 1L
    } else {
      alignment_context$row_maps[[sample_index]]
    }
    matched <- match(columns, colnames(matrices[[sample_index]]))
    column_maps[[sample_index]] <- ifelse(is.na(matched), -1L, matched - 1L)
  }
  structural_time <- grepl("(^time$|^times$|time_vec$)", columns)
  elements_per_block <- pst_distribution_cpp_elements_per_block(
    block_bytes, sample_count
  )
  reduced <- tryCatch(
    pst_distribution_cpp_align_reduce_tt(
      ordered_matrices, row_maps, column_maps, structural_time,
      canonical_time, probs, elements_per_block
    ),
    error = function(error) {
      stop(path, ": ", conditionMessage(error), call. = FALSE)
    }
  )
  # Restore canonical matrix dimensions and labels for every returned field.
  for (field in names(reduced)) {
    reduced[[field]] <- matrix(
      reduced[[field]], nrow = length(canonical_time), ncol = length(columns),
      dimnames = list(NULL, columns)
    )
  }
  storage.mode(reduced$n_observed) <- "integer"
  reduced
}
