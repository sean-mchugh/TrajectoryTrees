# Recursive SS reduction and step-aligned TT reduction.

#' Reduce aligned numeric arrays elementwise.
#'
#' @param arrays List of equal-dimension numeric arrays.
#' @param probs Lower and upper quantile probabilities.
#' @return Four arrays preserving the first input's dimensions and labels.
#' @noRd
pst_distribution_reduce_arrays <- function(arrays, probs, backend = "R",
                                           block_bytes = 16 * 1024^2,
                                           path = "distribution") {
  # Route whole numeric leaves to C++ only after R has established canonical shape.
  if (identical(backend, "cpp")) {
    return(pst_distribution_reduce_arrays_cpp(arrays, probs, block_bytes, path))
  }
  template <- arrays[[1L]]
  stacked <- do.call(cbind, lapply(arrays, as.numeric))
  summaries <- lapply(seq_len(nrow(stacked)), function(index) {
    pst_distribution_reduce_numeric(stacked[index, ], probs)
  })
  # Reconstruct each statistic with the exact canonical dimensions and labels.
  out <- lapply(c("lower", "mean", "median", "upper"), function(field) {
    values <- vapply(summaries, `[[`, numeric(1), field)
    array(values, dim = dim(template), dimnames = dimnames(template))
  })
  names(out) <- c("lower", "mean", "median", "upper")
  observed <- rowSums(!is.na(stacked))
  out$n_observed <- array(as.integer(observed), dim = dim(template),
                          dimnames = dimnames(template))
  out
}

#' Union keyed SS data frames and reduce their numeric cells.
#'
#' @param frames Parallel SS data frames.
#' @param probs Quantile probabilities.
#' @param path Diagnostic path.
#' @return Four keyed summary frames plus numeric observation counts.
#' @noRd
pst_distribution_reduce_data_frames <- function(frames, probs, path,
                                                backend = "R",
                                                block_bytes = 16 * 1024^2) {
  key_columns <- sort(unique(unlist(lapply(frames, function(frame) {
    names(frame)[!vapply(frame, is.numeric, logical(1))]
  }), use.names = FALSE)))
  numeric_columns <- sort(unique(unlist(lapply(frames, function(frame) {
    names(frame)[vapply(frame, is.numeric, logical(1))]
  }), use.names = FALSE)))
  if (!length(key_columns) || !length(numeric_columns)) {
    stop(path, " data frames require semantic key and numeric columns", call. = FALSE)
  }
  # Form stable composite row keys from every nonnumeric semantic identifier.
  keyed <- lapply(frames, function(frame) {
    missing_keys <- setdiff(key_columns, names(frame))
    if (length(missing_keys)) stop(path, " data-frame keys do not match", call. = FALSE)
    keys <- do.call(paste, c(lapply(frame[key_columns], as.character), sep = "\037"))
    if (anyDuplicated(keys)) stop(path, " contains duplicate semantic rows", call. = FALSE)
    list(frame = frame, keys = keys)
  })
  union_keys <- sort(unique(unlist(lapply(keyed, `[[`, "keys"), use.names = FALSE)))
  key_source <- do.call(rbind, lapply(keyed, function(item) {
    data.frame(.key = item$keys, item$frame[key_columns], stringsAsFactors = FALSE)
  }))
  key_source <- key_source[!duplicated(key_source$.key), , drop = FALSE]
  key_values <- key_source[match(union_keys, key_source$.key), key_columns, drop = FALSE]
  arrays <- vector("list", length(frames))
  # Align every sample onto union rows/columns; absent rows and columns stay zero.
  for (sample_index in seq_along(frames)) {
    frame <- keyed[[sample_index]]$frame
    values <- matrix(0, nrow = length(union_keys), ncol = length(numeric_columns),
                     dimnames = list(union_keys, numeric_columns))
    present_numeric <- intersect(numeric_columns, names(frame))
    values[keyed[[sample_index]]$keys, present_numeric] <- as.matrix(frame[present_numeric])
    arrays[[sample_index]] <- values
  }
  reduced <- pst_distribution_reduce_arrays(
    arrays, probs, backend = backend, block_bytes = block_bytes, path = path
  )
  out <- lapply(reduced[c("lower", "mean", "median", "upper")], function(values) {
    data.frame(key_values, setNames(as.data.frame(values), numeric_columns),
               check.names = FALSE, stringsAsFactors = FALSE)
  })
  out$n_observed <- data.frame(key_values,
                               setNames(as.data.frame(reduced$n_observed), numeric_columns),
                               check.names = FALSE, stringsAsFactors = FALSE)
  out
}

#' Recursively union and reduce an SS node.
#'
#' @param nodes Parallel SS nodes from all samples.
#' @param probs Quantile probabilities.
#' @param path Diagnostic object path.
#' @return Preserved structure with numeric leaves replaced by four summaries.
#' @noRd
pst_distribution_reduce_ss_node <- function(nodes, probs, path = "SS",
                                            backend = "R",
                                            block_bytes = 16 * 1024^2) {
  non_null <- nodes[!vapply(nodes, is.null, logical(1))]
  # A field absent from every sample carries no distribution information.
  if (!length(non_null)) return(NULL)
  first <- non_null[[1L]]
  # Missing numeric/list categories contribute zero or an empty semantic node.
  if (is.data.frame(first)) {
    nodes <- lapply(nodes, function(node) if (is.null(node)) first[0, , drop = FALSE] else node)
  } else if (is.list(first)) {
    nodes <- lapply(nodes, function(node) if (is.null(node)) list() else node)
  } else if (is.numeric(first)) {
    nodes <- lapply(nodes, function(node) {
      if (!is.null(node)) return(node)
      if (is.null(dim(first))) {
        if (is.null(names(first))) return(0)
        return(setNames(numeric(0), character(0)))
      }
      array(0, dim = dim(first), dimnames = dimnames(first))
    })
  }
  # Keyed SS tables union semantic rows and numeric columns before reduction.
  if (is.data.frame(first)) {
    return(pst_distribution_reduce_data_frames(
      nodes, probs, path, backend = backend, block_bytes = block_bytes
    ))
  }
  # Lists recurse over the semantic union of field names.
  if (is.list(first) && !is.data.frame(first)) {
    fields <- sort(unique(unlist(lapply(nodes, names), use.names = FALSE)))
    return(setNames(lapply(fields, function(field) {
      child <- lapply(nodes, function(node) node[[field]])
      pst_distribution_reduce_ss_node(
        child, probs, paste0(path, "$", field), backend, block_bytes
      )
    }), fields))
  }
  # Named numeric vectors union labels and fill absent categories with zero.
  if (is.numeric(first) && is.null(dim(first)) && !is.null(names(first))) {
    labels <- sort(unique(unlist(lapply(nodes, names), use.names = FALSE)))
    aligned <- lapply(nodes, function(node) {
      values <- setNames(numeric(length(labels)), labels)
      if (!is.null(node)) values[names(node)] <- as.numeric(node)
      values
    })
    reduced <- pst_distribution_reduce_arrays(
      lapply(aligned, function(x) array(x, dim = length(x), dimnames = list(names(x)))),
      probs, backend, block_bytes, path
    )
    reduced <- lapply(reduced, function(x) setNames(as.numeric(x), labels))
    reduced$n_observed <- setNames(as.integer(reduced$n_observed), labels)
    return(reduced)
  }
  # One-dimensional tables are semantic named vectors despite retaining a dim.
  if (is.numeric(first) && length(dim(first)) == 1L && !is.null(dimnames(first)[[1L]])) {
    labels <- sort(unique(unlist(lapply(nodes, function(node) dimnames(node)[[1L]]),
                                 use.names = FALSE)))
    aligned <- lapply(nodes, function(node) {
      values <- setNames(numeric(length(labels)), labels)
      node_labels <- dimnames(node)[[1L]]
      values[node_labels] <- as.numeric(node)
      array(values, dim = length(labels), dimnames = list(labels))
    })
    reduced <- pst_distribution_reduce_arrays(aligned, probs, backend, block_bytes, path)
    reduced <- lapply(reduced, function(x) setNames(as.numeric(x), labels))
    reduced$n_observed <- setNames(as.integer(reduced$n_observed), labels)
    return(reduced)
  }
  # Labelled SS matrices union semantic rows and columns independently.
  if (is.matrix(first) && is.numeric(first) &&
      all(vapply(nodes, function(node) !is.null(rownames(node)) && !is.null(colnames(node)), logical(1)))) {
    row_labels <- sort(unique(unlist(lapply(nodes, rownames), use.names = FALSE)))
    column_labels <- sort(unique(unlist(lapply(nodes, colnames), use.names = FALSE)))
    aligned <- lapply(nodes, function(node) {
      values <- matrix(0, nrow = length(row_labels), ncol = length(column_labels),
                       dimnames = list(row_labels, column_labels))
      values[rownames(node), colnames(node)] <- node
      values
    })
    return(pst_distribution_reduce_arrays(aligned, probs, backend, block_bytes, path))
  }
  # Numeric scalars and equal-shaped arrays reduce directly.
  if (is.numeric(first)) {
    if (any(vapply(nodes, function(node) !identical(dim(node), dim(first)), logical(1)))) {
      stop(path, " numeric dimensions do not match", call. = FALSE)
    }
    if (is.null(dim(first))) return(pst_distribution_reduce_numeric(unlist(nodes), probs))
    return(pst_distribution_reduce_arrays(nodes, probs, backend, block_bytes, path))
  }
  # Non-numeric structure must be invariant across posterior samples.
  if (!all(vapply(nodes, identical, logical(1), first))) {
    stop(path, " structural values do not match", call. = FALSE)
  }
  first
}

#' Build exact SS distribution surfaces.
#'
#' @param samples Trajectory objects containing SS payloads.
#' @param probs Quantile probabilities.
#' @return SS hierarchy with numeric leaves replaced by four summaries.
#' @noRd
pst_distribution_build_ss <- function(samples, probs = c(0.025, 0.975),
                                      backend = "R",
                                      block_bytes = 16 * 1024^2) {
  pst_distribution_reduce_ss_node(
    lapply(samples, `[[`, "ss"), probs, backend = backend,
    block_bytes = block_bytes
  )
}

#' Expand one TT matrix onto a canonical step-function grid.
#'
#' @param matrix_value Numeric TT matrix.
#' @param source_time Source time coordinates.
#' @param canonical_time Union time coordinates.
#' @param columns Canonical semantic columns.
#' @param tolerance Time comparison tolerance.
#' @return Canonically aligned numeric matrix with absent columns zero-filled.
#' @noRd
pst_distribution_align_tt_matrix <- function(matrix_value, source_time,
                                             canonical_time, columns, tolerance) {
  order_index <- order(source_time)
  source_time <- source_time[order_index]
  matrix_value <- matrix_value[order_index, , drop = FALSE]
  out <- matrix(0, nrow = length(canonical_time), ncol = length(columns),
                dimnames = list(NULL, columns))
  # Carry the last committed row forward at every canonical time.
  for (time_index in seq_along(canonical_time)) {
    candidates <- which(source_time <= canonical_time[[time_index]] + tolerance)
    if (!length(candidates)) {
      stop("TT canonical time precedes a sample range", call. = FALSE)
    }
    out[time_index, colnames(matrix_value)] <- matrix_value[max(candidates), ]
  }
  # Embedded public time columns are coordinates, not posterior measurements.
  time_columns <- columns[grepl("(^time$|^times$|time_vec$)", columns)]
  if (length(time_columns)) {
    out[, time_columns] <- canonical_time
  }
  out
}

#' Recursively align and reduce TT matrices.
#'
#' @param nodes Parallel TT nodes.
#' @param source_times Per-sample time vectors.
#' @param canonical_time Canonical union time vector.
#' @param probs Quantile probabilities.
#' @param tolerance Time tolerance.
#' @param path Diagnostic path.
#' @return TT hierarchy with matrix leaves replaced by four summary matrices.
#' @noRd
pst_distribution_reduce_tt_node <- function(nodes, source_times, canonical_time,
                                            probs, tolerance, path = "TT",
                                            backend = "R",
                                            block_bytes = 16 * 1024^2,
                                            alignment_context = NULL) {
  first <- nodes[[1L]]
  # Nested TT families recurse over the union of public field names.
  if (is.list(first) && !is.data.frame(first)) {
    fields <- sort(unique(unlist(lapply(nodes, names), use.names = FALSE)))
    fields <- setdiff(fields, c("metadata", "time"))
    return(setNames(lapply(fields, function(field) {
      pst_distribution_reduce_tt_node(lapply(nodes, `[[`, field), source_times,
                                      canonical_time, probs, tolerance,
                                      paste0(path, "$", field), backend,
                                      block_bytes, alignment_context)
    }), fields))
  }
  # TT matrices use row position as time and semantic column labels as identities.
  if (is.matrix(first) && is.numeric(first)) {
    if (any(vapply(seq_along(nodes), function(index) nrow(nodes[[index]]) != length(source_times[[index]]), logical(1)))) {
      stop(path, " rows are not aligned to TT time", call. = FALSE)
    }
    columns <- sort(unique(unlist(lapply(nodes, colnames), use.names = FALSE)))
    # The C++ backend fuses step alignment and reduction without aligned R copies.
    if (identical(backend, "cpp")) {
      return(pst_distribution_align_reduce_tt_cpp(
        nodes, source_times, canonical_time, columns, probs, tolerance,
        block_bytes, path, alignment_context
      ))
    }
    aligned <- lapply(seq_along(nodes), function(index) {
      pst_distribution_align_tt_matrix(nodes[[index]], source_times[[index]],
                                       canonical_time, columns, tolerance)
    })
    return(pst_distribution_reduce_arrays(aligned, probs, backend, block_bytes, path))
  }
  # Invariant non-matrix TT structure is preserved only when exactly equal.
  if (!all(vapply(nodes, identical, logical(1), first))) {
    stop(path, " unsupported or mismatched TT structure", call. = FALSE)
  }
  first
}

#' Canonicalize and preserve TT lookup metadata across samples.
#'
#' @param samples Trajectory samples.
#' @param tolerance Time alignment tolerance.
#' @return Reducer metadata with canonical state/path lookups and source records.
#' @noRd
pst_distribution_build_tt_metadata <- function(samples, tolerance) {
  source_metadata <- lapply(samples, function(sample) sample$tt$metadata)
  states <- sort(unique(unlist(lapply(source_metadata, function(metadata) {
    if (is.null(metadata$state_lookup$state)) character() else as.character(metadata$state_lookup$state)
  }), use.names = FALSE)))
  paths <- sort(unique(unlist(lapply(source_metadata, function(metadata) {
    if (is.null(metadata$path_lookup$label)) character() else as.character(metadata$path_lookup$label)
  }), use.names = FALSE)))
  state_lookup <- data.frame(state_id = seq_along(states), state = states,
                             stringsAsFactors = FALSE)
  parent_labels <- vapply(paths, function(label) {
    parts <- strsplit(label, "|", fixed = TRUE)[[1L]]
    parts <- parts[nzchar(parts)]
    if (length(parts) <= 1L) return(NA_character_)
    paste0("|", paste(parts[-length(parts)], collapse = "|"))
  }, character(1))
  added_states <- vapply(paths, function(label) {
    parts <- strsplit(label, "|", fixed = TRUE)[[1L]]
    parts <- parts[nzchar(parts)]
    if (!length(parts)) NA_character_ else parts[[length(parts)]]
  }, character(1))
  path_lookup <- data.frame(
    path_id = seq_along(paths), label = paths,
    parent_path_id = match(parent_labels, paths),
    added_state_id = match(added_states, states),
    stringsAsFactors = FALSE
  )
  list(
    n_samples = length(samples), time_tolerance = tolerance,
    alignment = "right_continuous_step",
    state_lookup = state_lookup, path_lookup = path_lookup,
    source_metadata = source_metadata
  )
}

#' Build step-aligned exact TT distribution surfaces.
#'
#' @param samples Trajectory objects containing TT payloads.
#' @param probs Quantile probabilities.
#' @param tolerance Time tolerance for endpoints and simultaneous events.
#' @return TT hierarchy with canonical time and four-summary numeric matrices.
#' @noRd
pst_distribution_build_tt <- function(samples, probs = c(0.025, 0.975),
                                      tolerance = 1e-8, backend = "R",
                                      block_bytes = 16 * 1024^2) {
  source_times <- lapply(samples, function(sample) as.numeric(sample$tt$time))
  ranges <- lapply(source_times, range)
  reference_range <- ranges[[1L]]
  # Shared-phylogeny mode rejects endpoint differences instead of filling them.
  for (sample_index in seq_along(ranges)) {
    if (any(abs(ranges[[sample_index]] - reference_range) > tolerance)) {
      stop("sample ", sample_index, " TT range does not match", call. = FALSE)
    }
  }
  canonical_time <- pst_distribution_canonical_times(source_times, tolerance)
  alignment_context <- NULL
  # Compute C++ step-row maps once for the complete TT object, not once per family.
  if (identical(backend, "cpp")) {
    order_indices <- lapply(source_times, order)
    row_maps <- lapply(seq_along(source_times), function(sample_index) {
      sorted_time <- source_times[[sample_index]][order_indices[[sample_index]]]
      interval <- findInterval(canonical_time + tolerance, sorted_time)
      if (any(interval < 1L)) {
        stop("TT canonical time precedes sample ", sample_index, call. = FALSE)
      }
      order_indices[[sample_index]][interval] - 1L
    })
    alignment_context <- list(order_indices = order_indices, row_maps = row_maps)
  }
  nodes <- lapply(samples, function(sample) sample$tt[setdiff(names(sample$tt), c("time", "metadata"))])
  out <- c(list(time = canonical_time),
           pst_distribution_reduce_tt_node(nodes, source_times, canonical_time,
                                           probs, tolerance, backend = backend,
                                           block_bytes = block_bytes,
                                           alignment_context = alignment_context))
  out$metadata <- pst_distribution_build_tt_metadata(samples, tolerance)
  class(out) <- c("pst_distribution_tt", "list")
  out
}
