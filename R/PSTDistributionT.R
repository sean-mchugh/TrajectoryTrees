# Union transition-tree aggregation.

#' List biological terminal labels descended from one transition-tree node.
#'
#' @param tree Transition tree.
#' @param node Node identifier.
#' @return Sorted terminal labels below `node`.
#' @noRd
pst_distribution_t_descendant_tips <- function(tree, node) {
  ntip <- length(tree$tip.label)
  # A terminal node owns its public represented-tip identity directly.
  if (node <= ntip) return(tree$tip.label[[node]])
  frontier <- node
  tips <- character()
  # Traverse descendants breadth-first so sibling edges with equal path labels
  # remain distinguishable by their disjoint represented-tip sets.
  while (length(frontier)) {
    children <- tree$edge[tree$edge[, 1L] %in% frontier, 2L]
    tips <- c(tips, tree$tip.label[children[children <= ntip]])
    frontier <- children[children > ntip]
  }
  sort(unique(tips))
}

#' Extract rooted semantic edge records from one transition tree.
#'
#' @param tree Transition tree with edge, maps, edge lengths, and tip labels.
#' @param sample_index One-based sample index retained for reduction.
#' @return Edge-record data frame whose identities do not depend on positive length.
#' @noRd
pst_distribution_t_edge_records <- function(tree, sample_index) {
  ntip <- length(tree$tip.label)
  children <- tree$edge[, 2L]
  parents <- tree$edge[, 1L]
  root <- setdiff(unique(parents), children)
  # A transition tree must have one rooted topology for semantic path identities.
  if (length(root) != 1L) {
    stop("sample ", sample_index, " transition tree has no unique root", call. = FALSE)
  }
  node_paths <- setNames("ROOT", as.character(root))
  pending <- seq_len(nrow(tree$edge))
  records <- vector("list", nrow(tree$edge))
  segment_labels <- vapply(tree$maps, function(map) paste(names(map), collapse = "|"), character(1))
  roles <- ifelse(children <= ntip, "terminal", "internal")
  sibling_bases <- paste(parents, segment_labels, roles, sep = "\037")
  sibling_counts <- table(sibling_bases)
  ambiguous_internal <- roles == "internal" & sibling_counts[sibling_bases] > 1L
  # Public V31 T trees expose no persistent lineage/scenario id that could align
  # same-parent, same-path internal siblings across posterior samples. Reject
  # this ambiguity rather than overwrite or invent sample-dependent identities.
  if (any(ambiguous_internal)) {
    stop("sample ", sample_index,
         " has ambiguous duplicate internal transition edges", call. = FALSE)
  }
  # Resolve edges only after their parent node path exists, accumulating rooted keys.
  while (length(pending)) {
    progressed <- FALSE
    next_pending <- integer()
    # Visit unresolved edges and serialize those whose parent is already known.
    for (edge_index in pending) {
      parent_path <- node_paths[[as.character(parents[[edge_index]])]]
      # Defer descendants until an earlier pass has serialized their parent.
      if (is.null(parent_path)) {
        next_pending <- c(next_pending, edge_index)
        next
      }
      child <- children[[edge_index]]
      terminal <- child <= ntip
      terminal_identity <- if (terminal) tree$tip.label[[child]] else ""
      # Internal descendant sets are not stable cross-sample identities and are
      # unnecessary once ambiguous same-label siblings are rejected above.
      descendant_identity <- if (terminal) terminal_identity else ""
      segment_label <- segment_labels[[edge_index]]
      role <- if (terminal) "terminal" else "internal"
      edge_key <- paste(parent_path, segment_label, role,
                        if (terminal) terminal_identity else "",
                        sep = "\037")
      records[[edge_index]] <- data.frame(
        sample_index = sample_index, edge_index = edge_index,
        parent_key = parent_path, edge_key = edge_key,
        segment_label = segment_label, edge_role = role,
        terminal_identity = terminal_identity,
        descendant_identity = descendant_identity,
        size = as.numeric(tree$edge.length[[edge_index]]),
        stringsAsFactors = FALSE
      )
      node_paths[[as.character(child)]] <- edge_key
      progressed <- TRUE
    }
    # Lack of progress proves a cycle or disconnected topology.
    if (!progressed && length(next_pending)) {
      stop("sample ", sample_index, " transition tree cannot be rooted", call. = FALSE)
    }
    pending <- next_pending
  }
  do.call(rbind, records)
}

#' Serialize transition-tree structure for within-call record caching.
#'
#' @param tree Transition tree.
#' @return Exact character key excluding sample-specific edge sizes.
#' @noRd
pst_distribution_t_structure_key <- function(tree) {
  paste(
    paste(as.integer(tree$edge), collapse = ","),
    paste(tree$tip.label, collapse = "\036"),
    paste(vapply(tree$maps, function(map) paste(names(map), collapse = "|"),
                 character(1)), collapse = "\036"),
    sep = "\035"
  )
}

#' Reduce transition-sequence lineage and scenario support by rooted path.
#'
#' This is the distribution-owned equivalent of the manuscript transition-tree
#' summary. It retains the quantities needed by the plotter so plotting never
#' needs the original trajectory sample.
#'
#' @param samples Trajectory objects containing V32 transition path-size maps.
#' @param probs Lower and upper quantile probabilities.
#' @return Path labels, state labels, and internal/terminal support summaries.
#' @noRd
pst_distribution_t_path_summary <- function(samples, probs) {
  sample_paths <- lapply(samples, function(sample) {
    maps <- sample$trans$path_size_maps$lineage_size_path_maps
    unique(names(unlist(maps, use.names = TRUE)))
  })
  paths <- sort(setdiff(unique(unlist(sample_paths, use.names = FALSE)), "ROOT"))
  columns <- c("lineages", "scenarios", "postprob")
  empty_matrix <- function() {
    matrix(
      0,
      nrow = length(paths),
      ncol = length(columns),
      dimnames = list(paths, columns)
    )
  }
  sample_matrices <- lapply(seq_along(samples), function(sample_index) {
    tree <- samples[[sample_index]]$trans
    lineage_maps <- tree$path_size_maps$lineage_size_path_maps
    scenario_maps <- tree$path_size_maps$scenario_size_path_maps
    if (is.null(lineage_maps) || is.null(scenario_maps) ||
        length(lineage_maps) != length(tree$maps) ||
        length(scenario_maps) != length(tree$maps) ||
        any(lengths(lineage_maps) != lengths(tree$maps)) ||
        any(lengths(scenario_maps) != lengths(tree$maps))) {
      stop(
        "sample ", sample_index,
        " transition tree has malformed path size maps",
        call. = FALSE
      )
    }
    internal <- empty_matrix()
    terminal <- empty_matrix()
    all_steps <- unlist(tree$maps, use.names = TRUE)
    for (edge_index in seq_along(tree$maps)) {
      is_parent_edge <- tree$edge[edge_index, 2L] %in% tree$edge[, 1L]
      for (segment_index in seq_along(tree$maps[[edge_index]])) {
        step <- tree$maps[[edge_index]][[segment_index]]
        path <- names(lineage_maps[[edge_index]])[[segment_index]]
        if (identical(path, "ROOT")) next
        lineage <- as.numeric(lineage_maps[[edge_index]][[segment_index]])
        scenario <- as.numeric(scenario_maps[[edge_index]][[segment_index]])
        values <- c(lineages = lineage, scenarios = scenario, postprob = 1)
        if (step > 0) {
          internal[path, ] <- values
          zero_step_present <- any(
            all_steps == 0 & names(all_steps) == path
          )
          if (!zero_step_present && !is_parent_edge) {
            terminal[path, ] <- values
          }
        } else {
          terminal[path, ] <- values
        }
      }
    }
    list(internal_edges = internal, tip_edges = terminal)
  })
  reduce_kind <- function(kind) {
    result <- lapply(
      c("lower", "mean", "median", "upper"),
      function(field) empty_matrix()
    )
    names(result) <- c("lower", "mean", "median", "upper")
    for (path_index in seq_along(paths)) {
      for (column_index in seq_along(columns)) {
        values <- vapply(
          sample_matrices,
          function(sample) sample[[kind]][path_index, column_index],
          numeric(1)
        )
        # Manuscript size summaries condition on the transition sequence being
        # present; posterior probabilities retain the explicit zero replicates.
        if (columns[[column_index]] != "postprob") {
          values <- values[values != 0]
        }
        reduced <- if (length(values)) {
          pst_distribution_reduce_numeric(values, probs)
        } else {
          list(lower = 0, mean = 0, median = 0, upper = 0)
        }
        for (field in names(result)) {
          result[[field]][path_index, column_index] <- reduced[[field]]
        }
      }
    }
    result
  }
  internal <- reduce_kind("internal_edges")
  terminal <- reduce_kind("tip_edges")
  statistics <- lapply(names(internal), function(field) {
    list(
      internal_edges = internal[[field]],
      tip_edges = terminal[[field]]
    )
  })
  names(statistics) <- names(internal)
  list(
    paths = paths,
    states = sort(unique(colnames(samples[[1L]]$phylo$mapped.edge))),
    statistics = statistics,
    size_conditioning = "present_replicates",
    interval = probs
  )
}

#' Build a union transition tree with unconditional edge-size summaries.
#'
#' @param samples Trajectory objects containing transition trees.
#' @param probs Lower and upper quantile probabilities.
#' @param backend Exact numeric backend, `"R"` or `"cpp"`.
#' @param block_bytes Positive C++ samples-by-edge block target in bytes.
#' @return A `pst_distribution_transition_tree` with union topology and edge summaries.
#' @noRd
pst_distribution_build_t <- function(samples, probs = c(0.025, 0.975),
                                     backend = "R",
                                     block_bytes = 16 * 1024^2) {
  transition_summary <- pst_distribution_t_path_summary(samples, probs)
  record_cache <- new.env(parent = emptyenv(), hash = TRUE)
  sample_records <- vector("list", length(samples))
  # Reuse structural records for repeated T topologies while applying each sample's sizes.
  for (sample_index in seq_along(samples)) {
    tree <- samples[[sample_index]]$trans
    structure_key <- pst_distribution_t_structure_key(tree)
    if (exists(structure_key, envir = record_cache, inherits = FALSE)) {
      records <- get(structure_key, envir = record_cache, inherits = FALSE)
    } else {
      records <- pst_distribution_t_edge_records(tree, 0L)
      assign(structure_key, records, envir = record_cache)
    }
    records$sample_index <- sample_index
    records$size <- as.numeric(tree$edge.length[records$edge_index])
    sample_records[[sample_index]] <- records
  }
  records <- do.call(rbind, sample_records)
  keys <- unique(records$edge_key)
  # Parent keys contain their complete ancestry, so separator count gives depth.
  depth <- nchar(keys) - nchar(gsub("\037", "", keys, fixed = TRUE))
  keys <- keys[order(depth, keys)]
  cpp_summaries <- NULL
  # Build and release bounded edge-size blocks for the C++ backend.
  if (identical(backend, "cpp")) {
    pst_distribution_load_cpp()
    fields <- c("lower", "mean", "median", "upper")
    cpp_summaries <- setNames(lapply(fields, function(field) numeric(length(keys))), fields)
    key_indices <- match(records$edge_key, keys)
    elements_per_block <- pst_distribution_cpp_elements_per_block(
      block_bytes, length(samples)
    )
    # Fill one samples-by-union-edge block, reduce it, and release it before advancing.
    for (block_start in seq.int(1L, length(keys), by = elements_per_block)) {
      block_end <- min(length(keys), block_start + elements_per_block - 1L)
      aligned_sizes <- matrix(
        0, nrow = length(samples), ncol = block_end - block_start + 1L
      )
      selected <- which(key_indices >= block_start & key_indices <= block_end)
      aligned_sizes[cbind(
        records$sample_index[selected], key_indices[selected] - block_start + 1L
      )] <- records$size[selected]
      reduced <- pst_distribution_cpp_reduce_exact(aligned_sizes, probs)
      for (field in fields) {
        cpp_summaries[[field]][block_start:block_end] <- reduced[[field]]
      }
    }
  }
  edge_rows <- vector("list", length(keys))
  # Reduce every union edge over all samples, filling structural absence with zero.
  for (key_index in seq_along(keys)) {
    key <- keys[[key_index]]
    matched <- records[records$edge_key == key, , drop = FALSE]
    summaries <- if (is.null(cpp_summaries)) {
      values <- numeric(length(samples))
      values[matched$sample_index] <- matched$size
      pst_distribution_reduce_numeric(values, probs)
    } else {
      lapply(c("lower", "mean", "median", "upper"), function(field) {
        cpp_summaries[[field]][[key_index]]
      }) |> setNames(c("lower", "mean", "median", "upper"))
    }
    edge_rows[[key_index]] <- data.frame(
      edge_key = key, parent_key = matched$parent_key[[1L]],
      segment_label = matched$segment_label[[1L]],
      edge_role = matched$edge_role[[1L]],
      terminal_identity = matched$terminal_identity[[1L]],
      descendant_identity = matched$descendant_identity[[1L]],
      lower = summaries$lower, mean = summaries$mean,
      median = summaries$median, upper = summaries$upper,
      present_count = nrow(matched),
      presence_probability = nrow(matched) / length(samples),
      explicit_zero_count = sum(matched$size == 0),
      stringsAsFactors = FALSE
    )
  }
  edges <- do.call(rbind, edge_rows)
  terminal_rows <- which(edges$edge_role == "terminal")
  internal_rows <- which(edges$edge_role == "internal")
  ntip <- length(terminal_rows)
  child_ids <- integer(nrow(edges))
  child_ids[terminal_rows] <- seq_len(ntip)
  root_id <- ntip + 1L
  child_ids[internal_rows] <- ntip + 1L + seq_along(internal_rows)
  internal_node_by_key <- setNames(child_ids[internal_rows], edges$edge_key[internal_rows])
  parent_ids <- vapply(edges$parent_key, function(parent_key) {
    # Root children attach to the unique root; descendants attach to their
    # parent internal edge's child node.
    if (identical(parent_key, "ROOT")) return(root_id)
    parent_id <- internal_node_by_key[[parent_key]]
    if (is.null(parent_id)) stop("union T parent edge is not internal", call. = FALSE)
    parent_id
  }, integer(1))
  map_list <- lapply(seq_len(nrow(edges)), function(edge_index) {
    setNames(edges$mean[[edge_index]], edges$segment_label[[edge_index]])
  })
  state_labels <- sort(unique(edges$segment_label))
  mapped_edge <- matrix(0, nrow = nrow(edges), ncol = length(state_labels),
                        dimnames = list(NULL, state_labels))
  # Populate the mean mapped-edge view while retaining all four size vectors.
  for (edge_index in seq_len(nrow(edges))) {
    mapped_edge[edge_index, edges$segment_label[[edge_index]]] <- edges$mean[[edge_index]]
  }
  out <- list(
    edge = cbind(parent = parent_ids, child = child_ids),
    edge.length = edges$mean,
    tip.label = make.unique(paste(edges$terminal_identity[terminal_rows],
                                  edges$segment_label[terminal_rows], sep = "@")),
    Nnode = length(internal_rows) + 1L,
    maps = map_list,
    mapped.edge = mapped_edge,
    edge.length_distribution = list(lower = edges$lower, mean = edges$mean,
                                    median = edges$median, upper = edges$upper),
    transition_summary = transition_summary,
    edges = edges
  )
  class(out) <- c("pst_distribution_transition_tree", "phylo")
  out
}
