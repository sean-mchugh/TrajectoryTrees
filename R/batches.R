#' Construct trajectory trees in resumable disk batches
#'
#' @inheritParams make_trajectory_objects
#' @param output_dir Directory in which to write `manifest.rds` and batch files.
#' @param batch_size Positive integer or `"auto"`. Automatic batches contain
#'   20 trees at 500 or more tips, 50 trees at 251--499 tips, and otherwise one
#'   complete input set.
#' @param resume Continue a compatible incomplete manifest.
#' @param progress Emit progress messages.
#' @param compress Compression setting passed to [saveRDS()].
#'
#' @return A `trajectory_batch_manifest`.
#' @export
make_trajectory_batches <- function(
    simmaps,
    output_dir,
    batch_size = "auto",
    resume = FALSE,
    progress = TRUE,
    compress = FALSE,
    root_state = NULL,
    include_tt = TRUE,
    include_ss = TRUE,
    include_scenario_mats = TRUE,
    include_path_maps = TRUE,
    time_tolerance = NULL,
    fields = NULL,
    ...) {
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("output_dir must be one nonempty path", call. = FALSE)
  }
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume)) {
    stop("resume must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("progress must be TRUE or FALSE", call. = FALSE)
  }
  normalized <- trajectory_normalize_simmaps(simmaps)
  batch_size <- trajectory_normalize_batch_size(batch_size)
  if (!dir.exists(output_dir) &&
      !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("could not create output_dir: ", output_dir, call. = FALSE)
  }
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  manifest_path <- file.path(output_dir, "manifest.rds")
  signatures <- lapply(normalized$sets, function(set) {
    vapply(set, digest::digest, character(1), algo = "xxhash64")
  })

  if (file.exists(manifest_path)) {
    if (!isTRUE(resume)) {
      stop("manifest already exists; use resume=TRUE to continue", call. = FALSE)
    }
    manifest <- readRDS(manifest_path)
    if (!inherits(manifest, "trajectory_batch_manifest") ||
        !identical(manifest$input_signatures, signatures)) {
      stop("existing manifest is incompatible with these simmaps", call. = FALSE)
    }
  } else {
    set_records <- lapply(seq_along(normalized$sets), function(set_index) {
      set <- normalized$sets[[set_index]]
      max_tips <- max(vapply(set, function(tree) length(tree$tip.label), integer(1)))
      size <- if (identical(batch_size, "auto")) {
        trajectory_auto_batch_size(max_tips)
      } else {
        batch_size
      }
      groups <- trajectory_batch_indices(length(set), size)
      files <- file.path(
        sprintf("set_%06d", set_index),
        sprintf("batch_%06d.rds", seq_along(groups))
      )
      list(
        n_objects = length(set),
        source_groups = groups,
        files = files,
        complete = rep(FALSE, length(groups))
      )
    })
    manifest <- list(
      format_version = 1L,
      status = "incomplete",
      output_dir = output_dir,
      shape = normalized$shape,
      outer_names = normalized$outer_names,
      inner_names = normalized$inner_names,
      input_signatures = signatures,
      configuration = list(
        batch_size = batch_size,
        root_state = root_state,
        include_tt = include_tt,
        include_ss = include_ss,
        include_scenario_mats = include_scenario_mats,
        include_path_maps = include_path_maps,
        time_tolerance = time_tolerance,
        fields = fields,
        compress = compress
      ),
      sets = set_records,
      files = unlist(lapply(set_records, `[[`, "files"), use.names = FALSE)
    )
    class(manifest) <- c("trajectory_batch_manifest", "list")
    saveRDS(manifest, manifest_path)
  }

  total <- sum(vapply(manifest$sets, function(set) length(set$source_groups), integer(1)))
  completed <- sum(vapply(manifest$sets, function(set) sum(set$complete), integer(1)))
  for (set_index in seq_along(manifest$sets)) {
    record <- manifest$sets[[set_index]]
    set_dir <- dirname(file.path(output_dir, record$files[[1L]]))
    dir.create(set_dir, recursive = TRUE, showWarnings = FALSE)
    for (batch_index in seq_along(record$source_groups)) {
      if (isTRUE(manifest$sets[[set_index]]$complete[[batch_index]])) {
        next
      }
      indices <- record$source_groups[[batch_index]]
      objects <- make_trajectory_objects(
        normalized$sets[[set_index]][indices],
        root_state = root_state,
        include_tt = include_tt,
        include_ss = include_ss,
        include_scenario_mats = include_scenario_mats,
        include_path_maps = include_path_maps,
        time_tolerance = time_tolerance,
        fields = fields,
        ...
      )
      batch_path <- file.path(output_dir, record$files[[batch_index]])
      if (file.exists(batch_path)) {
        stop("pending batch file already exists: ", batch_path, call. = FALSE)
      }
      saveRDS(
        list(source_indices = indices, objects = objects),
        batch_path,
        compress = compress
      )
      manifest$sets[[set_index]]$complete[[batch_index]] <- TRUE
      completed <- completed + 1L
      saveRDS(manifest, manifest_path)
      if (isTRUE(progress)) {
        message("Saved trajectory batch ", completed, "/", total, ": ", batch_path)
      }
    }
  }
  manifest$status <- "complete"
  saveRDS(manifest, manifest_path)
  manifest
}

#' Load trajectory trees from a batch manifest
#'
#' @param manifest A `trajectory_batch_manifest` or path to `manifest.rds`.
#'
#' @return One trajectory tree or a shape-preserving list of trajectory trees.
#' @export
load_trajectory_batches <- function(manifest) {
  if (is.character(manifest) && length(manifest) == 1L && !is.na(manifest)) {
    path <- normalizePath(manifest, mustWork = TRUE)
    manifest <- readRDS(path)
    manifest$output_dir <- dirname(path)
  }
  if (!inherits(manifest, "trajectory_batch_manifest") ||
      !identical(manifest$format_version, 1L)) {
    stop("manifest must be a supported trajectory_batch_manifest or path", call. = FALSE)
  }
  if (!identical(manifest$status, "complete")) {
    stop("trajectory batch manifest is incomplete", call. = FALSE)
  }
  sets <- lapply(manifest$sets, function(record) vector("list", record$n_objects))
  for (set_index in seq_along(manifest$sets)) {
    record <- manifest$sets[[set_index]]
    for (batch_index in seq_along(record$files)) {
      path <- file.path(manifest$output_dir, record$files[[batch_index]])
      payload <- readRDS(path)
      sets[[set_index]][payload$source_indices] <- payload$objects
    }
  }
  normalized <- list(
    shape = manifest$shape,
    outer_names = manifest$outer_names,
    inner_names = manifest$inner_names
  )
  trajectory_restore_shape(sets, normalized)
}

trajectory_normalize_batch_size <- function(batch_size) {
  if (identical(batch_size, "auto")) {
    return(batch_size)
  }
  if (!is.numeric(batch_size) || length(batch_size) != 1L ||
      !is.finite(batch_size) || batch_size < 1 ||
      batch_size != as.integer(batch_size)) {
    stop("batch_size must be \"auto\" or one positive integer", call. = FALSE)
  }
  as.integer(batch_size)
}

trajectory_auto_batch_size <- function(n_tips) {
  if (n_tips >= 500L) return(20L)
  if (n_tips >= 251L) return(50L)
  Inf
}

trajectory_batch_indices <- function(n_objects, size) {
  indices <- seq_len(n_objects)
  if (is.infinite(size) || size >= n_objects) {
    return(list(indices))
  }
  unname(split(indices, ceiling(indices / size)))
}

#' @export
print.trajectory_batch_manifest <- function(x, ...) {
  cat("<trajectory_batch_manifest>\n")
  cat("  status:", x$status, "\n")
  cat("  objects:", sum(vapply(x$sets, `[[`, integer(1), "n_objects")), "\n")
  cat("  batches:", length(x$files), "\n")
  cat("  directory:", x$output_dir, "\n")
  invisible(x)
}
