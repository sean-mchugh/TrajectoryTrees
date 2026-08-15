# Public orchestration for PST posterior distribution objects.

#' Count numeric cells whose observed posterior denominator is reduced.
#'
#' @param node Distribution surface node.
#' @param n_samples Full posterior sample count.
#' @return Number of cells with `n_observed < n_samples`.
#' @noRd
pst_distribution_count_under_observed <- function(node, n_samples) {
  if (!is.list(node)) return(0L)
  count <- 0L
  # Recurse through every distribution node and count explicit observation arrays.
  for (field in names(node)) {
    value <- node[[field]]
    if (identical(field, "n_observed") && is.numeric(value)) {
      count <- count + sum(value < n_samples, na.rm = TRUE)
    } else if (is.list(value)) {
      count <- count + pst_distribution_count_under_observed(value, n_samples)
    }
  }
  as.integer(count)
}

#' Build one posterior distribution object from trajectory samples.
#'
#' Validates a same-phylogeny ensemble, builds an interval-level MAP phylogeny,
#' unions transition-tree edges, and reduces SS/TT numeric surfaces exactly.
#'
#' @param samples Nonempty in-memory list of completed trajectory objects.
#' @param probs Ordered lower and upper empirical quantile probabilities.
#' @param tolerance Nonnegative tolerance for mapped and TT time coordinates.
#' @param include Any of `P`, `T`, `SS`, and `TT`.
#' @param backend Exact reduction backend, `"R"` or `"cpp"`.
#' @param block_bytes Positive C++ working-buffer target in bytes.
#' @return An object of class `pst_distribution`.
#' @noRd
get_pst_distribution <- function(samples, probs = c(0.025, 0.975),
                                 tolerance = 1e-8,
                                 include = c("P", "T", "SS", "TT"),
                                 backend = "cpp",
                                 block_bytes = 16 * 1024^2) {
  total_started <- proc.time()[["elapsed"]]
  include <- unique(toupper(include))
  # Validate global numeric controls before any surface-specific work.
  if (length(probs) != 2L || any(!is.finite(probs)) ||
      probs[[1L]] < 0 || probs[[2L]] > 1 || probs[[1L]] > probs[[2L]]) {
    stop("probs must contain ordered lower and upper probabilities", call. = FALSE)
  }
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be one nonnegative finite value", call. = FALSE)
  }
  # The initial release deliberately exposes only its verified exact R backend.
  if (length(backend) != 1L || !backend %in% c("R", "cpp")) {
    stop("backend must be \"R\" or \"cpp\"", call. = FALSE)
  }
  if (length(block_bytes) != 1L || !is.finite(block_bytes) || block_bytes <= 0) {
    stop("block_bytes must be one positive finite value", call. = FALSE)
  }
  if (!length(include) || any(!include %in% c("P", "T", "SS", "TT"))) {
    stop("include may contain only P, T, SS, and TT", call. = FALSE)
  }
  validation <- pst_distribution_validate_samples(
    samples, tolerance, require_tt = "TT" %in% include
  )
  # Required surfaces are checked before any expensive canonicalization begins.
  for (surface in include) {
    field <- switch(surface, P = "phylo", T = "trans", SS = "ss", TT = "tt")
    missing_samples <- which(vapply(samples, function(sample) is.null(sample[[field]]), logical(1)))
    if (length(missing_samples)) {
      stop("sample ", missing_samples[[1L]], " is missing ", surface, call. = FALSE)
    }
  }
  topology_id <- paste(c(as.vector(samples[[1L]]$phylo$edge),
                         samples[[1L]]$phylo$tip.label), collapse = ":")
  # Build requested surfaces independently so timing and memory diagnostics
  # expose the stage that a later C++ backend would need to accelerate.
  stage_seconds <- setNames(rep(0, 4L), c("P", "T", "SS", "TT"))
  stage_bytes <- setNames(rep(0, 4L), c("P", "T", "SS", "TT"))
  built <- setNames(vector("list", 4L), c("P", "T", "SS", "TT"))
  builders <- list(
    P = function() pst_distribution_build_p(samples, tolerance),
    T = function() pst_distribution_build_t(samples, probs, backend, block_bytes),
    SS = function() pst_distribution_build_ss(samples, probs, backend, block_bytes),
    TT = function() pst_distribution_build_tt(
      samples, probs, tolerance, backend, block_bytes
    )
  )
  # Materialize each requested surface once and capture its elapsed/object size.
  for (surface in intersect(c("P", "T", "SS", "TT"), include)) {
    stage_started <- proc.time()[["elapsed"]]
    built[[surface]] <- builders[[surface]]()
    stage_seconds[[surface]] <- proc.time()[["elapsed"]] - stage_started
    stage_bytes[[surface]] <- as.numeric(object.size(built[[surface]]))
  }
  diagnostics <- list(
    p_intervals = if (is.null(built$P)) 0L else nrow(attr(built$P, "posterior_intervals")),
    t_union_edges = if (is.null(built$T)) 0L else nrow(built$T$edges),
    t_absent_edge_observations = if (is.null(built$T)) 0L else
      sum(validation$n_samples - built$T$edges$present_count),
    t_explicit_zero_observations = if (is.null(built$T)) 0L else
      sum(built$T$edges$explicit_zero_count),
    ss_under_observed_cells = if (is.null(built$SS)) 0L else
      pst_distribution_count_under_observed(built$SS, validation$n_samples),
    tt_under_observed_cells = if (is.null(built$TT)) 0L else
      pst_distribution_count_under_observed(built$TT, validation$n_samples)
  )
  stage_seconds <- c(stage_seconds,
                     total = proc.time()[["elapsed"]] - total_started)
  # Retain a stable five-field public object shape regardless of include flags.
  out <- list(
    metadata = list(
      n_samples = validation$n_samples, probs = as.numeric(probs),
      topology_id = topology_id, time_alignment = "shared_phylogeny",
      time_tolerance = tolerance, map_tie_policy = "lexical_state_label",
      missing_numeric_category = "zero", missing_value_policy = "exclude_NA_report_n_observed",
      backend = backend, block_bytes = block_bytes, stage_seconds = stage_seconds,
      stage_object_bytes = stage_bytes, diagnostics = diagnostics
    ),
    P = built$P, T = built$T, SS = built$SS, TT = built$TT
  )
  class(out) <- c("pst_distribution", "list")
  out
}

#' Build a trajectory-object distribution with an exact R or C++ backend.
#'
#' Primary public entry point for posterior P/T/SS/TT aggregation. The historical
#' `get_pst_distribution()` name remains available with identical behavior.
#' @param samples Nonempty in-memory list of completed trajectory objects.
#' @param probs Ordered lower and upper empirical quantile probabilities.
#' @param tolerance Nonnegative time and branch-length tolerance.
#' @param include Any of `P`, `T`, `SS`, and `TT`.
#' @param backend Exact reduction backend, `"cpp"` or `"R"`.
#' @param block_bytes Positive C++ working-buffer target in bytes.
#' @return An object of class `pst_distribution`.
#' @noRd
get_trjectory_obj_dist <- function(
    samples,
    probs = c(0.025, 0.975),
    tolerance = 1e-8,
    include = c("P", "T", "SS", "TT"),
    backend = "cpp",
    block_bytes = 16 * 1024^2) {
  get_pst_distribution(
    samples = samples, probs = probs, tolerance = tolerance,
    include = include, backend = backend, block_bytes = block_bytes
  )
}
