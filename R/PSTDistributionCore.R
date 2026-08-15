# Exact reduction and shared-phylogeny validation for PST distributions.

#' Reduce numeric posterior samples to four exact summaries.
#'
#' @param values Finite numeric values for one aligned distribution element.
#' @param probs Length-two lower and upper empirical quantile probabilities.
#' @return A list containing `lower`, `mean`, `median`, and `upper` scalars.
#' @noRd
pst_distribution_reduce_numeric <- function(values, probs = c(0.025, 0.975)) {
  values <- as.numeric(values)
  # Infinite values indicate malformed numeric input; documented NA/NaN values
  # are excluded and reported separately by array/table reducers.
  if (!length(values) || any(is.infinite(values))) {
    stop("distribution values must be nonempty and cannot be infinite", call. = FALSE)
  }
  values <- values[!is.na(values)]
  # An entirely undefined cell remains explicitly undefined in all summaries.
  if (!length(values)) {
    return(list(lower = NA_real_, mean = NA_real_, median = NA_real_, upper = NA_real_))
  }
  # Require an ordered interior probability pair so field names retain their
  # declared lower/upper meaning.
  if (length(probs) != 2L || any(!is.finite(probs)) ||
      probs[[1L]] < 0 || probs[[2L]] > 1 || probs[[1L]] > probs[[2L]]) {
    stop("probs must contain ordered lower and upper probabilities", call. = FALSE)
  }
  list(
    lower = unname(stats::quantile(values, probs[[1L]], names = FALSE, type = 7)),
    mean = mean(values),
    median = stats::median(values),
    upper = unname(stats::quantile(values, probs[[2L]], names = FALSE, type = 7))
  )
}

#' Form a deterministic tolerance-collapsed union of time coordinates.
#'
#' @param time_vectors List of numeric time vectors.
#' @param tolerance Maximum absolute separation merged into one time coordinate.
#' @return Sorted numeric vector using the first sorted value in each cluster.
#' @noRd
pst_distribution_canonical_times <- function(time_vectors, tolerance = 1e-8) {
  values <- sort(unique(as.numeric(unlist(time_vectors, use.names = FALSE))))
  # A missing or non-finite time grid cannot define TT step functions.
  if (!length(values) || any(!is.finite(values))) {
    stop("time vectors must be nonempty and finite", call. = FALSE)
  }
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be one nonnegative finite value", call. = FALSE)
  }
  canonical <- values[[1L]]
  # Retain one representative per sorted tolerance cluster; compare against the
  # retained representative so clustering cannot drift through chained values.
  for (value in values[-1L]) {
    if (abs(value - canonical[[length(canonical)]]) > tolerance) {
      canonical <- c(canonical, value)
    }
  }
  canonical
}

#' Validate common topology and TT range for a trajectory ensemble.
#'
#' @param samples Nonempty list of completed trajectory objects.
#' @param tolerance Numeric tolerance for shared TT start and end coordinates.
#' @param require_tt Whether TT time ranges must be present and equal.
#' @return Validation metadata containing sample count and common TT range.
#' @noRd
pst_distribution_validate_samples <- function(samples, tolerance = 1e-8,
                                              require_tt = TRUE) {
  # Reject empty input before indexing the reference sample.
  if (!is.list(samples) || !length(samples)) {
    stop("samples must be a nonempty list", call. = FALSE)
  }
  reference <- samples[[1L]]
  # The common phylogeny and TT range are mandatory for this first reducer mode.
  if (is.null(reference$phylo$edge) || is.null(reference$phylo$tip.label) ||
      is.null(reference$phylo$edge.length) || is.null(reference$phylo$maps)) {
    stop("sample 1 is missing phylo topology", call. = FALSE)
  }
  reference_lengths <- as.numeric(reference$phylo$edge.length)
  reference_range <- NULL
  # TT range validation is required only when the requested distribution includes TT.
  if (isTRUE(require_tt)) {
    if (is.null(reference$tt$time) || !length(reference$tt$time)) {
      stop("sample 1 is missing TT time", call. = FALSE)
    }
    reference_range <- range(as.numeric(reference$tt$time))
  }
  # Validate every sample directly against the first structural reference.
  for (sample_index in seq_along(samples)) {
    sample <- samples[[sample_index]]
    if (!identical(sample$phylo$edge, reference$phylo$edge) ||
        !identical(sample$phylo$tip.label, reference$phylo$tip.label)) {
      stop("sample ", sample_index, " phylo topology does not match", call. = FALSE)
    }
    sample_lengths <- as.numeric(sample$phylo$edge.length)
    if (length(sample_lengths) != length(reference_lengths) ||
        any(!is.finite(sample_lengths)) ||
        any(abs(sample_lengths - reference_lengths) > tolerance)) {
      stop("sample ", sample_index, " phylo branch lengths do not match", call. = FALSE)
    }
    if (length(sample$phylo$maps) != length(sample_lengths)) {
      stop("sample ", sample_index, " phylo maps are not edge-aligned", call. = FALSE)
    }
    map_lengths <- vapply(sample$phylo$maps, function(map) {
      sum(as.numeric(map))
    }, numeric(1))
    if (any(!is.finite(map_lengths)) ||
        any(abs(map_lengths - sample_lengths) > tolerance)) {
      stop("sample ", sample_index,
           " phylo map totals do not match branch lengths", call. = FALSE)
    }
    if (isTRUE(require_tt)) {
      if (is.null(sample$tt$time) || !length(sample$tt$time)) {
        stop("sample ", sample_index, " is missing TT time", call. = FALSE)
      }
      sample_range <- range(as.numeric(sample$tt$time))
      if (any(!is.finite(sample_range)) ||
          any(abs(sample_range - reference_range) > tolerance)) {
        stop("sample ", sample_index, " TT range does not match", call. = FALSE)
      }
    }
  }
  list(n_samples = length(samples), tt_range = unname(reference_range))
}
