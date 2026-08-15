# Posterior MAP aggregation for mapped phylogenies.

#' Build an interval-level posterior MAP phylogeny.
#'
#' @param samples Trajectory objects sharing one phylogeny topology.
#' @param tolerance Boundary tolerance used to collapse mapped-state times.
#' @return A simmap-compatible MAP tree with `posterior_intervals` metadata.
#' @noRd
pst_distribution_build_p <- function(samples, tolerance = 1e-8) {
  reference <- samples[[1L]]$phylo
  n_edges <- nrow(reference$edge)
  rows <- vector("list", n_edges)
  map_list <- vector("list", n_edges)
  # Aggregate each biological edge independently on the union of map boundaries.
  for (edge_index in seq_len(n_edges)) {
    sample_maps <- lapply(samples, function(sample) sample$phylo$maps[[edge_index]])
    boundaries <- lapply(sample_maps, function(map) c(0, cumsum(as.numeric(map))))
    grid <- pst_distribution_canonical_times(boundaries, tolerance)
    edge_rows <- vector("list", length(grid) - 1L)
    edge_map <- numeric()
    # Evaluate categorical occupancy at each canonical interval midpoint.
    for (interval_index in seq_len(length(grid) - 1L)) {
      start <- grid[[interval_index]]
      end <- grid[[interval_index + 1L]]
      midpoint <- (start + end) / 2
      states <- vapply(seq_along(sample_maps), function(sample_index) {
        map <- sample_maps[[sample_index]]
        ends <- cumsum(as.numeric(map))
        names(map)[which(midpoint <= ends + tolerance)[[1L]]]
      }, character(1))
      counts <- table(states)
      winners <- sort(names(counts)[counts == max(counts)])
      winner <- winners[[1L]]
      edge_map <- c(edge_map, setNames(end - start, winner))
      edge_rows[[interval_index]] <- data.frame(
        edge = edge_index, start = start, end = end, map_state = winner,
        map_support = unname(max(counts) / length(samples)),
        tied = length(winners) > 1L, stringsAsFactors = FALSE
      )
      edge_rows[[interval_index]]$state_probabilities <- list(as.numeric(counts) / length(samples))
      names(edge_rows[[interval_index]]$state_probabilities[[1L]]) <- names(counts)
      edge_rows[[interval_index]]$tied_states <- list(winners)
    }
    rows[[edge_index]] <- do.call(rbind, edge_rows)
    map_list[[edge_index]] <- edge_map
  }
  out <- reference
  out$maps <- map_list
  out$edge.length <- vapply(map_list, function(map) sum(as.numeric(map)), numeric(1))
  states <- sort(unique(unlist(lapply(map_list, names), use.names = FALSE)))
  out$mapped.edge <- matrix(0, nrow = n_edges, ncol = length(states),
                            dimnames = list(NULL, states))
  # Rebuild mapped-edge totals from the MAP intervals without altering topology.
  for (edge_index in seq_len(n_edges)) {
    totals <- tapply(as.numeric(map_list[[edge_index]]), names(map_list[[edge_index]]), sum)
    out$mapped.edge[edge_index, names(totals)] <- totals
  }
  attr(out, "posterior_intervals") <- do.call(rbind, rows)
  class(out) <- unique(c("pst_distribution_map", class(reference)))
  out
}
