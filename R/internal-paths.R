# Internal path and root helpers used by the V35 input boundary.

pst_observed_states <- function(simmap) {
  unique(unlist(lapply(simmap$maps, names), use.names = FALSE))
}

pst_root_edge_indices <- function(simmap) {
  root <- length(simmap$tip.label) + 1L
  which(simmap$edge[, 1L] == root)
}
