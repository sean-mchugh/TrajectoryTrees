# V35 input adapter.
#
# Role:
#   Convert a user-facing `simmap` phylogeny into the compact integer/double
#   payload that the V35 C++ kernel consumes.
#
# Shape note:
#   The returned payload is an internal list with one row per phylogeny edge and
#   one row per map segment. Integer ids are 1-based at the R boundary.
#   Debug artifact example, once captured:
#   payload <- readRDS("testing/correctness/reference/v35_intermediates/five_tip_cpp_example/input_payload.rds")

#' Validate that a user input tree has the simmap fields required by V35.
#'
#' Invisibly returns the input tree; throws on missing or misaligned fields.
#' @param phylo_charmap Input `phylo`/`simmap` object whose `maps` and `mapped.edge` fields describe mapped character histories.
#' @return Invisibly returns the input tree; throws on missing or misaligned fields.
#' @noRd
pst_v35_validate_simmap <- function(phylo_charmap) {
  # Validate the minimum simmap fields before any C++ bridge work can occur.
  required <- c("edge", "edge.length", "tip.label", "maps", "mapped.edge")
  missing <- setdiff(required, names(phylo_charmap))
  # Reject inputs before traversal when required simmap fields are absent.
  if (length(missing)) {
    stop(
      "V35 requires a mapped simmap with fields: ",
      paste(required, collapse = ", "),
      "; missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  # Validate edge/map cardinality so map offsets can index every edge safely.
  nedge <- nrow(phylo_charmap$edge)
  # Reject edge-parallel fields that would break C++ offset indexing.
  if (length(phylo_charmap$maps) != nedge) {
    stop("V35 simmap maps must be parallel to the edge matrix", call. = FALSE)
  }
  # Reject edge-parallel fields that would break C++ offset indexing.
  if (length(phylo_charmap$edge.length) != nedge) {
    stop("V35 edge.length must be parallel to the edge matrix", call. = FALSE)
  }

  invisible(phylo_charmap)
}

#' Choose the root-anchor state policy used to seed V35 path histories.
#'
#' A root-policy list describing anchor state, provenance, daughter states, and synthetic-root status.
#' @param phylo_charmap Input `phylo`/`simmap` object whose `maps` and `mapped.edge` fields describe mapped character histories.
#' @param root_state Optional explicit root-anchor state; when `NULL`, V35 infers or synthesizes one from root-descending edges.
#' @return A root-policy list describing anchor state, provenance, daughter states, and synthetic-root status.
#' @noRd
pst_v35_resolve_root_policy <- function(phylo_charmap, root_state = NULL) {
  # Keep root identity consistent with V17, but avoid scanning every map unless
  # a synthetic root name must be collision-checked.
  pst_v35_validate_simmap(phylo_charmap)
  root_edges <- pst_root_edge_indices(phylo_charmap)
  # Reject trees whose root cannot seed a root-descending traversal.
  if (!length(root_edges)) {
    stop("The simmap has no root-descending edges", call. = FALSE)
  }
  daughter_states <- vapply(
    root_edges,
    function(i) names(phylo_charmap$maps[[i]])[[1L]],
    character(1)
  )

  # Use the caller-supplied root anchor instead of inferred daughter-state policy.
  if (!is.null(root_state)) {
    # Validate explicit root anchors before they enter the state lookup.
    if (length(root_state) != 1L || is.na(root_state) || !nzchar(root_state)) {
      stop("root_state must be one nonempty, non-missing state", call. = FALSE)
    }
    return(list(
      anchor = as.character(root_state),
      provenance = "explicit",
      daughter_states = daughter_states,
      synthetic = FALSE
    ))
  }

  shared <- unique(daughter_states)
  # Infer a real root anchor when all root-descending edges share the same starting state.
  if (length(shared) == 1L) {
    return(list(
      anchor = shared,
      provenance = "shared-root-edge-inference",
      daughter_states = daughter_states,
      synthetic = FALSE
    ))
  }

  observed <- colnames(phylo_charmap$mapped.edge)
  # Fall back to scanning maps when `mapped.edge` lacks state column names.
  if (is.null(observed) || !length(observed)) {
    observed <- pst_observed_states(phylo_charmap)
  }
  anchor <- "ROOT"
  suffix <- 0L
  # Advance the synthetic root label until it cannot collide with an observed mapped state.
  while (anchor %in% observed) {
    suffix <- suffix + 1L
    anchor <- paste0("ROOT#", suffix)
  }
  list(
    anchor = anchor,
    provenance = "synthetic-root",
    daughter_states = daughter_states,
    synthetic = TRUE
  )
}

#' Build the integer state lookup shared by R metadata and C++ traversal.
#'
#' A data frame with `state_id` and `state` columns.
#' @param phylo_charmap Input `phylo`/`simmap` object whose `maps` and `mapped.edge` fields describe mapped character histories.
#' @param root_policy Resolved root-anchor policy list.
#' @return A data frame with `state_id` and `state` columns.
#' @noRd
pst_v35_state_lookup <- function(phylo_charmap, root_policy) {
  # Gather all observed map states and the root anchor into a stable lookup.
  map_states <- colnames(phylo_charmap$mapped.edge)
  # Fall back to map names only when `mapped.edge` does not declare its state domain.
  if (is.null(map_states) || !length(map_states)) {
    map_states <- unique(unlist(lapply(phylo_charmap$maps, names), use.names = FALSE))
  }
  state_levels <- sort(unique(c(as.character(root_policy$anchor), map_states)))

  # Return both table and named id vector for fast R-side encoding.
  data.frame(
    state_id = seq_along(state_levels),
    state = state_levels,
    stringsAsFactors = FALSE
  )
}

#' Encode a simmap into the full compact vector payload consumed by the C++ smoke and traversal kernels.
#'
#' A list of edge, map, state, root-policy, label, and option vectors.
#' @param phylo_charmap Input `phylo`/`simmap` object whose `maps` and `mapped.edge` fields describe mapped character histories.
#' @param root_state Optional explicit root-anchor state; when `NULL`, V35 infers or synthesizes one from root-descending edges.
#' @param tt_mats Logical flag requesting through-time matrix construction.
#' @param ss Logical flag requesting empirical summary-statistics construction.
#' @param debug Logical flag that attaches low-level C++ round-trip diagnostics to the result.
#' @param debug_snapshots Logical flag recording whether future V35 debug snapshots are requested.
#' @param debug_snapshot_scope Snapshot scope; V35 currently accepts only `"batch"`.
#' @param debug_snapshot_steps Optional traversal-step selectors retained for later debug snapshots.
#' @param debug_snapshot_times Optional time selectors retained for later debug snapshots.
#' @param time_tolerance Positive numeric resolution used for every V35 time
#'   comparison, including event batches and scenario-tree cleanup.
#' @return A list of edge, map, state, root-policy, label, and option vectors.
#' @noRd
pst_v35_encode_input <- function(
    phylo_charmap,
    root_state = NULL,
    tt_mats = TRUE,
    ss = TRUE,
    debug = FALSE,
    debug_snapshots = FALSE,
    debug_snapshot_scope = "batch",
    debug_snapshot_steps = NULL,
    debug_snapshot_times = NULL,
    time_tolerance = 1e-6) {
  # Validate the simmap before constructing any integer id payload.
  pst_v35_validate_simmap(phylo_charmap)

  # Resolve root-state policy in R so the C++ kernel receives only ids.
  root_policy <- pst_v35_resolve_root_policy(phylo_charmap, root_state = root_state)
  state_lookup <- pst_v35_state_lookup(phylo_charmap, root_policy)
  state_ids <- setNames(state_lookup$state_id, state_lookup$state)

  # Encode map segments as flat vectors with an nedge + 1 offset table.
  map_lengths <- lengths(phylo_charmap$maps)
  map_edge_id <- rep.int(seq_along(phylo_charmap$maps), map_lengths)
  map_duration_named <- unlist(phylo_charmap$maps, use.names = TRUE)
  map_state_id <- unname(state_ids[names(map_duration_named)])
  map_duration <- unname(map_duration_named)
  map_edge_offset <- as.integer(c(1L, cumsum(map_lengths) + 1L))

  # Package the R/C++ boundary object in one list with documented id spaces.
  list(
    ntips = as.integer(length(phylo_charmap$tip.label)),
    nedge = as.integer(nrow(phylo_charmap$edge)),
    edge_parent = as.integer(phylo_charmap$edge[, 1L]),
    edge_child = as.integer(phylo_charmap$edge[, 2L]),
    edge_length = as.numeric(phylo_charmap$edge.length),
    map_edge_id = as.integer(map_edge_id),
    map_state_id = as.integer(map_state_id),
    map_duration = as.numeric(map_duration),
    map_edge_offset = as.integer(map_edge_offset),
    tip_label = as.character(phylo_charmap$tip.label),
    state_lookup = state_lookup,
    root_anchor_state_ids = as.integer(unname(state_ids[as.character(root_policy$anchor)])),
    root_policy = root_policy,
    options = list(
      tt_mats = isTRUE(tt_mats),
      ss = isTRUE(ss),
      debug = isTRUE(debug),
      debug_snapshots = list(
        enabled = isTRUE(debug_snapshots),
        scope = debug_snapshot_scope,
        steps = debug_snapshot_steps,
        times = debug_snapshot_times
      ),
      time_tolerance = as.numeric(time_tolerance)
    )
  )
}

#' Build the lightweight metadata payload used when C++ parses the simmap tree directly.
#'
#' A list containing counts, labels, lookup tables, root anchors, observed states, and options.
#' @param phylo_charmap Input `phylo`/`simmap` object whose `maps` and `mapped.edge` fields describe mapped character histories.
#' @param root_state Optional explicit root-anchor state; when `NULL`, V35 infers or synthesizes one from root-descending edges.
#' @param tt_mats Logical flag requesting through-time matrix construction.
#' @param ss Logical flag requesting empirical summary-statistics construction.
#' @param debug Logical flag that attaches low-level C++ round-trip diagnostics to the result.
#' @param debug_snapshots Logical flag recording whether future V35 debug snapshots are requested.
#' @param debug_snapshot_scope Snapshot scope; V35 currently accepts only `"batch"`.
#' @param debug_snapshot_steps Optional traversal-step selectors retained for later debug snapshots.
#' @param debug_snapshot_times Optional time selectors retained for later debug snapshots.
#' @param time_tolerance Positive numeric resolution used for every V35 time
#'   comparison, including event batches and scenario-tree cleanup.
#' @return A list containing counts, labels, lookup tables, root anchors, observed states, and options.
#' @noRd
pst_v35_encode_input_metadata <- function(
    phylo_charmap,
    root_state = NULL,
    tt_mats = TRUE,
    ss = TRUE,
    debug = FALSE,
    debug_snapshots = FALSE,
    debug_snapshot_scope = "batch",
    debug_snapshot_steps = NULL,
    debug_snapshot_times = NULL,
    time_tolerance = 1e-6) {
  # Lightweight metadata path for C++ direct simmap parsing.
  pst_v35_validate_simmap(phylo_charmap)
  root_policy <- pst_v35_resolve_root_policy(phylo_charmap, root_state = root_state)
  state_lookup <- pst_v35_state_lookup(phylo_charmap, root_policy)
  state_ids <- setNames(state_lookup$state_id, state_lookup$state)
  observed_state_ids <- as.integer(unname(state_ids[colnames(phylo_charmap$mapped.edge)]))
  observed_state_ids <- observed_state_ids[!is.na(observed_state_ids)]

  list(
    ntips = as.integer(length(phylo_charmap$tip.label)),
    nedge = as.integer(nrow(phylo_charmap$edge)),
    tip_label = as.character(phylo_charmap$tip.label),
    state_lookup = state_lookup,
    observed_state_ids = sort(unique(observed_state_ids)),
    root_anchor_state_ids = as.integer(unname(state_ids[as.character(root_policy$anchor)])),
    root_policy = root_policy,
    options = list(
      tt_mats = isTRUE(tt_mats),
      ss = isTRUE(ss),
      debug = isTRUE(debug),
      debug_snapshots = list(
        enabled = isTRUE(debug_snapshots),
        scope = debug_snapshot_scope,
        steps = debug_snapshot_steps,
        times = debug_snapshot_times
      ),
      time_tolerance = as.numeric(time_tolerance)
    )
  )
}
