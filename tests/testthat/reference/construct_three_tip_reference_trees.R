# Hand-constructed three-tip mapped character trees for ScenarioMats inspection.
#
# Every tree has the same rooted topology:
#   root -> tip_1            length 5
#   root -> internal node    length 2
#   internal -> tip_2        length 3
#   internal -> tip_3        length 3
#
# Non-ultrametric cases shorten exactly one terminal edge. All mapped changes
# occur at integer times, so histories can be read directly against a 0:5 grid.

make_three_tip_reference_simmap <- function(
    root_to_tip_1_map,
    root_to_internal_map,
    internal_to_tip_2_map,
    internal_to_tip_3_map,
    case_id,
    root_state = "A") {
  edge_maps <- list(
    root_to_tip_1_map,
    root_to_internal_map,
    internal_to_tip_2_map,
    internal_to_tip_3_map
  )

  for (edge_index in seq_along(edge_maps)) {
    edge_map <- edge_maps[[edge_index]]
    if (!is.numeric(edge_map) || !length(edge_map) ||
        is.null(names(edge_map)) || any(!nzchar(names(edge_map))) ||
        any(!is.finite(edge_map)) || any(edge_map <= 0)) {
      stop(
        "Every edge map must be a named, positive numeric vector: ",
        case_id,
        call. = FALSE
      )
    }
  }

  edge_matrix <- matrix(
    c(
      4L, 1L,
      4L, 5L,
      5L, 2L,
      5L, 3L
    ),
    ncol = 2L,
    byrow = TRUE
  )
  edge_lengths <- vapply(edge_maps, sum, numeric(1))
  observed_states <- sort(unique(unlist(lapply(edge_maps, names))))
  mapped_edge <- matrix(
    0,
    nrow = length(edge_maps),
    ncol = length(observed_states),
    dimnames = list(NULL, observed_states)
  )

  for (edge_index in seq_along(edge_maps)) {
    for (state in unique(names(edge_maps[[edge_index]]))) {
      mapped_edge[edge_index, state] <- sum(
        edge_maps[[edge_index]][names(edge_maps[[edge_index]]) == state]
      )
    }
  }

  terminal_edge_indices <- c(1L, 3L, 4L)
  tip_states <- vapply(
    edge_maps[terminal_edge_indices],
    function(edge_map) tail(names(edge_map), 1L),
    character(1)
  )
  names(tip_states) <- c("tip_1", "tip_2", "tip_3")

  structure(
    list(
      edge = edge_matrix,
      edge.length = edge_lengths,
      tip.label = c("tip_1", "tip_2", "tip_3"),
      Nnode = 2L,
      maps = edge_maps,
      mapped.edge = mapped_edge,
      states = tip_states,
      case_id = case_id,
      reference_root_state = root_state
    ),
    class = c("simmap", "phylo")
  )
}

.three_tip_next_states <- function(parent_state) {
  state_order <- c("A", "B", "C")
  parent_position <- match(parent_state, state_order)
  if (is.na(parent_position)) {
    stop("Unknown parent state: ", parent_state, call. = FALSE)
  }
  c(
    state_order[parent_position %% length(state_order) + 1L],
    state_order[(parent_position + 1L) %% length(state_order) + 1L]
  )
}

.three_tip_daughter_states <- function(parent_state, cladogenesis_type) {
  new_states <- .three_tip_next_states(parent_state)
  switch(
    cladogenesis_type,
    both_daughters_inherit = c(parent_state, parent_state),
    one_daughter_inherits = c(parent_state, new_states[[1L]]),
    neither_daughter_inherits_same_state = c(new_states[[1L]], new_states[[1L]]),
    neither_daughter_inherits_different_states = c(new_states[[1L]], new_states[[2L]]),
    stop("Unknown cladogenesis type: ", cladogenesis_type, call. = FALSE)
  )
}

.three_tip_trim_terminal_map <- function(edge_map, shortened_length) {
  remaining_length <- shortened_length
  trimmed_durations <- numeric()
  trimmed_states <- character()

  for (segment_index in seq_along(edge_map)) {
    if (remaining_length <= 0) break
    retained_duration <- min(edge_map[[segment_index]], remaining_length)
    trimmed_durations <- c(trimmed_durations, retained_duration)
    trimmed_states <- c(trimmed_states, names(edge_map)[[segment_index]])
    remaining_length <- remaining_length - retained_duration
  }

  names(trimmed_durations) <- trimmed_states
  trimmed_durations
}

.three_tip_path_label <- function(first_edge_map, second_edge_map = NULL) {
  states <- names(first_edge_map)
  if (!is.null(second_edge_map)) states <- c(states, names(second_edge_map))
  states <- states[c(TRUE, states[-1L] != states[-length(states)])]
  paste0("|", paste(states, collapse = "|"))
}

.three_tip_case_record <- function(
    case_id,
    description,
    group,
    root_cladogenesis,
    internal_cladogenesis,
    root_to_tip_1_map,
    root_to_internal_map,
    internal_to_tip_2_map,
    internal_to_tip_3_map,
    root_state = "A",
    shortened_tip = NA_character_,
    shortening = 0) {
  list(
    case_id = case_id,
    description = description,
    group = group,
    root_state = root_state,
    root_cladogenesis = root_cladogenesis,
    internal_cladogenesis = internal_cladogenesis,
    shortened_tip = shortened_tip,
    shortening = shortening,
    maps = list(
      root_to_tip_1_map,
      root_to_internal_map,
      internal_to_tip_2_map,
      internal_to_tip_3_map
    ),
    expected_tip_paths = c(
      tip_1 = .three_tip_path_label(root_to_tip_1_map),
      tip_2 = .three_tip_path_label(
        root_to_internal_map,
        internal_to_tip_2_map
      ),
      tip_3 = .three_tip_path_label(
        root_to_internal_map,
        internal_to_tip_3_map
      )
    )
  )
}

.three_tip_reference_case_records <- function() {
  cladogenesis_types <- c(
    "both_daughters_inherit",
    "one_daughter_inherits",
    "neither_daughter_inherits_same_state",
    "neither_daughter_inherits_different_states"
  )
  short_cladogenesis_names <- c(
    both_daughters_inherit = "both_inherit",
    one_daughter_inherits = "one_inherits",
    neither_daughter_inherits_same_state = "neither_same",
    neither_daughter_inherits_different_states = "neither_different"
  )

  records <- list()
  record_index <- 0L

  # The first 16 cases cross all root and internal cladogenetic categories.
  # Their edges contain no anagenetic change, which keeps each cladogenetic
  # relationship visually isolated.
  for (root_cladogenesis in cladogenesis_types) {
    root_daughter_states <- .three_tip_daughter_states(
      "A",
      root_cladogenesis
    )
    internal_parent_state <- root_daughter_states[[2L]]

    for (internal_cladogenesis in cladogenesis_types) {
      internal_daughter_states <- .three_tip_daughter_states(
        internal_parent_state,
        internal_cladogenesis
      )
      record_index <- record_index + 1L
      case_id <- paste0(
        "core__root_", short_cladogenesis_names[[root_cladogenesis]],
        "__internal_", short_cladogenesis_names[[internal_cladogenesis]]
      )
      records[[record_index]] <- .three_tip_case_record(
        case_id = case_id,
        description = paste(
          "Core cladogenesis case:",
          root_cladogenesis,
          "at the root and",
          internal_cladogenesis,
          "at time 2."
        ),
        group = "cladogenesis_core",
        root_cladogenesis = root_cladogenesis,
        internal_cladogenesis = internal_cladogenesis,
        root_to_tip_1_map = setNames(5, root_daughter_states[[1L]]),
        root_to_internal_map = setNames(2, root_daughter_states[[2L]]),
        internal_to_tip_2_map = setNames(3, internal_daughter_states[[1L]]),
        internal_to_tip_3_map = setNames(3, internal_daughter_states[[2L]])
      )
    }
  }

  # These nine cases add anagenetic paths and repeated-state overlaps. Every
  # duration and transition time remains an integer.
  long_path_records <- list(
    .three_tip_case_record(
      "paths__unit_cycles_both_both",
      "Unit cycles with equal inheritance at both cladogenetic events.",
      "long_paths",
      "both_daughters_inherit",
      "both_daughters_inherit",
      c(A = 1, B = 1, C = 1, A = 1, B = 1),
      c(A = 1, B = 1),
      c(B = 1, C = 1, A = 1),
      c(B = 1, C = 1, A = 1)
    ),
    .three_tip_case_record(
      "paths__shifted_cycles_neither_different_both",
      "Shifted cyclic paths reach the same states by offset histories.",
      "long_paths",
      "neither_daughter_inherits_different_states",
      "both_daughters_inherit",
      c(B = 1, C = 1, A = 1, B = 1, C = 1),
      c(C = 1, A = 1),
      c(A = 1, B = 1, C = 1),
      c(A = 1, B = 1, C = 1)
    ),
    .three_tip_case_record(
      "paths__two_unit_state_one_one",
      "A two-unit ancestral state followed by partial inheritance twice.",
      "long_paths",
      "one_daughter_inherits",
      "one_daughter_inherits",
      c(A = 2, B = 1, C = 1, A = 1),
      c(B = 2),
      c(B = 2, C = 1),
      c(C = 1, A = 2)
    ),
    .three_tip_case_record(
      "paths__two_unit_middle_both_neither_same",
      "A two-unit middle state with equal novel daughter states at time 2.",
      "long_paths",
      "both_daughters_inherit",
      "neither_daughter_inherits_same_state",
      c(A = 1, B = 2, C = 1, A = 1),
      c(A = 1, C = 1),
      c(B = 1, C = 1, A = 1),
      c(B = 1, C = 2)
    ),
    .three_tip_case_record(
      "paths__two_unit_return_both_neither_different",
      "Different novel daughters enter shifted return paths.",
      "long_paths",
      "both_daughters_inherit",
      "neither_daughter_inherits_different_states",
      c(A = 1, B = 1, C = 2, A = 1),
      c(A = 2),
      c(B = 1, C = 1, A = 1),
      c(C = 1, A = 1, B = 1)
    ),
    .three_tip_case_record(
      "paths__different_routes_neither_same_one",
      "Sisters reach A and B through different post-cladogenesis routes.",
      "long_paths",
      "neither_daughter_inherits_same_state",
      "one_daughter_inherits",
      c(B = 1, C = 1, A = 1, B = 1, C = 1),
      c(B = 1, C = 1),
      c(C = 1, A = 1, B = 1),
      c(B = 1, A = 1, B = 1)
    ),
    .three_tip_case_record(
      "paths__different_routes_one_neither_different",
      "Partial root inheritance precedes different novel sister states.",
      "long_paths",
      "one_daughter_inherits",
      "neither_daughter_inherits_different_states",
      c(A = 1, B = 1, C = 1, A = 2),
      c(B = 2),
      c(C = 1, A = 1, B = 1),
      c(A = 1, C = 1, B = 1)
    ),
    .three_tip_case_record(
      "paths__mixed_returns_neither_different_neither_same",
      "Different root daughters precede shared novel sister inheritance.",
      "long_paths",
      "neither_daughter_inherits_different_states",
      "neither_daughter_inherits_same_state",
      c(B = 1, C = 1, A = 1, B = 2),
      c(C = 1, A = 1),
      c(B = 1, C = 1, A = 1),
      c(B = 2, C = 1)
    ),
    .three_tip_case_record(
      "paths__conserved_homoplasy_same_route",
      paste0(
        "Tip 2 retains ancestral B while tip 3 leaves B for A and ",
        "reacquires B through the ancestral A-to-B route."
      ),
      "long_paths",
      "one_daughter_inherits",
      "one_daughter_inherits",
      c(C = 5),
      c(A = 1, B = 1),
      c(B = 3),
      c(A = 1, B = 2)
    )
  )
  records <- c(records, long_path_records)

  # Each non-ultrametric case is a transparent copy of one long-path case with
  # one terminal map truncated by exactly one or two time units.
  shortening_specs <- data.frame(
    source_position = seq_len(8L),
    shortened_tip = c(
      "tip_1", "tip_1", "tip_2", "tip_2",
      "tip_3", "tip_3", "tip_2", "tip_1"
    ),
    shortening = c(1, 2, 1, 2, 1, 2, 1, 1),
    stringsAsFactors = FALSE
  )

  for (shortening_index in seq_len(nrow(shortening_specs))) {
    source_record <- long_path_records[[
      shortening_specs$source_position[[shortening_index]]
    ]]
    shortened_tip <- shortening_specs$shortened_tip[[shortening_index]]
    shortening <- shortening_specs$shortening[[shortening_index]]
    shortened_maps <- source_record$maps
    shortened_edge_index <- c(tip_1 = 1L, tip_2 = 3L, tip_3 = 4L)[[
      shortened_tip
    ]]
    shortened_maps[[shortened_edge_index]] <- .three_tip_trim_terminal_map(
      shortened_maps[[shortened_edge_index]],
      sum(shortened_maps[[shortened_edge_index]]) - shortening
    )

    record_index <- length(records) + 1L
    records[[record_index]] <- .three_tip_case_record(
      case_id = paste0(
        "nonultra__", sub("^paths__", "", source_record$case_id),
        "__short_", shortened_tip, "_by_", shortening
      ),
      description = paste0(
        source_record$description,
        " The ", shortened_tip, " terminal edge is shortened by ",
        shortening, "."
      ),
      group = "non_ultrametric",
      root_cladogenesis = source_record$root_cladogenesis,
      internal_cladogenesis = source_record$internal_cladogenesis,
      root_to_tip_1_map = shortened_maps[[1L]],
      root_to_internal_map = shortened_maps[[2L]],
      internal_to_tip_2_map = shortened_maps[[3L]],
      internal_to_tip_3_map = shortened_maps[[4L]],
      root_state = source_record$root_state,
      shortened_tip = shortened_tip,
      shortening = shortening
    )
  }

  records
}

three_tip_reference_case_specs <- function() {
  records <- .three_tip_reference_case_records()
  rows <- vector("list", length(records))

  for (record_index in seq_along(records)) {
    record <- records[[record_index]]
    row <- data.frame(
      case_id = record$case_id,
      description = record$description,
      group = record$group,
      root_state = record$root_state,
      root_cladogenesis = record$root_cladogenesis,
      internal_cladogenesis = record$internal_cladogenesis,
      shortened_tip = record$shortened_tip,
      shortening = record$shortening,
      ultrametric = is.na(record$shortened_tip),
      stringsAsFactors = FALSE
    )
    row$root_to_tip_1_map <- I(list(record$maps[[1L]]))
    row$root_to_internal_map <- I(list(record$maps[[2L]]))
    row$internal_to_tip_2_map <- I(list(record$maps[[3L]]))
    row$internal_to_tip_3_map <- I(list(record$maps[[4L]]))
    row$expected_tip_paths <- I(list(record$expected_tip_paths))
    rows[[record_index]] <- row
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

build_three_tip_reference_simmaps <- function() {
  case_specs <- three_tip_reference_case_specs()
  trees <- vector("list", nrow(case_specs))
  names(trees) <- case_specs$case_id

  for (case_index in seq_len(nrow(case_specs))) {
    trees[[case_index]] <- make_three_tip_reference_simmap(
      root_to_tip_1_map = case_specs$root_to_tip_1_map[[case_index]],
      root_to_internal_map = case_specs$root_to_internal_map[[case_index]],
      internal_to_tip_2_map = case_specs$internal_to_tip_2_map[[case_index]],
      internal_to_tip_3_map = case_specs$internal_to_tip_3_map[[case_index]],
      case_id = case_specs$case_id[[case_index]],
      root_state = case_specs$root_state[[case_index]]
    )
  }

  trees
}

three_tip_reference_pairwise_state_overlaps <- function(simmap_tree) {
  map_intervals <- function(edge_map, edge_start_time) {
    interval_end <- edge_start_time + cumsum(edge_map)
    data.frame(
      state = names(edge_map),
      start = c(edge_start_time, head(interval_end, -1L)),
      end = interval_end,
      stringsAsFactors = FALSE
    )
  }

  tip_intervals <- list(
    tip_1 = map_intervals(simmap_tree$maps[[1L]], 0),
    tip_2 = rbind(
      map_intervals(simmap_tree$maps[[2L]], 0),
      map_intervals(simmap_tree$maps[[3L]], 2)
    ),
    tip_3 = rbind(
      map_intervals(simmap_tree$maps[[2L]], 0),
      map_intervals(simmap_tree$maps[[4L]], 2)
    )
  )
  overlap_rows <- list()
  overlap_row_index <- 0L

  for (first_tip_index in 1:2) {
    for (second_tip_index in (first_tip_index + 1L):3L) {
      first_tip_name <- names(tip_intervals)[[first_tip_index]]
      second_tip_name <- names(tip_intervals)[[second_tip_index]]
      first_history <- tip_intervals[[first_tip_index]]
      second_history <- tip_intervals[[second_tip_index]]
      shared_states <- intersect(first_history$state, second_history$state)

      for (state in shared_states) {
        first_state_intervals <- first_history[
          first_history$state == state,
          ,
          drop = FALSE
        ]
        second_state_intervals <- second_history[
          second_history$state == state,
          ,
          drop = FALSE
        ]
        overlap_duration <- 0

        for (first_interval_index in seq_len(nrow(first_state_intervals))) {
          for (second_interval_index in seq_len(nrow(second_state_intervals))) {
            overlap_duration <- overlap_duration + max(
              0,
              min(
                first_state_intervals$end[[first_interval_index]],
                second_state_intervals$end[[second_interval_index]]
              ) -
                max(
                  first_state_intervals$start[[first_interval_index]],
                  second_state_intervals$start[[second_interval_index]]
                )
            )
          }
        }

        if (overlap_duration > 0) {
          overlap_row_index <- overlap_row_index + 1L
          overlap_rows[[overlap_row_index]] <- data.frame(
            tip_i = first_tip_name,
            tip_j = second_tip_name,
            state = state,
            overlap_duration = overlap_duration,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (!length(overlap_rows)) {
    return(data.frame(
      tip_i = character(),
      tip_j = character(),
      state = character(),
      overlap_duration = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, overlap_rows)
}

