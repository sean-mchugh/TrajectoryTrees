# V35 TT serialization.
#
# The C++ traversal owns every biological count and score. This file only
# packages committed matrices and builds the approved dependent transition
# views (`tt[[family]]$trans`) on the shared online time axis.

#' Build the native V35 TT object from canonical summary matrices and optional scenario tree context.
#'
#' A `pst_v35_tt` list containing LTT, STT, TT, transition-event families, branch lengths, and metadata.
#' @param summary Canonical V35 traversal summary with lookup tables, scenario matrices, and edge records.
#' @param time_tolerance Shared V35 resolution for time-axis comparisons.
#' @param scenario_tree Optional scenario tree used to align TT or size-support identities to public tree ids.
#' @return A `pst_v35_tt` list containing LTT, STT, TT, transition-event families, branch lengths, and metadata.
#' @noRd
pst_v35_build_tt_canonical <- function(summary, time_tolerance) {
  # Consume the traversal-owned counter payload. Its rows were captured from
  # live state at each committed event time and were not reconstructed by
  # scanning `scenario_summary` after traversal.
  tt_core <- summary$online_tt_core
  # Reject incomplete internal summaries instead of silently falling back to
  # the obsolete dense-matrix reconstruction path.
  if (is.null(tt_core)) {
    stop("V35 traversal did not return online TT counters", call. = FALSE)
  }
  # Each row is already one validated same-time transaction; serialization
  # must not collapse or reinterpret traversal events.
  time_vec <- as.numeric(tt_core$time_vec)
  # Path matrices are the sole numeric source for occupancy and exposure
  # families. C++ aggregates their state and total views from terminal path
  # identities, avoiding parallel traversal matrices with duplicate values.
  ltt <- pst_v35_cpp_tt_path_counter_views_payload(
    tt_core$lin_path,
    summary$path_lookup,
    summary$state_lookup
  )
  cltt <- pst_v35_cpp_tt_path_counter_views_payload(
    tt_core$cumulative_lin_path,
    summary$path_lookup,
    summary$state_lookup
  )
  stt <- pst_v35_cpp_tt_path_counter_views_payload(
    tt_core$scn_path,
    summary$path_lookup,
    summary$state_lookup
  )
  cstt <- pst_v35_cpp_tt_path_counter_views_payload(
    tt_core$cumulative_scn_path,
    summary$path_lookup,
    summary$state_lookup
  )
  sequence_tt <- pst_v35_cpp_tt_path_counter_views_payload(
    tt_core$seq_path,
    summary$path_lookup,
    summary$state_lookup
  )
  tlen <- pst_v35_cpp_tt_path_counter_views_payload(
    tt_core$tlen_path,
    summary$path_lookup,
    summary$state_lookup
  )
  change_records <- tt_core$change_records
  # Native LDIF records are created only after atomic batch commit resolves a
  # distinct public Leave. Root/Stay edges, whole-scenario steps, and provisional
  # same-time splits never enter this record family; duplicate biological
  # representations are already deduplicated.
  lineage_diff_records <- tt_core$lineage_differentiating_records

  path_records <- tt_core$path_records
  # Transition-event view contract:
  # `ctt`, `ldif`, and `puniq` are separate event families. Their live canonical
  # path matrices are the numerical source; records retain only legacy event
  # vocabularies for deterministic dependent formatting.
  ctt <- pst_v35_cpp_tt_transition_views_payload(
    change_records,
    tt_core$ctt_path,
    summary,
    time_vec,
    time_tolerance
  )
  ldif <- pst_v35_cpp_tt_transition_views_payload(
    lineage_diff_records,
    tt_core$ldif_path,
    summary,
    time_vec,
    time_tolerance
  )
  puniq <- pst_v35_cpp_tt_transition_views_payload(
    path_records,
    tt_core$puniq_path,
    summary,
    time_vec,
    time_tolerance
  )

  # Construct the V35 public TT shape directly from native C++ matrices;
  # `event_family` is intentionally not a V35 field.
  out <- list(
    time = time_vec,
    ltt = ltt,
    cltt = cltt,
    stt = stt,
    cstt = cstt,
    tt = sequence_tt,
    ctt = ctt,
    ldif = ldif,
    puniq = puniq,
    tlen = tlen,
    metadata = list(
      state_lookup = summary$state_lookup,
      # `terminal_state_id` is internal dependency metadata used to derive
      # state views from `trans`. Preserve the established four-column public
      # lookup contract in TT metadata.
      path_lookup = summary$path_lookup[
        , c("path_id", "label", "parent_path_id", "added_state_id"),
        drop = FALSE
      ],
      change_records = change_records,
      lineage_differentiating_records = lineage_diff_records,
      path_unique_records = path_records,
      path_through_counts = tt_core$path_through_counts[
        , c("path_id", "path", "lin_through", "scn_through"), drop = FALSE
      ],
      tt_event_tracker = attr(summary, "pst_v35_tt_event_tracker", exact = TRUE),
      construction_mode = tt_core$construction_mode,
      batch_semantics = tt_core$batch_semantics,
      scoring_phase = tt_core$scoring_phase,
      finalization_mode = tt_core$finalization_mode
    )
  )
  class(out) <- c("pst_v35_tt", class(out))
  out
}
