# V34 native summary-stat bridge.
#
# V34 keeps PST, TT, and SS construction C++ owned. This R file is deliberately
# thin: it loads the compiled C++ symbol and attaches the returned payload to the
# canonical trajectory object. It must not rebuild, reshape, or patch summary
# statistics in R.

#' Build V34 canonical empirical summary statistics for public PST outputs.
#'
#' A summary-statistics list, usually delegated to the C++ implementation.
#' @param tt Optional native or legacy TT object used by summary-statistics builders.
#' @param tree_parts C++ metrics captured as P/S/T trees were materialized.
#' @param path_size_maps Traversal-backed transition path support.
#' @return A summary-statistics list, usually delegated to the C++ implementation.
#' @noRd
pst_v34_build_summary_stats_canonical <- function(
    tt,
    tree_parts,
    path_size_maps) {
  # Missing parts indicate an ownership defect. V34 has no tree-scanning or
  # legacy-summary fallback that could conceal it.
  if (is.null(tree_parts)) {
    stop(
      "V34 SS assembly requires traversal/materialization-owned summary parts",
      call. = FALSE
    )
  }
  pst_v34_cpp_summary_stats_from_parts_payload(
    tree_parts = tree_parts,
    path_size_maps = path_size_maps,
    tt = tt
  )
}
