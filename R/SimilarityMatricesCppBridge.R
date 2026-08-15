# Package-native bridge for the ScenarioMats similarity calculation.

#' Verify that the similarity native routine is loaded.
#'
#' @return Invisibly returns `TRUE`; throws if the package DLL is unavailable.
#' @noRd
scenario_mats_v3_load_cpp <- function() {
  if (!is.loaded("_TrajectoryTrees_scenario_mats_v3_calculate_cpp")) {
    stop(
      "TrajectoryTrees similarity native routine is not loaded; reinstall the package",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
