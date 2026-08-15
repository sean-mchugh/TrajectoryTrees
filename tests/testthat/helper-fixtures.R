small_simmap <- function() {
  readRDS(system.file(
    "extdata", "generated_000001.rds",
    package = "TrajectoryTrees",
    mustWork = TRUE
  ))
}

public_function <- function(name) {
  tryCatch(
    getExportedValue("TrajectoryTrees", name),
    error = function(error) error
  )
}
