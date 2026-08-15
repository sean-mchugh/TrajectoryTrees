#' TrajectoryTrees: evolutionary trajectory trees from stochastic maps
#'
#' The package compiles its native implementation during installation and
#' registers every native entry point through Rcpp.
#'
#' @keywords internal
#' @useDynLib TrajectoryTrees, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
