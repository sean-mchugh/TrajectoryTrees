#ifndef PST_V34_TT_OUTPUT_HPP
#define PST_V34_TT_OUTPUT_HPP

#include "pst_v34_types.hpp"

/** Build cumulative event, arrival, departure, and path views on the shared batch axis. */
Rcpp::List pst_v34_transition_views_cpp(
  Rcpp::DataFrame records,
  Rcpp::NumericMatrix canonical_trans,
  Rcpp::List summary,
  Rcpp::NumericVector time_vec,
  double time_tolerance);

/** Derive public total and state views from one canonical path matrix. */
Rcpp::List pst_v34_path_counter_views_cpp(
  Rcpp::NumericMatrix path_matrix,
  Rcpp::DataFrame path_lookup,
  Rcpp::DataFrame state_lookup);
#endif
