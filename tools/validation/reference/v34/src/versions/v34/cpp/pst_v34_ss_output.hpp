#ifndef PST_V34_SS_OUTPUT_HPP
#define PST_V34_SS_OUTPUT_HPP

#include "pst_v34_types.hpp"

#include <string>

/**
 * Capture one finalized PST tree's compatibility metrics during C++ tree
 * materialization so later SS assembly never has to reread that public tree.
 */
Rcpp::List pst_v34_capture_tree_summary_parts_cpp(
  Rcpp::List tree,
  const std::string& suffix,
  const Rcpp::CharacterVector& state_levels,
  bool include_edge_boundaries);

/**
 * Assemble SS from materialization-time tree parts, traversal-owned path-size
 * support, and online TT counters without scanning public PST trees.
 */
Rcpp::List pst_v34_build_summary_stats_from_parts_cpp(
  Rcpp::List tree_parts,
  SEXP path_size_maps,
  Rcpp::Nullable<Rcpp::List> tt = R_NilValue);

#endif
