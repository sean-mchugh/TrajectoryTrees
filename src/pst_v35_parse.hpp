#ifndef PST_V35_PARSE_HPP
#define PST_V35_PARSE_HPP

#include "pst_v35_types.hpp"

/**
 * Parse the compact R payload into a validated traversal input.
 *
 * Throws on missing fields, edge/map length disagreement, or invalid root
 * policy; no traversal state exists when failure occurs.
 */
V35Input pst_v35_parse_input(const Rcpp::List& payload);

/**
 * Parse a public simmap tree directly into compact C++ traversal storage.
 *
 * The caller supplies the canonical state lookup and root policy. The returned
 * value preserves 1-based biological ids and owns every edge/map vector needed
 * after the R tree is no longer inspected.
 */
V35Input pst_v35_parse_tree_input(
    const Rcpp::List& tree,
    const Rcpp::DataFrame& state_lookup,
    const Rcpp::IntegerVector& root_anchor_state_ids,
    bool root_is_synthetic);

#endif
