#ifndef PST_V34_TYPES_HPP
#define PST_V34_TYPES_HPP

#include <Rcpp.h>
#include <vector>

/**
 * One contiguous mapped-state segment on an input phylogeny edge.
 *
 * Ownership: values are copied from the R simmap into C++ traversal storage.
 * Invariant: ids are positive and 1-based, duration is nonnegative, and
 * `end_offset` is the cumulative distance from the edge's parent endpoint.
 */
struct EdgeMapSegment {
  int edge_id;
  int state_id;
  double duration;
  double end_offset;
};

/**
 * Validated, compact ownership boundary for one V34 traversal.
 *
 * Topology vectors are edge parallel; flattened map vectors are partitioned by
 * `map_edge_offset`; state and root ids remain in the public 1-based id space.
 * Rcpp vectors own protected R storage for the duration of the C++ bridge call.
 * Parse functions reject inconsistent lengths or unknown state labels before
 * constructing this value.
 */
struct V34Input {
  int ntips;
  int nedge;
  Rcpp::IntegerVector edge_parent;
  Rcpp::IntegerVector edge_child;
  Rcpp::NumericVector edge_length;
  Rcpp::IntegerVector map_edge_id;
  Rcpp::IntegerVector map_state_id;
  Rcpp::NumericVector map_duration;
  Rcpp::IntegerVector map_edge_offset;
  Rcpp::CharacterVector tip_label;
  Rcpp::DataFrame state_lookup;
  Rcpp::IntegerVector root_anchor_state_ids;
  // A synthetic root is bookkeeping rather than an observed state. Online
  // PUNIQ/CTT scoring suppresses its initial daughter paths but retains real
  // departures from an explicitly supplied biological root.
  bool root_is_synthetic;
};

#endif
