#ifndef PST_V35_DEBUG_HPP
#define PST_V35_DEBUG_HPP

#include <Rcpp.h>

/**
 * Clip one rooted simmap tree at an absolute root-relative time.
 *
 * Inputs: `tree` is a V35 public phylogeny/scenario/transition tree and
 * `cutoff_time` is a committed TT time. Output: a valid rooted simmap whose
 * reached topology, map prefixes, size maps, and path-size maps end at the
 * cutoff. Active frontier edges terminate in deterministic provisional tips.
 * Ownership: the returned R list owns all new matrices and vectors; `tree` is
 * read-only. Failure: malformed topology, negative cutoffs, or inconsistent map
 * lengths raise an R error before a partial tree is returned.
 */
Rcpp::List pst_v35_debug_clip_tree_cpp(
  const Rcpp::List& tree,
  double cutoff_time,
  double time_tolerance
);

#endif
