#ifndef PST_V35_TRAVERSAL_HPP
#define PST_V35_TRAVERSAL_HPP

#include "pst_v35_types.hpp"

/** Return structural input counts without scheduling biological events. */
Rcpp::List pst_v35_traversal_smoke_summary(const V35Input& input);

/**
 * Execute chronological V35 traversal and return its committed online state.
 *
 * Output owns scenario segments, terminal records, canonical paths, TT/SS
 * counters, and materialization metadata. Normal execution throws on failure;
 * snapshot-mode recovery can instead serialize the checkpoint preceding the
 * failed batch.
 */
Rcpp::List pst_v35_extract_active_summary(
  const V35Input& input,
  bool recover_debug_failure = false,
  int debug_fail_batch_id = NA_INTEGER,
  double time_tolerance = 1e-6);

#endif
