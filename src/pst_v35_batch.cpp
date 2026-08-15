#include "pst_v35_batch.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

/** Convert the public time resolution to one finite positive key scale. */
double pst_v35_batch_time_scale(double time_tolerance) {
  if (!std::isfinite(time_tolerance) || time_tolerance <= 0.0) {
    Rcpp::stop("V35 time tolerance must be finite and positive");
  }
  double time_scale = 1.0 / time_tolerance;
  if (!std::isfinite(time_scale) || time_scale <= 0.0) {
    Rcpp::stop("V35 time tolerance produces an invalid batch scale");
  }
  return time_scale;
}

/**
 * Quantize one root-relative event time for deterministic same-time grouping.
 * Finite input returns an integer key; callers validate biological time ranges.
 */
long long pst_v35_batch_time_key(double time, double time_scale) {
  if (!std::isfinite(time) || !std::isfinite(time_scale) || time_scale <= 0.0) {
    Rcpp::stop("V35 batch time and scale must be finite with positive scale");
  }
  long double scaled_time = static_cast<long double>(time) *
    static_cast<long double>(time_scale);
  if (scaled_time < static_cast<long double>(std::numeric_limits<long long>::min()) ||
      scaled_time > static_cast<long double>(std::numeric_limits<long long>::max())) {
    Rcpp::stop("V35 batch time exceeds the integer-key range");
  }
  return static_cast<long long>(std::llrint(time * time_scale));
}

/**
 * Normalize one event time through the exact transactional key representation.
 * This prevents public TT rows from exposing sub-tolerance floating noise.
 */
double pst_v35_batch_round_time(double time, double time_scale) {
  return static_cast<double>(pst_v35_batch_time_key(time, time_scale)) /
    time_scale;
}

/**
 * Choose a stable queue tie-breaker from one nonempty descendant-tip set.
 * Tip order has no biological meaning; this value only makes output repeatable.
 */
int pst_v35_batch_min_tip(const std::vector<int>& tips) {
  // Scheduling an edge without descendants violates topology ownership and
  // would make deterministic ordering undefined.
  if (tips.empty()) {
    Rcpp::stop("V35 batch scheduler received an empty descendant-tip set");
  }
  return *std::min_element(tips.begin(), tips.end());
}

/**
 * Select the earliest pending event without mutating queue state.
 *
 * Every comparison applies the same total order; independent same-time events
 * therefore cannot acquire iteration-order-dependent ownership.
 */
int pst_v35_batch_next_event_index(
    const std::vector<V35BatchEvent>& queue,
    double time_scale) {
  // A traversal loop calls this helper only while work remains. Reject an empty
  // queue explicitly so no caller can index an invented event zero.
  if (queue.empty()) {
    Rcpp::stop("V35 batch scheduler cannot select from an empty queue");
  }
  int best = 0;
  // Compare every pending event with the current winner. After each iteration
  // `best` identifies the first event under the stable transactional order.
  for (int event_id = 1; event_id < static_cast<int>(queue.size()); ++event_id) {
    const V35BatchEvent& candidate = queue[event_id];
    const V35BatchEvent& incumbent = queue[best];
    long long candidate_time_key = pst_v35_batch_time_key(
      candidate.time, time_scale
    );
    long long incumbent_time_key = pst_v35_batch_time_key(
      incumbent.time, time_scale
    );
    // Earlier times win. Within one time key, causal phase precedes descendant
    // tip and phylogeny edge so incoming changes, cladogenesis, zero-offset
    // daughters, and terminals cannot be reordered by raw edge storage.
    if (candidate_time_key < incumbent_time_key ||
        (candidate_time_key == incumbent_time_key &&
         candidate.phase < incumbent.phase) ||
        (candidate_time_key == incumbent_time_key &&
         candidate.phase == incumbent.phase &&
         candidate.priority_tip < incumbent.priority_tip) ||
        (candidate_time_key == incumbent_time_key &&
         candidate.phase == incumbent.phase &&
         candidate.priority_tip == incumbent.priority_tip &&
         candidate.edge_id < incumbent.edge_id)) {
      best = event_id;
    }
  }
  return best;
}

/**
 * Move pending causal ownership from a closed scenario edge to its replacement.
 * Frozen batch origins remain unchanged so pre-batch source support is stable.
 */
void pst_v35_batch_reassign_scenario(
    std::vector<V35BatchEvent>& queue,
    int from_edge_id,
    int to_edge_id) {
  // Visit all pending biological boundaries; after each iteration only events
  // owned by the replaced scenario have changed mutable ownership.
  for (V35BatchEvent& event : queue) {
    // Unrelated scenarios retain their owners. Matching events move forward to
    // the committed child but keep `batch_origin_scenario_edge_id` frozen.
    if (event.scenario_edge_id == from_edge_id) {
      event.scenario_edge_id = to_edge_id;
    }
  }
}

/**
 * Initialize one empty same-time source-removal transaction.
 *
 * The integer key and rounded time are immutable. Capacity/removal maps start
 * empty because sources are discovered lazily in causal queue order.
 */
V35BatchTransaction pst_v35_batch_begin(double time, double time_scale) {
  return V35BatchTransaction{
    pst_v35_batch_time_key(time, time_scale),
    pst_v35_batch_round_time(time, time_scale),
    std::map<int, int>(),
    std::map<int, int>(),
    0,
    0
  };
}

/**
 * Plan and validate one source lineage withdrawal before destination mutation.
 *
 * Repeated events from the same scenario must report the same frozen capacity.
 * The cumulative withdrawal count may never exceed that capacity. Phase order
 * is monotonic across the complete transaction.
 */
void pst_v35_batch_plan_source_removal(
    V35BatchTransaction& transaction,
    int scenario_edge_id,
    int source_capacity,
    int phase) {
  // Every biological event consumes one positive source scenario lineage.
  if (scenario_edge_id < 1 || source_capacity < 1) {
    Rcpp::stop("V35 batch source removal has invalid scenario capacity");
  }
  // Causal phases are closed in ascending order; returning to an earlier phase
  // would make later handlers inspect an incomplete virtual predecessor state.
  if (phase < V35_BATCH_INCOMING_ANAGENETIC ||
      phase > V35_BATCH_TERMINAL ||
      phase < transaction.last_phase) {
    Rcpp::stop("V35 same-time event violated causal phase order");
  }
  std::map<int, int>::iterator known_capacity =
    transaction.source_capacity.find(scenario_edge_id);
  // The first removal freezes pre-phase capacity for this source scenario.
  if (known_capacity == transaction.source_capacity.end()) {
    transaction.source_capacity[scenario_edge_id] = source_capacity;
  } else {
    // Later simultaneous removals must observe the exact same frozen capacity.
    if (known_capacity->second != source_capacity) {
      Rcpp::stop("V35 same-time source capacity changed within a batch");
    }
  }
  int planned = ++transaction.source_removals[scenario_edge_id];
  // All withdrawals are calculated before the current handler creates any
  // destination; exceeding capacity proves the source would be overdrawn.
  if (planned > transaction.source_capacity[scenario_edge_id]) {
    std::ostringstream diagnostic;
    diagnostic
      << "V35 same-time events overdraw one source scenario"
      << " [time_key=" << transaction.time_key
      << ", time=" << transaction.time
      << ", scenario_edge_id=" << scenario_edge_id
      << ", capacity=" << transaction.source_capacity[scenario_edge_id]
      << ", planned=" << planned
      << ", phase=" << phase
      << "]";
    Rcpp::stop(diagnostic.str());
  }
  ++transaction.event_count;
  transaction.last_phase = phase;
}

/**
 * Validate one fully planned batch before its final frontier is serialized.
 *
 * Every encountered source must have positive, bounded removals and at least
 * one event must have advanced the transaction.
 */
void pst_v35_batch_validate(const V35BatchTransaction& transaction) {
  // Empty transactions cannot own a public TT/scenario-matrix row.
  if (transaction.event_count < 1) {
    Rcpp::stop("V35 attempted to commit an empty same-time transaction");
  }
  // Check each frozen source after all causal phases. Every processed source
  // retains a positive bounded removal count.
  for (const std::pair<const int, int>& item : transaction.source_capacity) {
    std::map<int, int>::const_iterator removals =
      transaction.source_removals.find(item.first);
    // A frozen capacity without a planned withdrawal indicates incomplete
    // event planning; an excess withdrawal would violate source conservation.
    if (removals == transaction.source_removals.end() ||
        removals->second < 1 ||
        removals->second > item.second) {
      Rcpp::stop("V35 same-time source-removal validation failed");
    }
  }
}

/**
 * Freeze endpoint scenario/state/path ownership before causal phase handling.
 *
 * The input vectors are one-based and parallel. Inactive tips may carry NA
 * scenario ownership and are omitted; every active row receives one immutable
 * source identity and an initially unchanged survivor intent.
 */
V35ScenarioBatchIntent pst_v35_scenario_batch_begin(
    long long time_key,
    double time,
    const std::vector<int>& current_scenario_edge_id,
    const std::vector<int>& current_state_id,
    const std::vector<int>& current_path_id,
    const std::vector<int>& current_phylo_edge_id,
    const std::vector<int>& previously_terminal_tip_ids) {
  if (current_scenario_edge_id.size() != current_state_id.size() ||
      current_scenario_edge_id.size() != current_path_id.size() ||
      current_scenario_edge_id.size() != current_phylo_edge_id.size()) {
    Rcpp::stop("V35 scenario batch frontier vectors are not parallel");
  }
  V35ScenarioBatchIntent intent;
  intent.time_key = time_key;
  intent.time = time;
  std::set<int> previously_terminal(
    previously_terminal_tip_ids.begin(),
    previously_terminal_tip_ids.end()
  );
  for (int tip_id = 1;
       tip_id < static_cast<int>(current_scenario_edge_id.size());
       ++tip_id) {
    if (previously_terminal.count(tip_id)) {
      continue;
    }
    int scenario_edge_id = current_scenario_edge_id[tip_id];
    if (scenario_edge_id == NA_INTEGER) {
      continue;
    }
    int state_id = current_state_id[tip_id];
    int path_id = current_path_id[tip_id];
    int phylo_edge_id = current_phylo_edge_id[tip_id];
    if (scenario_edge_id < 1 || state_id < 1 || path_id < 1 ||
        phylo_edge_id < 1) {
      Rcpp::stop("V35 active scenario batch member has invalid ownership");
    }
    intent.tips.push_back(V35ScenarioBatchTipIntent{
      tip_id,
      scenario_edge_id,
      state_id,
      path_id,
      phylo_edge_id,
      false,
      state_id,
      path_id,
      phylo_edge_id
    });
  }
  return intent;
}

/** Return the unique mutable intent row for a valid active biological tip. */
V35ScenarioBatchTipIntent& pst_v35_scenario_batch_tip(
    V35ScenarioBatchIntent& intent,
    int tip_id) {
  V35ScenarioBatchTipIntent* match = nullptr;
  for (V35ScenarioBatchTipIntent& tip : intent.tips) {
    if (tip.tip_id != tip_id) {
      continue;
    }
    if (match != nullptr) {
      Rcpp::stop("V35 scenario batch contains a duplicate tip row");
    }
    match = &tip;
  }
  if (match == nullptr) {
    Rcpp::stop("V35 scenario batch references an unknown active tip");
  }
  return *match;
}

/** Update final survivor history while keeping the frozen source immutable. */
void pst_v35_scenario_batch_update_survivor(
    V35ScenarioBatchIntent& intent,
    int tip_id,
    int state_id,
    int path_id) {
  if (state_id < 1 || path_id < 1) {
    Rcpp::stop("V35 scenario survivor intent has invalid state or path");
  }
  V35ScenarioBatchTipIntent& tip =
    pst_v35_scenario_batch_tip(intent, tip_id);
  if (tip.terminal) {
    Rcpp::stop("V35 terminal scenario member cannot receive survivor history");
  }
  tip.survivor_state_id = state_id;
  tip.survivor_path_id = path_id;
}

/** Mark a frozen endpoint member terminal without changing its source ID. */
void pst_v35_scenario_batch_mark_terminal(
    V35ScenarioBatchIntent& intent,
    int tip_id) {
  V35ScenarioBatchTipIntent& tip =
    pst_v35_scenario_batch_tip(intent, tip_id);
  tip.terminal = true;
  intent.affected_source_scenario_ids.insert(tip.source_scenario_id);
}

/** Mark one valid frozen source scenario for the post-phase commit. */
void pst_v35_scenario_batch_mark_source_affected(
    V35ScenarioBatchIntent& intent,
    int scenario_edge_id) {
  if (scenario_edge_id < 1) {
    Rcpp::stop("V35 scenario batch affected source ID is invalid");
  }
  bool found = false;
  for (const V35ScenarioBatchTipIntent& tip : intent.tips) {
    if (tip.source_scenario_id == scenario_edge_id) {
      found = true;
      break;
    }
  }
  if (!found) {
    Rcpp::stop("V35 scenario batch affected source has no endpoint member");
  }
  intent.affected_source_scenario_ids.insert(scenario_edge_id);
}
