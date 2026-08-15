#include "pst_v34_traversal.hpp"
#include "pst_v34_batch.hpp"
#include "pst_v34_state.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cfloat>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <set>
#include <string>
#include <utility>
#include <vector>

/**
 * Return cheap structural counts for bridge smoke tests without running traversal.
 *
 * Input: a validated compact traversal input. Output: counts only; no
 * biological event is scheduled and no mutable traversal state is created.
 */
Rcpp::List pst_v34_traversal_smoke_summary(const V34Input& input) {
  // Return cheap structural counts used by Phase 2 tests.
  return Rcpp::List::create(
    Rcpp::Named("ntips") = input.ntips,
    Rcpp::Named("nedge") = input.nedge,
    Rcpp::Named("nsegments") = static_cast<int>(input.map_edge_id.size()),
    Rcpp::Named("nstates") = static_cast<int>(input.state_lookup.nrows()),
    Rcpp::Named("root_anchor_state_ids") = input.root_anchor_state_ids
  );
}

namespace {

enum class V34ScenarioEdgeType {
  Pending,
  Root,
  Stay,
  Leave
};

struct V34ScenarioEdgeRecord {
  std::vector<double> end_times;
  std::vector<double> durations;
  std::vector<int> state_ids;
  std::vector<int> path_ids;
  std::vector<int> sizes;
  std::vector<int> parent_nodes;
  std::vector<int> child_nodes;
  bool is_terminal;
  int initial_path_id;
  int terminal_path_id;
  int parent_scenario_edge_id;
  double formation_time;
  V34ScenarioEdgeType edge_type;
};

/**
 * Compact diagnostic identity for one committed same-time batch.
 *
 * This record is allocated only when snapshot recovery is enabled. Event ids,
 * types, phases, phylogeny edges, and source scenarios are position-aligned in
 * causal execution order. Terminal ids contain affected biological tips for
 * terminal events. The record owns diagnostics only and never participates in
 * PST, TT, SS, path, scenario, or scoring decisions.
 */
struct V34DebugStepRecord {
  std::vector<std::string> event_types;
  std::vector<std::string> causal_phases;
  std::vector<int> event_ids;
  std::vector<int> phylogeny_edges;
  std::vector<int> scenario_ids;
  std::vector<int> terminal_ids;
};

/**
 * Return the public diagnostic label for one transactional causal phase.
 *
 * Valid phase ids map to stable names used by debug snapshots. An invalid id
 * indicates a scheduler defect and fails before the corresponding batch can be
 * published. This formatting helper owns no traversal state.
 */
std::string v34_debug_phase_label(int phase) {
  // Incoming-edge state changes precede topology mutation at the same time.
  if (phase == V34_BATCH_INCOMING_ANAGENETIC) {
    return "incoming_edge_anagenetic";
  }
  // Cladogenesis removes one parent lineage and creates daughter ownership.
  if (phase == V34_BATCH_CLADOGENESIS) {
    return "cladogenesis";
  }
  // Zero-offset daughter transitions inspect the virtual post-split frontier.
  if (phase == V34_BATCH_OUTGOING_ZERO_OFFSET) {
    return "outgoing_edge_anagenetic";
  }
  // Terminal closure retains endpoint support in this row and removes it later.
  if (phase == V34_BATCH_TERMINAL) {
    return "terminal";
  }
  Rcpp::stop("V34 debug metadata received an invalid causal phase");
  return std::string();
}

/**
 * Append one integer identity only when it is not already represented.
 *
 * Debug scenario and terminal collections describe set membership while
 * preserving first causal appearance. The helper mutates only the supplied
 * diagnostic vector and cannot affect traversal ownership.
 */
void v34_debug_append_unique(std::vector<int>& values, int value) {
  // Missing ids do not identify a public biological or scenario object.
  if (value == NA_INTEGER) {
    return;
  }
  // Preserve the first causal occurrence and suppress repeated membership from
  // simultaneous events that reference the same source or terminal lineage.
  if (std::find(values.begin(), values.end(), value) == values.end()) {
    values.push_back(value);
  }
}

/**
 * Debug-only copy of every structure a same-time transaction may mutate.
 *
 * Normal traversal never allocates this object. Snapshot recovery creates one
 * before a batch, discards it after successful validation, and restores it if
 * any phase throws. Public slices need no copy because they are appended only
 * after validation.
 */
struct V34DebugBatchCheckpoint {
  std::vector<V34ScenarioEdgeRecord> records;
  int next_edge_id;
  int next_node_id;
  std::vector<int> current_state_id;
  std::vector<int> current_edge_node_id;
  std::vector<int> current_scenario_edge_id;
  std::vector<int> current_phylo_edge_id;
  std::vector<int> current_path_id;
  std::vector<std::vector<int> > path_components;
  std::vector<int> path_parent;
  std::vector<int> path_added_state_id;
  std::vector<std::vector<int> > transition_to_path;
  V34LiveTTState online_counters;
  std::vector<int> terminal_tip_ids;
  std::vector<int> terminal_phylo_edge_ids;
  std::vector<int> terminal_scenario_edge_ids;
  std::vector<int> terminal_path_ids;
  std::vector<double> terminal_times;
  std::vector<V34BatchEvent> queue;

  /**
   * Copy the complete pre-batch mutation surface for debug rollback.
   *
   * Every argument is traversal-owned committed state. The checkpoint owns its
   * copies, preserves all vector/id alignment invariants, and cannot fail after
   * ordinary container allocation succeeds. Normal nondebug traversal never
   * invokes this constructor.
   */
  V34DebugBatchCheckpoint(
      const std::vector<V34ScenarioEdgeRecord>& records_value,
      int next_edge_id_value,
      int next_node_id_value,
      const std::vector<int>& current_state_id_value,
      const std::vector<int>& current_edge_node_id_value,
      const std::vector<int>& current_scenario_edge_id_value,
      const std::vector<int>& current_phylo_edge_id_value,
      const std::vector<int>& current_path_id_value,
      const std::vector<std::vector<int> >& path_components_value,
      const std::vector<int>& path_parent_value,
      const std::vector<int>& path_added_state_id_value,
      const std::vector<std::vector<int> >& transition_to_path_value,
      const V34LiveTTState& online_counters_value,
      const std::vector<int>& terminal_tip_ids_value,
      const std::vector<int>& terminal_phylo_edge_ids_value,
      const std::vector<int>& terminal_scenario_edge_ids_value,
      const std::vector<int>& terminal_path_ids_value,
      const std::vector<double>& terminal_times_value,
      const std::vector<V34BatchEvent>& queue_value)
      : records(records_value),
        next_edge_id(next_edge_id_value),
        next_node_id(next_node_id_value),
        current_state_id(current_state_id_value),
        current_edge_node_id(current_edge_node_id_value),
        current_scenario_edge_id(current_scenario_edge_id_value),
        current_phylo_edge_id(current_phylo_edge_id_value),
        current_path_id(current_path_id_value),
        path_components(path_components_value),
        path_parent(path_parent_value),
        path_added_state_id(path_added_state_id_value),
        transition_to_path(transition_to_path_value),
        online_counters(online_counters_value),
        terminal_tip_ids(terminal_tip_ids_value),
        terminal_phylo_edge_ids(terminal_phylo_edge_ids_value),
        terminal_scenario_edge_ids(terminal_scenario_edge_ids_value),
        terminal_path_ids(terminal_path_ids_value),
        terminal_times(terminal_times_value),
        queue(queue_value) {}
};

/**
 * Allocate one scenario edge and its child topology node at event time.
 *
 * `records`, `next_edge_id`, and `next_node_id` are mutated together. The new
 * record starts with one zero-duration `state/path/size` segment and remains
 * terminal until a later split closes it. The returned edge id is positive and
 * indexes the resized record vector; allocation failure or invalid upstream ids
 * aborts the enclosing batch before commit.
 */
int v34_new_scenario_edge(std::vector<V34ScenarioEdgeRecord>& records,
                          int& next_edge_id,
                          int& next_node_id,
                          int state_id,
                          int size,
                          int parent_node,
                          double time,
                          int path_id,
                          int parent_scenario_edge_id) {
  int edge_id = next_edge_id++;
  int child_node = next_node_id++;
  // Root status is structural at allocation. Non-root edges remain pending
  // until every simultaneous event has resolved; assigning Stay/Leave here
  // would make biological classification depend on traversal iteration order.
  V34ScenarioEdgeType edge_type = parent_scenario_edge_id < 1 ?
    V34ScenarioEdgeType::Root : V34ScenarioEdgeType::Pending;
  // Scenario records are 1-based. Grow storage before allocating a newly
  // numbered edge so topology and record indexes remain parallel.
  if (edge_id >= static_cast<int>(records.size())) {
    records.resize(edge_id + 1);
  }
  V34ScenarioEdgeRecord rec;
  rec.end_times.push_back(time);
  rec.durations.push_back(0.0);
  rec.state_ids.push_back(state_id);
  rec.path_ids.push_back(path_id);
  rec.sizes.push_back(size);
  rec.parent_nodes.push_back(parent_node);
  rec.child_nodes.push_back(child_node);
  rec.is_terminal = true;
  rec.initial_path_id = path_id;
  rec.terminal_path_id = path_id;
  rec.parent_scenario_edge_id = parent_scenario_edge_id;
  rec.formation_time = time;
  rec.edge_type = edge_type;
  records[edge_id] = rec;
  return edge_id;
}

/**
 * Classify every scenario edge whose simultaneous formation batch has resolved.
 *
 * Root edges are fixed at allocation. A pending non-root edge compares its
 * initial path with its parent's terminal path only after all events at the
 * timestamp have been processed. The function mutates classification metadata
 * only; LDIF and TT accounting remain owned by the atomic counter commit.
 */
void v34_classify_pending_scenario_edges(
    std::vector<V34ScenarioEdgeRecord>& records) {
  // Creation order is parent-first, so parent terminal paths are available
  // before every child is inspected.
  for (int edge_id = 1; edge_id < static_cast<int>(records.size()); ++edge_id) {
    V34ScenarioEdgeRecord& record = records[edge_id];
    // Root and already resolved edges retain their biological formation class.
    if (record.edge_type != V34ScenarioEdgeType::Pending) {
      continue;
    }
    int parent_id = record.parent_scenario_edge_id;
    if (parent_id < 1 || parent_id >= static_cast<int>(records.size()) ||
        record.initial_path_id < 1 ||
        records[parent_id].terminal_path_id < 1) {
      Rcpp::stop("V34 pending scenario edge has invalid path ancestry");
    }
    record.edge_type =
      record.initial_path_id == records[parent_id].terminal_path_id ?
        V34ScenarioEdgeType::Stay : V34ScenarioEdgeType::Leave;
  }
}

/** Return the stable public label for one resolved traversal edge class. */
std::string v34_scenario_edge_type_label(V34ScenarioEdgeType edge_type) {
  if (edge_type == V34ScenarioEdgeType::Root) return "Root";
  if (edge_type == V34ScenarioEdgeType::Stay) return "Stay";
  if (edge_type == V34ScenarioEdgeType::Leave) return "Leave";
  Rcpp::stop("V34 scenario edge remained pending after batch resolution");
  return std::string();
}

/**
 * Close the currently open scenario segment at a nonterminal event boundary.
 *
 * The record and time are traversal-owned. Duration is computed from the open
 * segment's start marker, its endpoint becomes `time`, and `is_terminal` is
 * cleared because daughters or a changed continuation will own later exposure.
 * The caller guarantees a valid edge id and a nonempty segment vector.
 */
void v34_close_active_segment(std::vector<V34ScenarioEdgeRecord>& records,
                              int scenario_edge_id,
                              double time) {
  V34ScenarioEdgeRecord& rec = records[scenario_edge_id];
  int last = static_cast<int>(rec.end_times.size()) - 1;
  double start_time = rec.end_times[last];
  rec.durations[last] = time - start_time;
  rec.end_times[last] = time;
  rec.is_terminal = false;
}

/**
 * Extend and close the current scenario segment at a lineage terminal time.
 *
 * Endpoint inclusion may call this after the segment already has positive
 * duration. The function preserves its original start, extends only to the
 * later of the stored or supplied endpoint, and marks the edge terminal. The
 * caller guarantees valid scenario ownership; no public state is emitted here.
 */
void v34_close_terminal_segment(std::vector<V34ScenarioEdgeRecord>& records,
                                int scenario_edge_id,
                                double time) {
  V34ScenarioEdgeRecord& rec = records[scenario_edge_id];
  int last = static_cast<int>(rec.end_times.size()) - 1;
  double previous_end = rec.end_times[last];
  double previous_duration = rec.durations[last];
  double start_time = previous_end - previous_duration;
  double end_time = std::max(previous_end, time);
  rec.durations[last] = end_time - start_time;
  rec.end_times[last] = end_time;
  rec.is_terminal = true;
}

/**
 * Start a new open segment on an existing scenario edge after a committed event.
 *
 * State, path, size, and start time are appended in parallel, with zero initial
 * duration. The segment is provisionally terminal until another event closes
 * or extends it. The caller owns validation and invokes this only inside the
 * current batch transaction.
 */
void v34_append_segment(std::vector<V34ScenarioEdgeRecord>& records,
                        int scenario_edge_id,
                        int state_id,
                        int size,
                        double time,
                        int path_id) {
  V34ScenarioEdgeRecord& rec = records[scenario_edge_id];
  rec.end_times.push_back(time);
  rec.durations.push_back(0.0);
  rec.state_ids.push_back(state_id);
  rec.path_ids.push_back(path_id);
  rec.sizes.push_back(size);
  rec.is_terminal = true;
  rec.terminal_path_id = path_id;
}

/**
 * Commit final stable scenario identities after all same-key causal phases.
 *
 * Phase handlers may have changed state/path/phylogeny ownership, but the
 * supplied scenario records and tip scenario IDs are the frozen pre-batch
 * topology. This function is the only point that ends, continues, or splits
 * those stable scenarios for the completed survivor partition.
 */
void v34_commit_scenario_batch(
    V34ScenarioBatchIntent& intent,
    std::vector<V34ScenarioEdgeRecord>& records,
    int& next_edge_id,
    int& next_node_id,
    std::vector<int>& current_scenario_edge_id,
    const std::vector<int>& current_state_id,
    const std::vector<int>& current_path_id,
    const std::vector<int>& current_phylo_edge_id,
    std::vector<V34BatchEvent>& queue,
    V34LiveTTState& online_counters) {
  if (current_scenario_edge_id.size() != current_state_id.size() ||
      current_scenario_edge_id.size() != current_path_id.size() ||
      current_scenario_edge_id.size() != current_phylo_edge_id.size()) {
    Rcpp::stop("V34 scenario commit frontier vectors are not parallel");
  }

  typedef std::pair<int, int> V34StatePathKey;
  std::map<
    int,
    std::map<V34StatePathKey, std::vector<int> >
  > survivor_groups_by_source;
  std::map<int, std::set<int> > endpoint_phylo_edges_by_source;

  // Sync the final survivor histories from the completed causal frontier while
  // retaining each row's immutable pre-batch source identity.
  for (V34ScenarioBatchTipIntent& tip : intent.tips) {
    endpoint_phylo_edges_by_source[tip.source_scenario_id].insert(
      tip.endpoint_phylo_edge_id
    );
    if (tip.terminal) {
      current_scenario_edge_id[tip.tip_id] = tip.source_scenario_id;
      continue;
    }
    int state_id = current_state_id[tip.tip_id];
    int path_id = current_path_id[tip.tip_id];
    if (state_id < 1 || path_id < 1) {
      Rcpp::stop("V34 scenario commit survivor has invalid state or path");
    }
    tip.survivor_state_id = state_id;
    tip.survivor_path_id = path_id;
    tip.survivor_phylo_edge_id = current_phylo_edge_id[tip.tip_id];
    survivor_groups_by_source[tip.source_scenario_id][
      V34StatePathKey(state_id, path_id)
    ].push_back(tip.tip_id);
  }

  std::set<int> committed_sources;
  for (int source_scenario_id : intent.affected_source_scenario_ids) {
    if (!committed_sources.insert(source_scenario_id).second) {
      Rcpp::stop("V34 scenario batch committed one source more than once");
    }
    if (source_scenario_id < 1 ||
        source_scenario_id >= static_cast<int>(records.size()) ||
        records[source_scenario_id].sizes.empty() ||
        records[source_scenario_id].state_ids.empty() ||
        records[source_scenario_id].path_ids.empty()) {
      Rcpp::stop("V34 scenario batch source record is invalid");
    }
    int endpoint_count = static_cast<int>(
      endpoint_phylo_edges_by_source[source_scenario_id].size()
    );
    if (endpoint_count < 1 ||
        records[source_scenario_id].sizes.back() != endpoint_count) {
      Rcpp::stop("V34 scenario batch endpoint support disagrees with its source");
    }

    std::map<V34StatePathKey, std::vector<int> >& survivor_groups =
      survivor_groups_by_source[source_scenario_id];
    if (survivor_groups.empty()) {
      v34_close_terminal_segment(records, source_scenario_id, intent.time);
      continue;
    }

    V34ScenarioEdgeRecord& source = records[source_scenario_id];
    v34_close_active_segment(records, source_scenario_id, intent.time);
    if (survivor_groups.size() == 1) {
      const std::pair<const V34StatePathKey, std::vector<int> >& group =
        *survivor_groups.begin();
      std::set<int> survivor_phylo_edges;
      for (int tip_id : group.second) {
        survivor_phylo_edges.insert(current_phylo_edge_id[tip_id]);
      }
      v34_append_segment(
        records,
        source_scenario_id,
        group.first.first,
        static_cast<int>(survivor_phylo_edges.size()),
        intent.time,
        group.first.second
      );
      for (int tip_id : group.second) {
        current_scenario_edge_id[tip_id] = source_scenario_id;
      }
      continue;
    }

    int parent_node = source.child_nodes.back();
    for (const std::pair<const V34StatePathKey, std::vector<int> >& group :
         survivor_groups) {
      std::set<int> survivor_phylo_edges;
      for (int tip_id : group.second) {
        survivor_phylo_edges.insert(current_phylo_edge_id[tip_id]);
      }
      int child_scenario_id = v34_new_scenario_edge(
        records,
        next_edge_id,
        next_node_id,
        group.first.first,
        static_cast<int>(survivor_phylo_edges.size()),
        parent_node,
        intent.time,
        group.first.second,
        source_scenario_id
      );
      for (int tip_id : group.second) {
        current_scenario_edge_id[tip_id] = child_scenario_id;
      }
    }
  }

  // Scenario-dependent CTT identity is published only after every stable child
  // exists. LDIF then observes the same exact raw scenario IDs at row capture.
  for (const V34PendingScenarioTransition& transition : intent.transitions) {
    int tip_id = transition.representative_tip_id;
    if (tip_id < 1 ||
        tip_id >= static_cast<int>(current_scenario_edge_id.size()) ||
        current_scenario_edge_id[tip_id] < 1) {
      Rcpp::stop("V34 pending transition has no committed destination scenario");
    }
    online_counters.record_ctt_event(
      transition.time,
      transition.from_state_id,
      transition.to_state_id,
      transition.phylo_edge_id,
      transition.source_path_id,
      transition.destination_path_id,
      transition.source_scenario_id,
      current_scenario_edge_id[tip_id],
      transition.phylo_edge_is_terminal
    );
  }

  // Every queued boundary now lies after this completed key. Resolve both its
  // mutable and next-batch frozen owner from the surviving representative.
  for (V34BatchEvent& event : queue) {
    int tip_id = event.priority_tip;
    if (tip_id < 1 ||
        tip_id >= static_cast<int>(current_scenario_edge_id.size()) ||
        current_scenario_edge_id[tip_id] < 1) {
      Rcpp::stop("V34 future event has no committed scenario owner");
    }
    event.scenario_edge_id = current_scenario_edge_id[tip_id];
    event.batch_origin_scenario_edge_id = event.scenario_edge_id;
  }
}

/**
 * Group final tips by their committed root-to-tip scenario trajectory.
 *
 * Inputs are live scenario records and final tip ownership. The key uses
 * state ids, canonical path ids, and durations rounded to the caller's public
 * scenario cleanup tolerance. Lineage-size labels are deliberately ignored.
 * Output is one positive group id per tip. Invalid ownership or cyclic topology
 * stops before tree formatting can merge unrelated biological trajectories.
 */
std::vector<int> v34_terminal_trajectory_groups(
    const std::vector<V34ScenarioEdgeRecord>& records,
    const std::vector<int>& final_scenario_edge_by_tip,
    int ntips,
    double time_tolerance) {
  std::map<int, int> edge_by_child_node;
  // Index every committed scenario edge by its topology child so terminal
  // paths can move rootward without scanning a materialized public tree.
  for (int edge_id = 1; edge_id < static_cast<int>(records.size()); ++edge_id) {
    const V34ScenarioEdgeRecord& record = records[edge_id];
    // Provisional empty records own no topology and are excluded from ancestry.
    if (record.child_nodes.empty()) {
      continue;
    }
    edge_by_child_node[record.child_nodes.back()] = edge_id;
  }

  std::map<std::string, int> group_by_key;
  std::vector<int> group_by_tip(ntips + 1, 0);
  int next_group_id = 1;
  // Build one committed trajectory key for every biological tip. After each
  // iteration the tip owns exactly one positive equivalence group.
  for (int tip_id = 1; tip_id <= ntips; ++tip_id) {
    // Final ownership vectors are 1-based and must contain every biological tip.
    if (tip_id >= static_cast<int>(final_scenario_edge_by_tip.size())) {
      Rcpp::stop("V34 final scenario ownership is not tip-aligned");
    }
    int edge_id = final_scenario_edge_by_tip[tip_id];
    std::vector<int> reverse_path;
    std::set<int> visited_edges;
    // Walk the committed scenario ancestry from this tip to the synthetic root.
    // Each iteration prepends one unique edge to the eventual root-first path.
    while (edge_id >= 1 && edge_id < static_cast<int>(records.size())) {
      // Revisited edges prove a topology cycle and cannot define a trajectory.
      if (!visited_edges.insert(edge_id).second) {
        Rcpp::stop("V34 scenario record topology contains a cycle");
      }
      reverse_path.push_back(edge_id);
      const V34ScenarioEdgeRecord& record = records[edge_id];
      // A record without a parent endpoint cannot be connected to the root.
      if (record.parent_nodes.empty()) {
        Rcpp::stop("V34 scenario trajectory record has no parent node");
      }
      std::map<int, int>::const_iterator parent_edge =
        edge_by_child_node.find(record.parent_nodes.back());
      // Absence of an incoming edge means the synthetic root was reached.
      if (parent_edge == edge_by_child_node.end()) {
        break;
      }
      edge_id = parent_edge->second;
    }
    // A missing final edge cannot be converted into a biological trajectory.
    if (reverse_path.empty()) {
      Rcpp::stop("V34 tip has no committed terminal scenario trajectory");
    }
    std::reverse(reverse_path.begin(), reverse_path.end());

    std::vector<int> positive_states;
    std::vector<int> positive_paths;
    std::vector<double> positive_durations;
    std::vector<int> zero_state_sequence;
    std::vector<int> zero_path_sequence;
    // Concatenate committed edge segments in root-to-tip biological order.
    // Positive adjacent equal states are compacted; zero-state sequence is
    // retained only as a fallback for wholly zero-duration trajectories.
    for (int path_edge_id : reverse_path) {
      const V34ScenarioEdgeRecord& record = records[path_edge_id];
      // State and duration vectors must remain segment-parallel.
      if (record.state_ids.size() != record.durations.size() ||
          record.path_ids.size() != record.durations.size()) {
        Rcpp::stop("V34 scenario trajectory segments are not aligned");
      }
      // Visit every segment on this scenario edge in chronological order.
      for (int segment_id = 0;
           segment_id < static_cast<int>(record.state_ids.size());
           ++segment_id) {
        int state_id = record.state_ids[segment_id];
        int path_id = record.path_ids[segment_id];
        double duration = record.durations[segment_id];
        // Keep one copy of each adjacent state for a zero-duration fallback key.
        if (zero_state_sequence.empty() || zero_state_sequence.back() != state_id) {
          zero_state_sequence.push_back(state_id);
          zero_path_sequence.push_back(path_id);
        } else if (zero_path_sequence.back() != path_id) {
          zero_state_sequence.push_back(state_id);
          zero_path_sequence.push_back(path_id);
        }
        // Zero-length formatting remnants do not distinguish positive trajectories.
        if (duration <= time_tolerance) {
          continue;
        }
        // Adjacent positive runs in the same state form one continuous interval.
        if (!positive_states.empty() &&
            positive_states.back() == state_id &&
            positive_paths.back() == path_id) {
          positive_durations.back() += duration;
        } else {
          // A changed state starts a new positive trajectory interval.
          positive_states.push_back(state_id);
          positive_paths.push_back(path_id);
          positive_durations.push_back(duration);
        }
      }
    }

    std::ostringstream key;
    // Positive trajectories are keyed by compact state and rounded duration.
    if (!positive_states.empty()) {
      // Emit each compact interval once in biological order.
      for (int interval_id = 0;
           interval_id < static_cast<int>(positive_states.size());
           ++interval_id) {
        // Later intervals are delimited from the preceding state-duration pair.
        if (interval_id > 0) {
          key << "|";
        }
        double duration = positive_durations[interval_id];
        key << positive_states[interval_id] << ":" <<
          positive_paths[interval_id] << ":";
        long double scaled_duration = time_tolerance > 0.0 ?
          static_cast<long double>(duration) /
            static_cast<long double>(time_tolerance) :
          std::numeric_limits<long double>::infinity();
        // Ordinary tolerances use a stable integer bin. When the requested
        // threshold is zero or too small for a 64-bit bin, exact double text
        // preserves distinct durations without overflow.
        if (std::isfinite(scaled_duration) &&
            scaled_duration <= static_cast<long double>(
              std::numeric_limits<long long>::max()
            )) {
          key << "q" << static_cast<long long>(std::llrint(scaled_duration));
        } else {
          key << "x" << std::setprecision(17) << duration;
        }
      }
    } else {
      // A wholly zero-duration trajectory retains only its compact state order.
      for (int state_index = 0;
           state_index < static_cast<int>(zero_state_sequence.size());
           ++state_index) {
        // Delimit every state after the first fallback component.
        if (state_index > 0) {
          key << "|";
        }
        key << zero_state_sequence[state_index] << ":" <<
          zero_path_sequence[state_index] << ":0";
      }
      // An empty zero trajectory receives a stable sentinel key.
      if (zero_state_sequence.empty()) {
        key << "|";
      }
    }

    std::map<std::string, int>::const_iterator known_group =
      group_by_key.find(key.str());
    // The first occurrence of a trajectory creates one deterministic group id.
    if (known_group == group_by_key.end()) {
      group_by_key[key.str()] = next_group_id;
      group_by_tip[tip_id] = next_group_id++;
    } else {
      // Equal committed state-duration trajectories reuse the existing group.
      group_by_tip[tip_id] = known_group->second;
    }
  }
  return group_by_tip;
}

/**
 * Resolve or create the canonical child path for one destination state.
 *
 * A state continuation returns `current_path_id`. A previously registered
 * `(parent path, state)` pair reuses its id; otherwise all path-registry vectors
 * and the transition cache grow atomically and the new id is returned. Inputs
 * use the internal positive id domain; invalid indexing fails before the batch
 * can be committed.
 */
int v34_resolve_next_path_id(int current_path_id,
                             int next_state_id,
                             std::vector<std::vector<int> >& transition_to_path,
                             std::vector<std::vector<int> >& path_components,
                             std::vector<int>& path_parent,
                             std::vector<int>& path_added_state_id,
                             int n_state) {
  const std::vector<int>& current_components = path_components[current_path_id];
  // A continuation in the path's current terminal state creates no biological
  // transition and therefore reuses the existing canonical path id.
  if (!current_components.empty() && current_components.back() == next_state_id) {
    return current_path_id;
  }
  int cached = transition_to_path[current_path_id][next_state_id];
  // A previously observed `(parent path, destination state)` transition reuses
  // its canonical destination path, preventing duplicate PUNIQ identities.
  if (cached > 0) {
    return cached;
  }
  int next_path_id = static_cast<int>(path_components.size());
  std::vector<int> next_components = current_components;
  next_components.push_back(next_state_id);
  path_components.push_back(next_components);
  path_parent.push_back(current_path_id);
  path_added_state_id.push_back(next_state_id);
  transition_to_path[current_path_id][next_state_id] = next_path_id;
  transition_to_path.push_back(std::vector<int>(n_state + 1, 0));
  return next_path_id;
}

} // namespace

/**
 * Run the C++ traversal state machine and return R-shaped active summary matrices and record lists.
 *
 * Input: a fully parsed mapped phylogeny and resolved root policy. Output: all
 * traversal-owned segment records, online TT/SS counter payloads, terminal
 * records, and diagnostic timings. The function commits only complete event
 * phases and throws before returning if frontier ownership becomes invalid.
 */
Rcpp::List pst_v34_extract_active_summary(
    const V34Input& input,
    bool recover_debug_failure,
    int debug_fail_batch_id,
    double time_tolerance) {
  double time_scale = pst_v34_batch_time_scale(time_tolerance);
  typedef std::chrono::steady_clock V34Clock;
  V34Clock::time_point setup_start = V34Clock::now();
  int max_node = 0;
  int n_state = input.state_lookup.nrows();
  // Scan all phylogeny edges once to size node-indexed topology arrays.
  for (int i = 0; i < input.nedge; ++i) {
    max_node = std::max(max_node, input.edge_parent[i]);
    max_node = std::max(max_node, input.edge_child[i]);
  }

  std::vector<std::vector<int> > child_edges(max_node + 1);
  std::vector<int> edge_for_child(max_node + 1, 0);
  std::vector<int> parent_seen(max_node + 1, 0);
  std::vector<int> child_seen(max_node + 1, 0);
  // Index each input edge by parent and child while recording root-detection
  // membership; after each iteration one edge has complete topology ownership.
  for (int edge_id = 1; edge_id <= input.nedge; ++edge_id) {
    int parent = input.edge_parent[edge_id - 1];
    int child = input.edge_child[edge_id - 1];
    child_edges[parent].push_back(edge_id);
    edge_for_child[child] = edge_id;
    parent_seen[parent] = 1;
    child_seen[child] = 1;
  }

  int root_node = input.ntips + 1;
  // Find the sole node that owns outgoing edges but has no incoming edge.
  for (int node = 1; node <= max_node; ++node) {
    // This topology signature identifies the biological root; the first match
    // is sufficient for a valid rooted phylogeny.
    if (parent_seen[node] && !child_seen[node]) {
      root_node = node;
      break;
    }
  }

  std::vector<double> node_height(max_node + 1, R_NaReal);
  node_height[root_node] = 0.0;
  std::vector<int> node_queue(1, root_node);
  // Traverse internal nodes from the root, assigning absolute endpoint times to
  // every child reached from the current parent.
  for (int pos = 0; pos < static_cast<int>(node_queue.size()); ++pos) {
    int node = node_queue[pos];
    // Advance each daughter edge by its branch length from the known parent time.
    for (int edge_id : child_edges[node]) {
      int child = input.edge_child[edge_id - 1];
      node_height[child] = node_height[node] + input.edge_length[edge_id - 1];
      // Internal children have outgoing edges and must be expanded later;
      // biological tips terminate this height traversal.
      if (child > input.ntips) {
        node_queue.push_back(child);
      }
    }
  }

  double tree_height = 0.0;
  double min_tip_height = R_PosInf;
  double max_tip_height = R_NegInf;
  // Inspect all biological tip endpoints to determine tree height and whether
  // terminal closures are ultrametric or require independent times.
  for (int tip = 1; tip <= input.ntips; ++tip) {
    min_tip_height = std::min(min_tip_height, node_height[tip]);
    max_tip_height = std::max(max_tip_height, node_height[tip]);
  }
  // Ultrametric classification must use the same integer resolution as event
  // batching. A fixed machine-epsilon threshold can call two endpoints equal
  // even when the caller-selected tolerance puts them in distinct terminal
  // batches; closing both at their midpoint then truncates the later survivor's
  // scenario support. Same-key endpoints share the rounded transactional time,
  // while distinct-key endpoints retain independent terminal closures.
  bool is_ultrametric =
    pst_v34_batch_time_key(min_tip_height, time_scale) ==
    pst_v34_batch_time_key(max_tip_height, time_scale);
  tree_height = is_ultrametric ?
    pst_v34_batch_round_time(max_tip_height, time_scale) :
    max_tip_height;

  std::vector<std::vector<int> > edge_tips(input.nedge + 1);
  std::vector<std::vector<int> > node_tips(max_node + 1);
  std::function<std::vector<int>(int)> collect_tips = [&](int node) -> std::vector<int> {
    // A cached descendant set avoids recomputing the same subtree for each
    // ancestral edge that requests its biological membership.
    if (!node_tips[node].empty()) {
      return node_tips[node];
    }
    // A biological tip is its own singleton descendant set and ends recursion.
    if (node <= input.ntips) {
      node_tips[node].push_back(node);
      return node_tips[node];
    }
    std::vector<int> tips;
    // Union descendant-tip sets from every daughter edge of this internal node.
    for (int edge_id : child_edges[node]) {
      std::vector<int> child_tips = collect_tips(input.edge_child[edge_id - 1]);
      tips.insert(tips.end(), child_tips.begin(), child_tips.end());
    }
    std::sort(tips.begin(), tips.end());
    node_tips[node] = tips;
    return node_tips[node];
  };
  collect_tips(root_node);
  // Assign each phylogeny edge the cached descendant-tip set of its child node.
  for (int edge_id = 1; edge_id <= input.nedge; ++edge_id) {
    edge_tips[edge_id] = collect_tips(input.edge_child[edge_id - 1]);
  }

  std::vector<std::vector<int> > edge_state_ids(input.nedge + 1);
  std::vector<std::vector<double> > map_end_offsets(input.nedge + 1);
  std::vector<int> map_lengths(input.nedge + 1, 0);
  // Expand the compact parsed simmap arrays into edge-local state sequences and
  // cumulative event offsets used by the traversal scheduler.
  for (int edge_id = 1; edge_id <= input.nedge; ++edge_id) {
    int start = input.map_edge_offset[edge_id - 1];
    int end = input.map_edge_offset[edge_id] - 1;
    double offset = 0.0;
    // Append each mapped segment in biological order and accumulate its endpoint
    // offset from the parent node.
    for (int segment_pos = start; segment_pos <= end; ++segment_pos) {
      int idx = segment_pos - 1;
      edge_state_ids[edge_id].push_back(input.map_state_id[idx]);
      offset += input.map_duration[idx];
      map_end_offsets[edge_id].push_back(offset);
    }
    map_lengths[edge_id] = static_cast<int>(edge_state_ids[edge_id].size());
  }

  auto natural_event_phase = [&](int edge_id, int step) {
    // A remaining mapped segment is an incoming-edge state transition. Once
    // segments are exhausted, the child topology selects cladogenesis or a
    // terminal closure.
    if (map_lengths[edge_id] > step) {
      return static_cast<int>(V34_BATCH_INCOMING_ANAGENETIC);
    }
    int child_node_id = input.edge_child[edge_id - 1];
    // Internal children split before terminal closures at the same time key.
    if (child_node_id > input.ntips) {
      return static_cast<int>(V34_BATCH_CLADOGENESIS);
    }
    return static_cast<int>(V34_BATCH_TERMINAL);
  };

  std::vector<V34ScenarioEdgeRecord> records(1);
  int next_edge_id = 1;
  int next_node_id = 2;
  std::vector<int> current_state_id(input.ntips + 1, NA_INTEGER);
  std::vector<int> current_edge_node_id(input.ntips + 1, NA_INTEGER);
  std::vector<int> current_scenario_edge_id(input.ntips + 1, NA_INTEGER);
  std::vector<int> current_phylo_edge_id(input.ntips + 1, NA_INTEGER);
  std::vector<int> current_path_id(input.ntips + 1, NA_INTEGER);

  std::vector<std::vector<int> > path_components(1);
  std::vector<int> root_components;
  // Build the known biological root path from supplied anchor-state ids.
  for (int id : input.root_anchor_state_ids) {
    // Missing synthetic placeholders are not biological path components;
    // observed anchor ids are retained in order.
    if (id != NA_INTEGER) {
      root_components.push_back(id);
    }
  }
  path_components.push_back(root_components);
  std::vector<int> path_parent(2, NA_INTEGER);
  std::vector<int> path_added_state_id(2, NA_INTEGER);
  std::vector<std::vector<int> > transition_to_path(2, std::vector<int>(n_state + 1, 0));

  // Initialize the traversal-owned TT counter before any biological lineage is
  // placed on the frontier. The root row captured below is therefore the first
  // committed state seen by every online counter family.
  int root_anchor_state_id = input.root_anchor_state_ids.size() > 0 ?
    input.root_anchor_state_ids[0] : NA_INTEGER;
  V34LiveTTState online_counters(
    input.state_lookup,
    input.ntips,
    root_anchor_state_id,
    input.root_is_synthetic,
    time_scale
  );
  // Register the root path before any initial frontier can expose it. NA parent
  // and added-state values identify a traversal initial condition.
  online_counters.register_path(1, NA_INTEGER, NA_INTEGER);

  std::vector<int> terminal_tip_ids;
  std::vector<int> terminal_phylo_edge_ids;
  std::vector<int> terminal_scenario_edge_ids;
  std::vector<int> terminal_path_ids;
  std::vector<double> terminal_times;

  std::vector<std::vector<int> > state_slices;
  std::vector<std::vector<int> > edge_node_slices;
  std::vector<std::vector<int> > phylo_edge_slices;
  std::vector<std::vector<int> > scenario_edge_slices;
  std::vector<std::vector<int> > path_slices;
  std::vector<double> time_vec;
  std::vector<int> event_vec;
  std::vector<std::vector<int> > event_tips;
  std::vector<V34DebugStepRecord> debug_step_records;
  int next_debug_event_id = 1;
  std::vector<V34BatchEvent> queue;
  double batch_capture_seconds = 0.0;
  // These diagnostic-only accumulators split the V34 scenario-transaction
  // overhead out of the encompassing traversal event loop. They do not alter
  // event order or ownership; each timer brackets work that already existed.
  double scenario_batch_frontier_seconds = 0.0;
  double scenario_batch_history_copy_seconds = 0.0;
  double scenario_batch_event_phase_seconds = 0.0;
  double scenario_batch_topology_commit_seconds = 0.0;
  double scenario_batch_history_restore_seconds = 0.0;
  double scenario_batch_commit_core_seconds = 0.0;
  double scenario_batch_classification_seconds = 0.0;
  std::size_t estimated_slice_count =
    static_cast<std::size_t>(input.map_state_id.size() + input.nedge + 1);
  // Reserve the traversal's conservative event upper bound once. Each mapped
  // segment boundary, topology edge, and root row can contribute at most one
  // stored slice before same-time collapse.
  state_slices.reserve(estimated_slice_count);
  edge_node_slices.reserve(estimated_slice_count);
  phylo_edge_slices.reserve(estimated_slice_count);
  scenario_edge_slices.reserve(estimated_slice_count);
  path_slices.reserve(estimated_slice_count);
  time_vec.reserve(estimated_slice_count);
  event_vec.reserve(estimated_slice_count);
  event_tips.reserve(estimated_slice_count);
  // Snapshot mode keeps one compact diagnostic record per committed row.
  // Normal traversal leaves this vector unreserved and empty.
  if (recover_debug_failure) {
    debug_step_records.reserve(estimated_slice_count);
  }
  queue.reserve(static_cast<std::size_t>(input.nedge));

  auto append_public_slice = [&](
      double event_time,
      int event_code,
      const std::vector<int>& tips,
      const V34DebugStepRecord* debug_step) {
    std::vector<int> state(input.ntips);
    std::vector<int> edge_node(input.ntips);
    std::vector<int> phylo_edge(input.ntips);
    std::vector<int> scenario_edge(input.ntips);
    std::vector<int> path(input.ntips);
    // Copy every tip's committed frontier ownership into one immutable output
    // slice; all five arrays remain index-aligned after each iteration.
    for (int tip = 1; tip <= input.ntips; ++tip) {
      state[tip - 1] = current_state_id[tip];
      edge_node[tip - 1] = current_edge_node_id[tip];
      phylo_edge[tip - 1] = current_phylo_edge_id[tip];
      scenario_edge[tip - 1] = current_scenario_edge_id[tip];
      path[tip - 1] = current_path_id[tip];
    }
    bool replaces_same_key =
      !time_vec.empty() &&
      pst_v34_batch_time_key(time_vec.back(), time_scale) ==
        pst_v34_batch_time_key(event_time, time_scale);
    // A zero-offset biological batch can share the root initialization key.
    // Replace that provisional public row so each integer time key remains one
    // committed output slice without retaining duplicate same-time rows.
    if (replaces_same_key) {
      state_slices.back() = std::move(state);
      edge_node_slices.back() = std::move(edge_node);
      phylo_edge_slices.back() = std::move(phylo_edge);
      scenario_edge_slices.back() = std::move(scenario_edge);
      path_slices.back() = std::move(path);
      time_vec.back() = event_time;
      event_vec.back() = event_code;
      event_tips.back() = tips;
      // A same-key biological transaction replaces its provisional row and its
      // diagnostic identity atomically. Normal traversal stores no debug copy.
      if (recover_debug_failure) {
        // Snapshot mode must supply metadata for every externally visible row.
        if (debug_step == nullptr || debug_step_records.empty()) {
          Rcpp::stop("V34 debug slice replacement omitted batch metadata");
        }
        debug_step_records.back() = *debug_step;
      }
    } else {
      // A later integer key owns a new public slice. Moving the completed
      // buffers avoids a second full-tip copy for the committed batch.
      state_slices.push_back(std::move(state));
      edge_node_slices.push_back(std::move(edge_node));
      phylo_edge_slices.push_back(std::move(phylo_edge));
      scenario_edge_slices.push_back(std::move(scenario_edge));
      path_slices.push_back(std::move(path));
      time_vec.push_back(event_time);
      event_vec.push_back(event_code);
      event_tips.push_back(tips);
      // Append compact diagnostics only on the explicitly requested snapshot
      // path; ordinary traversal performs no event-metadata allocation.
      if (recover_debug_failure) {
        // Snapshot mode requires one diagnostic record per committed row.
        if (debug_step == nullptr) {
          Rcpp::stop("V34 debug slice append omitted batch metadata");
        }
        debug_step_records.push_back(*debug_step);
      }
    }
  };

  auto capture_counter_phase = [&](
      double event_time,
      int event_code) {
    // Intermediate causal phases update only the private online transaction.
    // The public output buffers are written after the whole time batch validates.
    V34Clock::time_point batch_capture_start = V34Clock::now();
    online_counters.capture_committed_row(
      event_time,
      event_code,
      current_state_id,
      current_edge_node_id,
      current_phylo_edge_id,
      current_scenario_edge_id,
      current_path_id
    );
    batch_capture_seconds += std::chrono::duration<double>(
      V34Clock::now() - batch_capture_start
    ).count();
  };

  auto capture_and_append_public_slice = [&](
      double event_time,
      int event_code,
      const std::vector<int>& tips,
      const V34DebugStepRecord* debug_step,
      bool is_root_initialization) {
    // Root initialization and standalone terminal extension rows have no
    // additional same-time causes, so one validated capture can be serialized.
    if (is_root_initialization) {
      // The literal live-TT root operation owns the sole initial capture and
      // delegates to the same validated frontier accounting as later batches.
      V34Clock::time_point root_capture_start = V34Clock::now();
      online_counters.initialize_root(
        event_time,
        event_code,
        current_state_id,
        current_edge_node_id,
        current_phylo_edge_id,
        current_scenario_edge_id,
        current_path_id
      );
      batch_capture_seconds += std::chrono::duration<double>(
        V34Clock::now() - root_capture_start
      ).count();
    } else {
      // A non-root closure retains the existing ordinary frontier-capture path.
      capture_counter_phase(event_time, event_code);
    }
    append_public_slice(event_time, event_code, tips, debug_step);
  };

  V34Clock::time_point traversal_start = V34Clock::now();
  std::map<int, std::vector<int> > root_state_groups;
  // Group root-descending phylogeny edges by their first observed state so each
  // distinct initial state/path receives one scenario group.
  for (int edge_id : child_edges[root_node]) {
    // A mapped root edge must expose an initial state; only such edges can enter
    // a biological root scenario group.
    if (!edge_state_ids[edge_id].empty()) {
      root_state_groups[edge_state_ids[edge_id][0]].push_back(edge_id);
    }
  }
  // Initialize one canonical path and scenario edge for each root-state group.
  for (auto const& group : root_state_groups) {
    int state_id = group.first;
    int root_path_id = v34_resolve_next_path_id(
      1, state_id, transition_to_path, path_components, path_parent, path_added_state_id, n_state
    );
    // Path resolution may create or reuse the canonical root daughter; the
    // idempotent online registry validates both cases before batch capture.
    online_counters.register_path(
      root_path_id,
      path_parent[root_path_id],
      path_added_state_id[root_path_id]
    );
    int scenario_edge_id = v34_new_scenario_edge(
      records, next_edge_id, next_node_id, state_id,
      static_cast<int>(group.second.size()), 1, 0.0, root_path_id,
      0
    );
    // Attach every root phylogeny edge in this group to the shared initial
    // scenario and schedule its first mapped boundary.
    for (int edge_id : group.second) {
      const std::vector<int>& tips = edge_tips[edge_id];
      // Assign all descendants of one active root edge the same state, path,
      // phylogeny edge, and scenario ownership.
      for (int tip : tips) {
        current_state_id[tip] = state_id;
        current_edge_node_id[tip] = input.edge_child[edge_id - 1];
        current_scenario_edge_id[tip] = scenario_edge_id;
        current_phylo_edge_id[tip] = edge_id;
        current_path_id[tip] = root_path_id;
      }
      queue.push_back(V34BatchEvent{
        edge_id,
        1,
        scenario_edge_id,
        scenario_edge_id,
        map_end_offsets[edge_id][0],
        pst_v34_batch_min_tip(tips),
        natural_event_phase(edge_id, 1)
      });
      // A direct departure from a known biological root is a real transition.
      // Synthetic-root daughters are initial conditions and are suppressed by
      // the explicit policy branch here rather than repaired after traversal.
      if (!input.root_is_synthetic &&
          root_anchor_state_id != NA_INTEGER &&
          state_id != root_anchor_state_id) {
        online_counters.record_ctt_event(
          0.0,
          root_anchor_state_id,
          state_id,
          edge_id,
          1,
          root_path_id,
          NA_INTEGER,
          scenario_edge_id,
          input.edge_child[edge_id - 1] <= input.ntips
        );
      }
    }
  }

  std::vector<int> all_tips(input.ntips);
  std::iota(all_tips.begin(), all_tips.end(), 1);
  V34DebugStepRecord root_debug_step;
  // Snapshot mode records root initialization as a diagnostic initial
  // condition, not as a scored biological event. Root phylogeny/scenario ids
  // still identify the initialized frontier for failure localization.
  if (recover_debug_failure) {
    root_debug_step.event_types.push_back("root_initialized");
    root_debug_step.causal_phases.push_back("root_initialization");
    root_debug_step.phylogeny_edges = child_edges[root_node];
    // Tip loop:
    //   Collect each initialized source scenario once from the committed root
    //   frontier. No event id is created because initialization is pre-event.
    for (int tip = 1; tip <= input.ntips; ++tip) {
      v34_debug_append_unique(
        root_debug_step.scenario_ids,
        current_scenario_edge_id[tip]
      );
    }
  }
  capture_and_append_public_slice(
    0.0,
    3,
    all_tips,
    recover_debug_failure ? &root_debug_step : nullptr,
    true
  );

  bool debug_failed = false;
  std::string debug_failure_message;
  int debug_attempted_batch_id = NA_INTEGER;
  double debug_attempted_time = NA_REAL;
  int batch_id = 0;

  // Consume one complete integer-time batch per outer iteration. The inner loop
  // applies every causal event at that key to private working state and validates
  // each phase; public buffers receive only the final post-batch frontier.
  while (!queue.empty()) {
    ++batch_id;
    int attempted_event_index = pst_v34_batch_next_event_index(
      queue, time_scale
    );
    debug_attempted_batch_id = batch_id;
    debug_attempted_time = queue[attempted_event_index].time;

    std::unique_ptr<V34DebugBatchCheckpoint> checkpoint;
    // Snapshot recovery pays for a full pre-batch copy; normal traversal keeps
    // its allocation profile unchanged and continues to throw on failure.
    if (recover_debug_failure) {
      checkpoint.reset(new V34DebugBatchCheckpoint(
        records,
        next_edge_id,
        next_node_id,
        current_state_id,
        current_edge_node_id,
        current_scenario_edge_id,
        current_phylo_edge_id,
        current_path_id,
        path_components,
        path_parent,
        path_added_state_id,
        transition_to_path,
        online_counters,
        terminal_tip_ids,
        terminal_phylo_edge_ids,
        terminal_scenario_edge_ids,
        terminal_path_ids,
        terminal_times,
        queue
      ));
    }

    try {
    V34Clock::time_point scenario_frontier_start = V34Clock::now();
    int first_event_index = pst_v34_batch_next_event_index(queue, time_scale);
    V34BatchTransaction batch_transaction = online_counters.begin_batch(
      queue[first_event_index].time
    );
    long long batch_time_key = batch_transaction.time_key;
    double batch_time = batch_transaction.time;
    V34ScenarioBatchIntent scenario_batch_intent =
      pst_v34_scenario_batch_begin(
        batch_time_key,
        batch_time,
        current_scenario_edge_id,
        current_state_id,
        current_path_id,
        current_phylo_edge_id,
        terminal_tip_ids
      );
    scenario_batch_frontier_seconds += std::chrono::duration<double>(
      V34Clock::now() - scenario_frontier_start
    ).count();

    // V34 currently preserves the complete pre-batch scenario history because
    // phase-local topology is provisional until the batch's final canonical
    // commit. Timing this deep copy makes any history-size regression visible.
    V34Clock::time_point scenario_history_copy_start = V34Clock::now();
    std::vector<V34ScenarioEdgeRecord> scenario_records_before_batch = records;
    int next_edge_id_before_batch = next_edge_id;
    int next_node_id_before_batch = next_node_id;
    std::vector<int> scenario_owner_before_batch =
      current_scenario_edge_id;
    scenario_batch_history_copy_seconds += std::chrono::duration<double>(
      V34Clock::now() - scenario_history_copy_start
    ).count();
    V34Clock::time_point scenario_event_phase_start = V34Clock::now();
    int batch_event_code = NA_INTEGER;
    std::vector<int> batch_event_tips;
    V34DebugStepRecord batch_debug_step;
    std::map<int, int> terminal_removals_by_scenario_edge;

    std::vector<V34BatchEvent> path_reservation_events;
    // Collect events already scheduled at this key without mutating queue or
    // frontier ownership. Path IDs are identity reservations only and create no
    // destination support, scenario, event score, or public path appearance.
    for (const V34BatchEvent& queued_event : queue) {
      // Only events at this transaction's stable time key participate in the
      // canonical identity reservation pass.
      if (pst_v34_batch_time_key(queued_event.time, time_scale) == batch_time_key) {
        path_reservation_events.push_back(queued_event);
      }
    }
    std::sort(
      path_reservation_events.begin(),
      path_reservation_events.end(),
      [](const V34BatchEvent& left, const V34BatchEvent& right) {
        // Preserve the established canonical ID order: descendant tip first,
        // then phylogeny edge, independent of causal execution phase.
        if (left.priority_tip != right.priority_tip) {
          return left.priority_tip < right.priority_tip;
        }
        return left.edge_id < right.edge_id;
      }
    );
    // Reserve all directly discoverable destination paths in canonical order.
    // Later phase handlers reuse these ids while applying biological deltas.
    for (const V34BatchEvent& planned_event : path_reservation_events) {
      int planned_edge_id = planned_event.edge_id;
      int planned_step = planned_event.step;
      const std::vector<int>& planned_tips = edge_tips[planned_edge_id];
      // Every scheduled event owns at least one descendant tip by queue contract.
      if (planned_tips.empty()) {
        Rcpp::stop("V34 path reservation event has no descendants");
      }
      int source_path_id = current_path_id[planned_tips[0]];
      // An incoming anagenetic boundary reserves its next state/path pair.
      if (map_lengths[planned_edge_id] > planned_step) {
        int destination_state_id =
          edge_state_ids[planned_edge_id][planned_step];
        int destination_path_id = v34_resolve_next_path_id(
          source_path_id,
          destination_state_id,
          transition_to_path,
          path_components,
          path_parent,
          path_added_state_id,
          n_state
        );
        online_counters.register_path(
          destination_path_id,
          path_parent[destination_path_id],
          path_added_state_id[destination_path_id]
        );
        continue;
      }
      int planned_child_node = input.edge_child[planned_edge_id - 1];
      // Terminal events create no destination path; internal children reserve
      // every changed daughter's initial state/path in daughter edge order.
      if (planned_child_node <= input.ntips) {
        continue;
      }
      // Reserve each internal child's daughter path in deterministic phylogeny
      // edge order before causal phases mutate the frontier.
      for (int daughter_edge_id : child_edges[planned_child_node]) {
        int daughter_state_id = edge_state_ids[daughter_edge_id][0];
        int daughter_path_id = v34_resolve_next_path_id(
          source_path_id,
          daughter_state_id,
          transition_to_path,
          path_components,
          path_parent,
          path_added_state_id,
          n_state
        );
        online_counters.register_path(
          daughter_path_id,
          path_parent[daughter_path_id],
          path_added_state_id[daughter_path_id]
        );
      }
    }

    // Repeatedly select the next deterministic event while it belongs to this
    // same-time transaction. Newly queued zero-offset daughter events join the
    // current batch because selection is repeated after every causal phase.
    while (!queue.empty()) {
      int event_index = pst_v34_batch_next_event_index(queue, time_scale);
      // A different integer key belongs to the next transaction and must not
      // inspect or mutate the current batch's virtual post-state.
      if (pst_v34_batch_time_key(queue[event_index].time, time_scale) != batch_time_key) {
        break;
      }
      V34BatchEvent event = queue[event_index];
      queue.erase(queue.begin() + event_index);

      int edge_id = event.edge_id;
      int step = event.step;
      int scenario_edge_id = event.scenario_edge_id;
      int batch_origin_scenario_edge_id = event.batch_origin_scenario_edge_id;
      double event_time = event.time;
      int child_node = input.edge_child[edge_id - 1];
      const std::vector<int>& descendant_tips = edge_tips[edge_id];
      int frozen_source_scenario_id =
        pst_v34_scenario_batch_tip(
          scenario_batch_intent,
          descendant_tips[0]
        ).source_scenario_id;
      pst_v34_scenario_batch_mark_source_affected(
        scenario_batch_intent,
        frozen_source_scenario_id
      );

      // Freeze and validate this event's source withdrawal before any handler
      // creates paths, scenarios, daughter ownership, or terminal records.
      if (batch_origin_scenario_edge_id < 1 ||
          batch_origin_scenario_edge_id >= static_cast<int>(records.size()) ||
          records[batch_origin_scenario_edge_id].sizes.empty()) {
        Rcpp::stop("V34 batch event has no valid source scenario");
      }
      std::map<int, int>::const_iterator frozen_capacity =
        batch_transaction.source_capacity.find(batch_origin_scenario_edge_id);
      // A scenario that continues in place may append a new size segment before
      // another same-time event consumes its old source. Once encountered, the
      // transaction's pre-mutation capacity remains authoritative.
      int source_capacity = frozen_capacity ==
          batch_transaction.source_capacity.end() ?
        records[batch_origin_scenario_edge_id].sizes.back() :
        frozen_capacity->second;
      pst_v34_batch_plan_source_removal(
        batch_transaction,
        batch_origin_scenario_edge_id,
        source_capacity,
        event.phase
      );
      // Snapshot diagnostics identify every queued biological boundary before
      // its handler mutates private working state. These compact vectors are
      // discarded with the transaction if validation or a handler fails.
      if (recover_debug_failure) {
        batch_debug_step.event_ids.push_back(next_debug_event_id++);
        batch_debug_step.phylogeny_edges.push_back(edge_id);
        batch_debug_step.causal_phases.push_back(
          v34_debug_phase_label(event.phase)
        );
        v34_debug_append_unique(
          batch_debug_step.scenario_ids,
          batch_origin_scenario_edge_id
        );
      }

      // A remaining mapped segment means this boundary is an anagenetic state
      // transition on the incoming phylogeny edge rather than its child node.
      if (map_lengths[edge_id] > step) {
        // Anagenetic events retain one explicit type even when independent or
        // post-cladogenetic changes share this same committed time key.
        if (recover_debug_failure) {
          batch_debug_step.event_types.push_back("anagenetic");
        }
        V34ScenarioEdgeRecord& active = records[scenario_edge_id];
        int previous_state_id = active.state_ids.back();
        int previous_size = active.sizes.back();
        int previous_path_id = current_path_id[descendant_tips[0]];
        int changed_state_id = edge_state_ids[edge_id][step];
        int changed_path_id = v34_resolve_next_path_id(
          previous_path_id, changed_state_id, transition_to_path,
          path_components, path_parent, path_added_state_id, n_state
        );
        // Register the anagenetic destination before its same-time batch can
        // expose or score the path.
        online_counters.register_path(
          changed_path_id,
          path_parent[changed_path_id],
          path_added_state_id[changed_path_id]
        );

        // A new destination path records separate gross loss and gain entries.
        if (previous_path_id != changed_path_id) {
          online_counters.record_lineage_loss(event_time, previous_path_id);
          online_counters.record_lineage_gain(event_time, changed_path_id);
        }

        // Move the affected phylogeny lineage's full descendant membership to its
        // destination state/path while preserving the same phylogeny edge.
        for (int tip : descendant_tips) {
          current_state_id[tip] = changed_state_id;
          current_edge_node_id[tip] = child_node;
          current_phylo_edge_id[tip] = edge_id;
          current_path_id[tip] = changed_path_id;
        }

        // A shared source scenario must split the changing singleton from the
        // unrelated lineages that retain the source state and path.
        if (previous_size > 1) {
          v34_close_active_segment(records, scenario_edge_id, event_time);
          int parent_node = records[scenario_edge_id].child_nodes.back();
          int conserved_size = previous_size - 1;
          int conserved_edge_id = NA_INTEGER;
          // Positive source support creates a conserved child scenario and moves
          // pending unrelated events to that replacement edge.
          if (conserved_size > 0) {
            conserved_edge_id = v34_new_scenario_edge(
              records, next_edge_id, next_node_id, previous_state_id,
              conserved_size, parent_node, event_time, previous_path_id,
              scenario_edge_id
            );
            pst_v34_batch_reassign_scenario(
              queue,
              scenario_edge_id,
              conserved_edge_id
            );
            std::vector<char> is_desc(input.ntips + 1, 0);
            // Mark descendants of the changing phylogeny edge so they are not
            // reassigned to the conserved source child.
            for (int tip : descendant_tips) is_desc[tip] = 1;
            // Reassign every unrelated tip still owned by the replaced source
            // scenario to its conserved child.
            for (int tip = 1; tip <= input.ntips; ++tip) {
              // Only source-scenario members outside the changing lineage retain
              // source ownership after the split.
              if (current_scenario_edge_id[tip] == scenario_edge_id && !is_desc[tip]) {
                current_scenario_edge_id[tip] = conserved_edge_id;
              }
            }
          }
          int changed_edge_id = v34_new_scenario_edge(
            records, next_edge_id, next_node_id, changed_state_id,
            1, parent_node, event_time, changed_path_id,
            scenario_edge_id
          );
          // Assign all descendants of the changed phylogeny edge to the new
          // singleton destination scenario.
          for (int tip : descendant_tips) {
            current_scenario_edge_id[tip] = changed_edge_id;
          }
          double next_event_time = node_height[input.edge_parent[edge_id - 1]] + map_end_offsets[edge_id][step];
          queue.push_back(V34BatchEvent{
            edge_id,
            step + 1,
            changed_edge_id,
            changed_edge_id,
            next_event_time,
            pst_v34_batch_min_tip(descendant_tips),
            natural_event_phase(edge_id, step + 1)
          });
        } else {
          // A singleton source scenario changes in place: its segment closes and a
          // destination-state segment begins on the same scenario edge.
          v34_close_active_segment(records, scenario_edge_id, event_time);
          // Scenario loss/gain and Stay/Leave identity are resolved only from
          // the complete post-batch semantic groups. This provisional in-place
          // mutation therefore changes traversal ownership but scores no TT.
          v34_append_segment(records, scenario_edge_id, changed_state_id, previous_size, event_time, changed_path_id);
          // Descendants retain the in-place scenario edge after its state/path
          // segment is updated.
          for (int tip : descendant_tips) {
            current_scenario_edge_id[tip] = scenario_edge_id;
          }
          double next_event_time = node_height[input.edge_parent[edge_id - 1]] + map_end_offsets[edge_id][step];
          queue.push_back(V34BatchEvent{
            edge_id,
            step + 1,
            scenario_edge_id,
            scenario_edge_id,
            next_event_time,
            pst_v34_batch_min_tip(descendant_tips),
            natural_event_phase(edge_id, step + 1)
          });
        }
        // Record the biological edge change only after scenario splitting has
        // resolved the destination scenario. CTT retains singleton in-place
        // changes, while LDIF later applies the stricter scenario-change rule.
        scenario_batch_intent.transitions.push_back(
          V34PendingScenarioTransition{
            event_time,
            previous_state_id,
            changed_state_id,
            edge_id,
            previous_path_id,
            changed_path_id,
            frozen_source_scenario_id,
            descendant_tips[0],
            child_node <= input.ntips
          }
        );

        batch_event_code = 1;
        batch_event_tips = descendant_tips;
        continue;
      }

      // Exhausting the incoming edge at an internal child begins cladogenesis;
      // terminal children are handled by the closure branch below.
      if (child_node > input.ntips) {
        const std::vector<int>& daughter_edges = child_edges[child_node];
        std::vector<int> daughter_descendants;
        std::vector<int> daughter_state_ids;
        std::vector<int> daughter_path_ids;
        int previous_path_id = current_path_id[descendant_tips[0]];
        // Resolve every daughter's initial state and canonical path before any
        // child scenario groups are allocated.
        for (int daughter_edge_id : daughter_edges) {
          const std::vector<int>& tips = edge_tips[daughter_edge_id];
          daughter_descendants.insert(daughter_descendants.end(), tips.begin(), tips.end());
          int daughter_state_id = edge_state_ids[daughter_edge_id][0];
          int daughter_path_id = current_path_id[tips[0]];
          // A daughter whose initial state differs from its parent receives a new
          // destination path; conserved daughters retain the source path.
          if (daughter_state_id != current_state_id[tips[0]]) {
            daughter_path_id = v34_resolve_next_path_id(
              daughter_path_id, daughter_state_id, transition_to_path,
              path_components, path_parent, path_added_state_id, n_state
            );
            // Register each changed daughter's canonical destination before the
            // cladogenetic transaction commits all daughters atomically.
            online_counters.register_path(
              daughter_path_id,
              path_parent[daughter_path_id],
              path_added_state_id[daughter_path_id]
            );
          }
          daughter_state_ids.push_back(daughter_state_id);
          daughter_path_ids.push_back(daughter_path_id);
          // Move every descendant of this daughter onto its outgoing phylogeny
          // edge and resolved state/path ownership.
          for (int tip : tips) {
            current_state_id[tip] = daughter_state_id;
            current_edge_node_id[tip] = input.edge_child[daughter_edge_id - 1];
            current_phylo_edge_id[tip] = daughter_edge_id;
            current_path_id[tip] = daughter_path_id;
          }
        }

        // Aggregate the signed path change within this one cladogenetic event.
        // Thus P -> P,P forms one net new lineage on P, while P -> Q,R loses
        // one P and forms one Q plus one R. Each simultaneous biological event
        // is recorded separately, so a loss in another event cannot cancel a
        // formation here before cumulative accounting.
        std::map<int, int> lineage_event_delta;
        lineage_event_delta[previous_path_id] -= 1;
        for (int daughter_path_id : daughter_path_ids) {
          lineage_event_delta[daughter_path_id] += 1;
        }
        for (const auto& path_delta : lineage_event_delta) {
          if (path_delta.second < 0) {
            online_counters.record_lineage_loss(
              event_time,
              path_delta.first,
              -path_delta.second
            );
          } else if (path_delta.second > 0) {
            online_counters.record_lineage_gain(
              event_time,
              path_delta.first,
              path_delta.second
            );
          }
        }

        V34ScenarioEdgeRecord& active = records[scenario_edge_id];
        int previous_state_id = active.state_ids.back();
        int previous_size = active.sizes.back();
        v34_close_active_segment(records, scenario_edge_id, event_time);
        int parent_node = records[scenario_edge_id].child_nodes.back();

        typedef std::pair<int, int> V34StatePathKey;
        std::map<V34StatePathKey, int> group_counts;
        // Group each daughter by both state and complete transition path. Equal
        // states with different histories are distinct biological scenarios and
        // must never share ownership or size support.
        for (int daughter_id = 0;
             daughter_id < static_cast<int>(daughter_state_ids.size());
             ++daughter_id) {
          V34StatePathKey daughter_key(
            daughter_state_ids[daughter_id],
            daughter_path_ids[daughter_id]
          );
          ++group_counts[daughter_key];
        }
        // Unaffected members of the source scenario retain both its state and
        // path after the splitting parent lineage is removed.
        if (previous_size > 1) {
          group_counts[V34StatePathKey(previous_state_id, previous_path_id)] +=
            previous_size - 1;
        }

        std::map<V34StatePathKey, int> child_edges_by_group;
        // One destination group means the scenario did not differentiate. Its
        // existing edge receives a new path segment even when cladogenesis
        // increased lineage multiplicity; this is a CSTT move, never LDIF.
        if (group_counts.size() == 1) {
          const std::pair<const V34StatePathKey, int>& item =
            *group_counts.begin();
          if (item.second <= 0) {
            Rcpp::stop("V34 cladogenesis produced an empty scenario group");
          }
          int state_value = item.first.first;
          int path_value = item.first.second;
          // A one-group cladogenetic result is a candidate scenario step. Its
          // gross scenario accounting is deferred until the whole same-time
          // batch has established the source scenario's complete child set.
          v34_append_segment(
            records,
            scenario_edge_id,
            state_value,
            item.second,
            event_time,
            path_value
          );
          child_edges_by_group[item.first] = scenario_edge_id;
        } else {
          // Multiple provisional `(state, path)` groups allocate traversal
          // children only. Batch commit classifies the resolved semantic child
          // groups after every simultaneous event has been applied.
          for (const std::pair<const V34StatePathKey, int>& item : group_counts) {
            // Empty groups can arise only from a malformed delta and are ignored
            // before scenario topology is allocated.
            if (item.second <= 0) {
              continue;
            }
            int state_value = item.first.first;
            int path_value = item.first.second;
            int child_edge_id = v34_new_scenario_edge(
              records, next_edge_id, next_node_id, state_value,
              item.second, parent_node, event_time, path_value,
              scenario_edge_id
            );
            child_edges_by_group[item.first] = child_edge_id;
          }
        }

        // A source scenario with unrelated support must reassign those lineages
        // and their queued events to the conserved `(state, path)` child.
        if (previous_size > 1) {
          V34StatePathKey source_key(previous_state_id, previous_path_id);
          int conserved_edge_id = child_edges_by_group[source_key];
          pst_v34_batch_reassign_scenario(
            queue,
            scenario_edge_id,
            conserved_edge_id
          );
          std::vector<char> is_daughter(input.ntips + 1, 0);
          // Mark all descendants of the splitting parent so they are excluded
          // while unrelated source-scenario members are reassigned.
          for (int tip : daughter_descendants) {
            is_daughter[tip] = 1;
          }
          // Reassign every unrelated descendant tip still owned by the source
          // scenario; after each iteration ownership remains one-to-one.
          for (int tip = 1; tip <= input.ntips; ++tip) {
            // Only live members of the replaced source scenario move to its
            // conserved child; daughter groups are assigned below.
            if (current_scenario_edge_id[tip] == scenario_edge_id && !is_daughter[tip]) {
              current_scenario_edge_id[tip] = conserved_edge_id;
            }
          }
        }

        // Assign each daughter to the child scenario matching both its state and
        // path, then queue its next mapped boundary with that ownership.
        for (int daughter_id = 0;
             daughter_id < static_cast<int>(daughter_edges.size());
             ++daughter_id) {
          int daughter_edge_id = daughter_edges[daughter_id];
          int daughter_state_id = daughter_state_ids[daughter_id];
          int daughter_path_id = daughter_path_ids[daughter_id];
          V34StatePathKey daughter_key(daughter_state_id, daughter_path_id);
          int daughter_scenario_edge_id = child_edges_by_group[daughter_key];
          const std::vector<int>& tips = edge_tips[daughter_edge_id];
          // Every descendant tip of one daughter inherits the same scenario edge.
          for (int tip : tips) {
            current_scenario_edge_id[tip] = daughter_scenario_edge_id;
          }
          double daughter_event_time =
            event_time + map_end_offsets[daughter_edge_id][0];
          int daughter_event_phase = natural_event_phase(daughter_edge_id, 1);
          // A zero-duration first daughter segment whose next boundary is a real
          // transition belongs to the post-cladogenetic outgoing phase.
          if (pst_v34_batch_time_key(daughter_event_time, time_scale) == batch_time_key &&
              map_lengths[daughter_edge_id] > 1) {
            daughter_event_phase = V34_BATCH_OUTGOING_ZERO_OFFSET;
          }
          queue.push_back(V34BatchEvent{
            daughter_edge_id,
            1,
            daughter_scenario_edge_id,
            daughter_scenario_edge_id,
            daughter_event_time,
            pst_v34_batch_min_tip(tips),
            daughter_event_phase
          });
        }

        // Score changed daughter edges after all daughter scenario groups have
        // been created. Each daughter remains a distinct CTT event; LDIF and
        // PUNIQ deduplicate later at their documented path-level granularities.
        for (int daughter_id = 0;
             daughter_id < static_cast<int>(daughter_edges.size());
             ++daughter_id) {
          int daughter_state_id = daughter_state_ids[daughter_id];
          // Conserved daughters are cladogenetic lineage additions, not state changes.
          if (daughter_state_id == previous_state_id) {
            continue;
          }
          int daughter_edge_id = daughter_edges[daughter_id];
          const std::vector<int>& changed_daughter_tips =
            edge_tips[daughter_edge_id];
          int daughter_source_scenario_id =
            pst_v34_scenario_batch_tip(
              scenario_batch_intent,
              changed_daughter_tips[0]
            ).source_scenario_id;
          scenario_batch_intent.transitions.push_back(
            V34PendingScenarioTransition{
              event_time,
              previous_state_id,
              daughter_state_id,
              daughter_edge_id,
              previous_path_id,
              daughter_path_ids[daughter_id],
              daughter_source_scenario_id,
              changed_daughter_tips[0],
              input.edge_child[daughter_edge_id - 1] <= input.ntips
            }
          );
        }

        int event_code = 2;
        // Daughter states classify conserved versus state-splitting cladogenesis
        // for public event metadata.
        if (!daughter_state_ids.empty()) {
          std::vector<int> unique_states = daughter_state_ids;
          std::sort(unique_states.begin(), unique_states.end());
          unique_states.erase(std::unique(unique_states.begin(), unique_states.end()), unique_states.end());
          // Multiple daughter states identify a cladogenetic state split; a single
          // daughter state uses the conserved split event code.
          if (unique_states.size() > 1) event_code = 21;
        }
        // Preserve conserved versus state-splitting cladogenesis as the event
        // type aligned with this queued topology boundary.
        if (recover_debug_failure) {
          batch_debug_step.event_types.push_back(
            event_code == 21 ?
              "cladogenetic_state_split" : "cladogenetic"
          );
        }
        batch_event_code = event_code;
        batch_event_tips = daughter_descendants;
        continue;
      }

      if (scenario_edge_id < 1 ||
          scenario_edge_id >= static_cast<int>(records.size()) ||
          records[scenario_edge_id].sizes.empty() ||
          records[scenario_edge_id].state_ids.empty() ||
          records[scenario_edge_id].path_ids.empty()) {
        Rcpp::stop("V34 terminal event has no active scenario segment");
      }
      // Scenario support is mutated only after every same-time terminal event
      // has been collected. Applying removals one by one would invent an
      // order-dependent intermediate support run inside an atomic batch.
      ++terminal_removals_by_scenario_edge[scenario_edge_id];
      online_counters.record_lineage_loss(
        event_time,
        current_path_id[descendant_tips[0]]
      );
      // Record every biological tip ending on this terminal edge and mark its
      // exact endpoint for TT inclusion/removal semantics.
      for (int tip : descendant_tips) {
        V34ScenarioBatchTipIntent& terminal_intent =
          pst_v34_scenario_batch_tip(scenario_batch_intent, tip);
        pst_v34_scenario_batch_mark_terminal(
          scenario_batch_intent,
          tip
        );
        terminal_tip_ids.push_back(tip);
        terminal_phylo_edge_ids.push_back(edge_id);
        terminal_scenario_edge_ids.push_back(
          terminal_intent.source_scenario_id
        );
        terminal_path_ids.push_back(current_path_id[tip]);
        terminal_times.push_back(event_time);
        // Mark the terminal endpoint before capturing its row. Endpoint inclusion
        // keeps this lineage in the current TT slice and excludes it only later.
        online_counters.mark_terminal(tip, event_time);
        // Terminal tip ids identify every endpoint closed by this event while
        // suppressing duplicate descendants in malformed repeated schedules.
        if (recover_debug_failure) {
          v34_debug_append_unique(batch_debug_step.terminal_ids, tip);
        }
      }
      // Terminal closure remains visible alongside simultaneous anagenetic or
      // cladogenetic events rather than replacing their diagnostic identities.
      if (recover_debug_failure) {
        batch_debug_step.event_types.push_back("terminal");
      }
      batch_event_code = 4;
      batch_event_tips = descendant_tips;
    }

    // Apply one atomic terminal support change per canonical scenario edge.
    // A positive remainder is the complete post-batch survivor support; zero
    // closes the last pre-batch run directly with no artificial size step.
    for (const std::pair<const int, int>& terminal_removal :
         terminal_removals_by_scenario_edge) {
      int scenario_edge_id = terminal_removal.first;
      int ending_lineage_count = terminal_removal.second;
      V34ScenarioEdgeRecord& terminal_scenario = records[scenario_edge_id];
      int terminal_source_size = terminal_scenario.sizes.back();
      if (ending_lineage_count < 1 ||
          terminal_source_size < ending_lineage_count) {
        Rcpp::stop("V34 terminal batch exceeds scenario lineage support");
      }
      if (terminal_source_size > ending_lineage_count) {
        int continuing_state_id = terminal_scenario.state_ids.back();
        int continuing_path_id = terminal_scenario.path_ids.back();
        v34_close_active_segment(records, scenario_edge_id, batch_time);
        v34_append_segment(
          records,
          scenario_edge_id,
          continuing_state_id,
          terminal_source_size - ending_lineage_count,
          batch_time,
          continuing_path_id
        );
      } else {
        double terminal_segment_end = is_ultrametric ? tree_height : batch_time;
        v34_close_terminal_segment(
          records,
          scenario_edge_id,
          terminal_segment_end
        );
      }
    }

    // A selected batch must process at least one event. Serialize only its
    // validated final frontier; no intermediate causal phase enters PST output.
    if (batch_event_code == NA_INTEGER) {
      Rcpp::stop("V34 same-time batch completed without a biological event");
    }
    scenario_batch_event_phase_seconds += std::chrono::duration<double>(
      V34Clock::now() - scenario_event_phase_start
    ).count();
    // Deterministic failure injection is private to the debug contract test and
    // runs after all partial mutations but before validation/public commit.
    if (debug_fail_batch_id != NA_INTEGER &&
        batch_id == debug_fail_batch_id) {
      Rcpp::stop("V34 injected debug batch failure");
    }
    pst_v34_batch_validate(batch_transaction);
    V34Clock::time_point scenario_topology_commit_start = V34Clock::now();
    // Discard provisional phase-local scenario topology while retaining final
    // state/path/phylogeny intent, gross counter deltas, terminal records, and
    // future biological boundaries. The single commit below owns final records.
    V34Clock::time_point scenario_history_restore_start = V34Clock::now();
    records = std::move(scenario_records_before_batch);
    next_edge_id = next_edge_id_before_batch;
    next_node_id = next_node_id_before_batch;
    current_scenario_edge_id = std::move(scenario_owner_before_batch);
    scenario_batch_history_restore_seconds += std::chrono::duration<double>(
      V34Clock::now() - scenario_history_restore_start
    ).count();

    V34Clock::time_point scenario_commit_core_start = V34Clock::now();
    v34_commit_scenario_batch(
      scenario_batch_intent,
      records,
      next_edge_id,
      next_node_id,
      current_scenario_edge_id,
      current_state_id,
      current_path_id,
      current_phylo_edge_id,
      queue,
      online_counters
    );
    scenario_batch_commit_core_seconds += std::chrono::duration<double>(
      V34Clock::now() - scenario_commit_core_start
    ).count();
    // Only the fully resolved post-batch scenario topology can identify Stay
    // versus Leave. This metadata assignment occurs before the sole TT capture.
    V34Clock::time_point scenario_classification_start = V34Clock::now();
    v34_classify_pending_scenario_edges(records);
    scenario_batch_classification_seconds += std::chrono::duration<double>(
      V34Clock::now() - scenario_classification_start
    ).count();
    scenario_batch_topology_commit_seconds += std::chrono::duration<double>(
      V34Clock::now() - scenario_topology_commit_start
    ).count();
    // Every debug event id has exactly one type, phase, and phylogeny edge.
    // Misalignment is a diagnostic construction failure and cannot be emitted.
    if (recover_debug_failure &&
        (batch_debug_step.event_types.size() !=
           batch_debug_step.event_ids.size() ||
         batch_debug_step.causal_phases.size() !=
           batch_debug_step.event_ids.size() ||
         batch_debug_step.phylogeny_edges.size() !=
           batch_debug_step.event_ids.size())) {
      Rcpp::stop("V34 debug batch metadata is not event-aligned");
    }
    // The validated complete post-batch frontier owns the sole full counter
    // capture and all public invariants before any slice is appended.
    capture_counter_phase(batch_time, batch_event_code);
    append_public_slice(
      batch_time,
      batch_event_code,
      batch_event_tips,
      recover_debug_failure ? &batch_debug_step : nullptr
    );
    } catch (const std::exception& error) {
      // Normal execution exposes the original traversal failure and never
      // returns a partial object.
      if (!recover_debug_failure) {
        throw;
      }

      // Snapshot mode restores every mutation-capable structure to the exact
      // checkpoint preceding the attempted transaction. No event, path,
      // scenario, terminal, or counter effect from the failed batch survives.
      records = checkpoint->records;
      next_edge_id = checkpoint->next_edge_id;
      next_node_id = checkpoint->next_node_id;
      current_state_id = checkpoint->current_state_id;
      current_edge_node_id = checkpoint->current_edge_node_id;
      current_scenario_edge_id = checkpoint->current_scenario_edge_id;
      current_phylo_edge_id = checkpoint->current_phylo_edge_id;
      current_path_id = checkpoint->current_path_id;
      path_components = checkpoint->path_components;
      path_parent = checkpoint->path_parent;
      path_added_state_id = checkpoint->path_added_state_id;
      transition_to_path = checkpoint->transition_to_path;
      online_counters = checkpoint->online_counters;
      terminal_tip_ids = checkpoint->terminal_tip_ids;
      terminal_phylo_edge_ids = checkpoint->terminal_phylo_edge_ids;
      terminal_scenario_edge_ids = checkpoint->terminal_scenario_edge_ids;
      terminal_path_ids = checkpoint->terminal_path_ids;
      terminal_times = checkpoint->terminal_times;
      queue = checkpoint->queue;
      debug_failed = true;
      debug_failure_message = error.what();
      break;
    }
  }

  // Non-ultrametric processing can end before the global maximum height; append
  // one final closure row so integrated exposure reaches the full tree extent.
  if (!debug_failed && !time_vec.empty() && time_vec.back() < tree_height) {
    V34DebugStepRecord closure_debug_step;
    // A non-ultrametric final extension integrates remaining exposure but does
    // not invent another biological event. Its phase labels the formatting
    // closure explicitly and its batch size remains zero.
    if (recover_debug_failure) {
      closure_debug_step.event_types.push_back("terminal_closure");
      closure_debug_step.causal_phases.push_back("terminal_closure");
    }
    capture_and_append_public_slice(
      tree_height,
      4,
      all_tips,
      recover_debug_failure ? &closure_debug_step : nullptr,
      false
    );
  }
  V34Clock::time_point traversal_end = V34Clock::now();

  // The batch loop writes exactly one public frontier per integer time key, so
  // output materialization consumes the committed buffers directly.
  int nslice = static_cast<int>(time_vec.size());
  Rcpp::IntegerMatrix state_id_matrix(input.ntips, nslice);
  Rcpp::IntegerMatrix edge_node_matrix(input.ntips, nslice);
  Rcpp::IntegerMatrix phylo_edge_matrix(input.ntips, nslice);
  Rcpp::IntegerMatrix scenario_edge_matrix(input.ntips, nslice);
  Rcpp::IntegerMatrix path_id_matrix(input.ntips, nslice);
  Rcpp::NumericVector out_time(nslice);
  Rcpp::IntegerVector out_event(nslice);
  Rcpp::List out_event_tips(nslice);
  Rcpp::List out_debug_step_metadata = R_NilValue;
  Rcpp::CharacterVector state_column = input.state_lookup["state"];
  Rcpp::IntegerVector state_id_column = input.state_lookup["state_id"];
  Rcpp::CharacterVector state_labels(n_state + 1, NA_STRING);
  // Cache each R state string by integer id once. Matrix materialization then
  // reuses its CHARSXP instead of allocating the same label for every cell.
  for (int i = 0; i < state_id_column.size(); ++i) {
    int state_id = state_id_column[i];
    // Valid lookup ids populate their direct cache slot; malformed ids are
    // rejected later when a live frontier tries to use them.
    if (state_id >= 0 && state_id < static_cast<int>(state_labels.size())) {
      state_labels[state_id] = state_column[i];
    }
  }
  int maximum_phylo_node_id = max_node;
  Rcpp::CharacterVector phylo_node_labels(
    maximum_phylo_node_id + 1,
    NA_STRING
  );
  // Cache every public phylogeny node label once. Edge-matrix cells reference
  // these immutable R strings rather than calling `std::to_string` millions of
  // times on large trees.
  for (int node_id = 1; node_id <= maximum_phylo_node_id; ++node_id) {
    phylo_node_labels[node_id] = std::to_string(node_id);
  }
  Rcpp::CharacterMatrix scenario_char_matrix(input.ntips, nslice);
  Rcpp::CharacterMatrix edge_node_char_matrix(input.ntips, nslice);
  // Materialize every retained atomic batch as one public output column.
  for (int col = 0; col < nslice; ++col) {
    int src = col;
    out_time[col] = pst_v34_batch_round_time(time_vec[src], time_scale);
    out_event[col] = event_vec[src];
    out_event_tips[col] = Rcpp::wrap(event_tips[src]);
    // Copy each biological tip's five aligned ownership ids and cached labels
    // into the retained batch column.
    for (int tip = 0; tip < input.ntips; ++tip) {
      state_id_matrix(tip, col) = state_slices[src][tip];
      edge_node_matrix(tip, col) = edge_node_slices[src][tip];
      phylo_edge_matrix(tip, col) = phylo_edge_slices[src][tip];
      scenario_edge_matrix(tip, col) = scenario_edge_slices[src][tip];
      path_id_matrix(tip, col) = path_slices[src][tip];
      int state_id = state_slices[src][tip];
      // A valid state id selects its cached R label; invalid frontier state is
      // represented as `NA` so malformed input cannot index outside the cache.
      if (state_id >= 0 && state_id < static_cast<int>(state_labels.size())) {
        scenario_char_matrix(tip, col) = state_labels[state_id];
      } else {
        // Invalid state ids have no biological label in the supplied lookup.
        scenario_char_matrix(tip, col) = NA_STRING;
      }
      int edge_node_id = edge_node_slices[src][tip];
      // Frontier node ids must belong to the input phylogeny; invalid ids are
      // preserved as `NA` instead of allocating a misleading label.
      if (edge_node_id >= 1 && edge_node_id <= maximum_phylo_node_id) {
        edge_node_char_matrix(tip, col) = phylo_node_labels[edge_node_id];
      } else {
        // Invalid frontier node ids have no topology label in the input phylogeny.
        edge_node_char_matrix(tip, col) = NA_STRING;
      }
    }
  }

  // Debug metadata materialization:
  //   Snapshot mode serializes the compact records committed beside public
  //   slices. Normal traversal keeps this object NULL and allocates no per-event
  //   copies, preserving the performance contract of `debug = FALSE`.
  if (recover_debug_failure) {
    // One-to-one alignment is required before R can attach snapshot attributes
    // without inferring event identity from scalar compatibility codes.
    if (static_cast<int>(debug_step_records.size()) != nslice) {
      Rcpp::stop("V34 debug step metadata is not aligned to committed rows");
    }
    out_debug_step_metadata = Rcpp::List(nslice);
    // Step loop:
    //   Serialize one read-only committed record per TT/scenario-matrix row.
    //   Vector order preserves causal event order within simultaneous batches.
    for (int step_id = 0; step_id < nslice; ++step_id) {
      const V34DebugStepRecord& step = debug_step_records[step_id];
      out_debug_step_metadata[step_id] = Rcpp::List::create(
        Rcpp::Named("event_types") = step.event_types,
        Rcpp::Named("causal_phases") = step.causal_phases,
        Rcpp::Named("event_ids") = step.event_ids,
        Rcpp::Named("phylogeny_edges") = step.phylogeny_edges,
        Rcpp::Named("scenario_ids") = step.scenario_ids,
        Rcpp::Named("terminal_ids") = step.terminal_ids,
        Rcpp::Named("batch_size") =
          static_cast<int>(step.event_ids.size())
      );
    }
  }

  std::vector<int> segment_edge_id;
  std::vector<int> segment_index;
  std::vector<double> segment_end_times;
  std::vector<double> segment_durations;
  std::vector<int> segment_state_ids;
  std::vector<int> segment_path_ids;
  std::vector<int> segment_sizes;
  std::vector<int> segment_parent_nodes;
  std::vector<int> segment_child_nodes;
  std::vector<int> segment_is_terminal;
  int n_records = static_cast<int>(records.size()) - 1;
  Rcpp::List record_end_times(n_records);
  Rcpp::List record_durations(n_records);
  Rcpp::List record_state_ids(n_records);
  Rcpp::List record_path_ids(n_records);
  Rcpp::List record_sizes(n_records);
  Rcpp::List record_parent_nodes(n_records);
  Rcpp::List record_child_nodes(n_records);
  Rcpp::List record_tip_ids(n_records);
  Rcpp::LogicalVector record_is_terminal(n_records);
  Rcpp::IntegerVector record_initial_path_id(n_records);
  Rcpp::IntegerVector record_terminal_path_id(n_records);
  Rcpp::IntegerVector record_parent_scenario_edge_id(n_records);
  Rcpp::NumericVector record_formation_time(n_records);
  Rcpp::CharacterVector record_edge_type(n_records);
  std::vector<std::vector<int> > final_tips_by_scenario_edge(
    static_cast<std::size_t>(n_records + 1)
  );
  long long maximum_terminal_key = std::numeric_limits<long long>::min();
  std::vector<long long> terminal_key_by_tip(
    static_cast<std::size_t>(input.ntips + 1),
    std::numeric_limits<long long>::min()
  );
  std::vector<int> terminal_path_by_tip(
    static_cast<std::size_t>(input.ntips + 1),
    NA_INTEGER
  );
  for (int terminal_index = 0;
       terminal_index < static_cast<int>(terminal_tip_ids.size());
       ++terminal_index) {
    int tip_id = terminal_tip_ids[terminal_index];
    long long terminal_key = pst_v34_batch_time_key(
      terminal_times[terminal_index],
      time_scale
    );
    terminal_key_by_tip[tip_id] = terminal_key;
    terminal_path_by_tip[tip_id] = terminal_path_ids[terminal_index];
    maximum_terminal_key = std::max(maximum_terminal_key, terminal_key);
  }
  // Capture the final committed tip membership of every scenario edge from the
  // live frontier. Early terminal tips retain exact terminal records but are
  // not final survivor-topology leaves; only maximum-key endpoints label the
  // public scenario tips.
  for (int tip_id = 1; tip_id <= input.ntips; ++tip_id) {
    int scenario_edge_id = current_scenario_edge_id[tip_id];
    // Only valid committed scenario ids may receive represented biological tips.
    if (scenario_edge_id >= 1 && scenario_edge_id <= n_records &&
        terminal_key_by_tip[tip_id] == maximum_terminal_key &&
        terminal_path_by_tip[tip_id] ==
          records[scenario_edge_id].terminal_path_id) {
      final_tips_by_scenario_edge[scenario_edge_id].push_back(tip_id);
    }
  }
  std::vector<int> scenario_child_count(
    static_cast<std::size_t>(n_records + 1),
    0
  );
  for (int edge_id = 1; edge_id <= n_records; ++edge_id) {
    int parent_id = records[edge_id].parent_scenario_edge_id;
    if (parent_id >= 1 && parent_id <= n_records) {
      ++scenario_child_count[parent_id];
    }
  }
  std::vector<bool> terminal_only_leaf(
    static_cast<std::size_t>(n_records + 1),
    false
  );
  for (int edge_id = 1; edge_id <= n_records; ++edge_id) {
    terminal_only_leaf[edge_id] =
      scenario_child_count[edge_id] == 0 &&
      final_tips_by_scenario_edge[edge_id].empty();
  }
  // A scenario lineage that becomes wholly extinct remains a biological leaf
  // even when it ends before the global horizon. Retain its endpoint tips only
  // when no stable child continues that source.
  for (int terminal_index = 0;
       terminal_index < static_cast<int>(terminal_tip_ids.size());
       ++terminal_index) {
    int scenario_edge_id = terminal_scenario_edge_ids[terminal_index];
    if (scenario_edge_id >= 1 &&
        scenario_edge_id <= n_records &&
        terminal_only_leaf[scenario_edge_id] &&
        terminal_path_ids[terminal_index] ==
          records[scenario_edge_id].terminal_path_id) {
      final_tips_by_scenario_edge[scenario_edge_id].push_back(
        terminal_tip_ids[terminal_index]
      );
    }
  }
  // Serialize every scenario edge record into edge-parallel R lists while also
  // flattening its segments for downstream C++ tree materialization.
  for (int edge_id = 1; edge_id < static_cast<int>(records.size()); ++edge_id) {
    const V34ScenarioEdgeRecord& rec = records[edge_id];
    record_end_times[edge_id - 1] = Rcpp::wrap(rec.end_times);
    record_durations[edge_id - 1] = Rcpp::wrap(rec.durations);
    record_state_ids[edge_id - 1] = Rcpp::wrap(rec.state_ids);
    record_path_ids[edge_id - 1] = Rcpp::wrap(rec.path_ids);
    record_sizes[edge_id - 1] = Rcpp::wrap(rec.sizes);
    record_parent_nodes[edge_id - 1] = Rcpp::wrap(rec.parent_nodes);
    record_child_nodes[edge_id - 1] = Rcpp::wrap(rec.child_nodes);
    record_tip_ids[edge_id - 1] = Rcpp::wrap(
      final_tips_by_scenario_edge[edge_id]
    );
    record_is_terminal[edge_id - 1] = rec.is_terminal;
    record_initial_path_id[edge_id - 1] = rec.initial_path_id;
    record_terminal_path_id[edge_id - 1] = rec.terminal_path_id;
    record_parent_scenario_edge_id[edge_id - 1] =
      rec.parent_scenario_edge_id;
    record_formation_time[edge_id - 1] = rec.formation_time;
    record_edge_type[edge_id - 1] =
      v34_scenario_edge_type_label(rec.edge_type);
    // Append each biological state/path segment in edge-local chronological
    // order to the flat segment table.
    for (int seg = 0; seg < static_cast<int>(rec.state_ids.size()); ++seg) {
      segment_edge_id.push_back(edge_id);
      segment_index.push_back(seg + 1);
      segment_end_times.push_back(rec.end_times[seg]);
      segment_durations.push_back(rec.durations[seg]);
      segment_state_ids.push_back(rec.state_ids[seg]);
      segment_path_ids.push_back(rec.path_ids[seg]);
      segment_sizes.push_back(rec.sizes[seg]);
      segment_parent_nodes.push_back(rec.parent_nodes.back());
      segment_child_nodes.push_back(rec.child_nodes.back());
      segment_is_terminal.push_back(rec.is_terminal ? 1 : 0);
    }
  }

  Rcpp::List path_component_list(path_components.size() - 1);
  Rcpp::IntegerVector path_parent_out(path_components.size() - 1);
  Rcpp::IntegerVector path_added_out(path_components.size() - 1);
  Rcpp::IntegerVector path_terminal_state_out(path_components.size() - 1);
  // Serialize every canonical path's component sequence, parent path, and added
  // state under the shared 1-based registry.
  for (int path_id = 1; path_id < static_cast<int>(path_components.size()); ++path_id) {
    path_component_list[path_id - 1] = Rcpp::wrap(path_components[path_id]);
    path_parent_out[path_id - 1] = path_parent[path_id];
    path_added_out[path_id - 1] = path_added_state_id[path_id];
    path_terminal_state_out[path_id - 1] =
      path_components[path_id].empty() ? NA_INTEGER :
      path_components[path_id].back();
  }

  // Close the final terminal transaction before any output work begins. This
  // commits semantic scenarios, PUNIQ, LDIF, through support, and exposure;
  // all subsequent operations are serialization or allowed dependent formats.
  V34Clock::time_point final_commit_start = V34Clock::now();
  online_counters.flush_pending_batch();
  batch_capture_seconds += std::chrono::duration<double>(
    V34Clock::now() - final_commit_start
  ).count();
  double batch_commit_seconds = online_counters.batch_commit_seconds();
  double batch_counter_capture_seconds = std::max(
    0.0,
    batch_capture_seconds - batch_commit_seconds
  );
  // Materialize TT buffers without biological replay or dense-matrix scans.
  Rcpp::List online_tt_core = online_counters.serialize_trans(path_components);
  // Expand traversal-committed per-tip transition runs into the compatibility
  // matrix shape consumed by deterministic transition-tree serialization.
  Rcpp::List online_transition_support = online_counters.transition_support();
  bool has_cladogenesis = false;
  // Inspect phylogeny edge destinations once. Any internal child is a
  // cladogenetic boundary capable of merging or splitting scenario histories.
  for (int child_node_id : input.edge_child) {
    // Internal node ids lie above the tip-id range in an R `phylo` object.
    if (child_node_id > input.ntips) {
      has_cladogenesis = true;
      break;
    }
  }
  // Serialize path-size counters that were committed with their event batches.
  // Path components provide labels only; no scenario/path matrix is inspected.
  Rcpp::List online_size_support = online_counters.size_support(
    path_components,
    has_cladogenesis
  );
  std::vector<int> terminal_trajectory_group_by_tip =
    v34_terminal_trajectory_groups(
      records,
      current_scenario_edge_id,
      input.ntips,
      time_tolerance
    );
  V34Clock::time_point materialization_end = V34Clock::now();

  return Rcpp::List::create(
    Rcpp::Named("scenario_summary") = Rcpp::List::create(
      Rcpp::Named("state_id_matrix") = state_id_matrix,
      Rcpp::Named("edge_node_ids") = edge_node_matrix,
      Rcpp::Named("phylo_edge_ids") = phylo_edge_matrix,
      Rcpp::Named("scenario_edge_ids") = scenario_edge_matrix,
      Rcpp::Named("path_id_matrix") = path_id_matrix,
      Rcpp::Named("time_vec") = out_time,
      Rcpp::Named("event_vec") = out_event,
      Rcpp::Named("event_tips") = out_event_tips
    ),
    Rcpp::Named("scenario_mats") = Rcpp::List::create(
      Rcpp::Named("scenarios") = scenario_char_matrix,
      Rcpp::Named("edges") = edge_node_char_matrix,
      Rcpp::Named("phylo_edge_ids") = phylo_edge_matrix,
      Rcpp::Named("scenario_edge_ids") = scenario_edge_matrix,
      Rcpp::Named("time_vec") = out_time,
      Rcpp::Named("event_vec") = out_event,
      Rcpp::Named("event_tips") = out_event_tips
    ),
    Rcpp::Named("scenario_edge_records") = Rcpp::List::create(
      Rcpp::Named("edge_end_times") = record_end_times,
      Rcpp::Named("edge_durations") = record_durations,
      Rcpp::Named("edge_state_ids") = record_state_ids,
      Rcpp::Named("edge_path_ids") = record_path_ids,
      Rcpp::Named("edge_sizes") = record_sizes,
      Rcpp::Named("edge_parent_nodes") = record_parent_nodes,
      Rcpp::Named("edge_child_nodes") = record_child_nodes,
      Rcpp::Named("edge_tip_ids") = record_tip_ids,
      Rcpp::Named("edge_is_terminal") = record_is_terminal,
      Rcpp::Named("edge_initial_path_id") = record_initial_path_id,
      Rcpp::Named("edge_terminal_path_id") = record_terminal_path_id,
      Rcpp::Named("edge_parent_scenario_edge_id") =
        record_parent_scenario_edge_id,
      Rcpp::Named("edge_formation_time") = record_formation_time,
      Rcpp::Named("edge_type") = record_edge_type
    ),
    Rcpp::Named("scenario_edge_segments") = Rcpp::DataFrame::create(
      Rcpp::Named("scenario_edge_id") = segment_edge_id,
      Rcpp::Named("segment_index") = segment_index,
      Rcpp::Named("end_time") = segment_end_times,
      Rcpp::Named("duration") = segment_durations,
      Rcpp::Named("state_id") = segment_state_ids,
      Rcpp::Named("path_id") = segment_path_ids,
      Rcpp::Named("lineage_size") = segment_sizes,
      Rcpp::Named("parent_node") = segment_parent_nodes,
      Rcpp::Named("child_node") = segment_child_nodes,
      Rcpp::Named("is_terminal") = segment_is_terminal
    ),
    Rcpp::Named("path_components") = path_component_list,
    Rcpp::Named("path_parent") = path_parent_out,
    Rcpp::Named("path_added_state_id") = path_added_out,
    Rcpp::Named("path_terminal_state_id") = path_terminal_state_out,
    Rcpp::Named("terminal_trajectory_group_id_by_tip") = Rcpp::wrap(
      std::vector<int>(
        terminal_trajectory_group_by_tip.begin() + 1,
        terminal_trajectory_group_by_tip.end()
      )
    ),
    Rcpp::Named("final_frontier") = Rcpp::List::create(
      Rcpp::Named("state_id") = Rcpp::wrap(std::vector<int>(
        current_state_id.begin() + 1,
        current_state_id.end()
      )),
      Rcpp::Named("path_id") = Rcpp::wrap(std::vector<int>(
        current_path_id.begin() + 1,
        current_path_id.end()
      )),
      Rcpp::Named("scenario_edge_id") = Rcpp::wrap(std::vector<int>(
        current_scenario_edge_id.begin() + 1,
        current_scenario_edge_id.end()
      )),
      Rcpp::Named("phylo_edge_id") = Rcpp::wrap(std::vector<int>(
        current_phylo_edge_id.begin() + 1,
        current_phylo_edge_id.end()
      ))
    ),
    Rcpp::Named("online_tt_core") = online_tt_core,
    Rcpp::Named("online_transition_support") = online_transition_support,
    Rcpp::Named("online_size_support") = online_size_support,
    Rcpp::Named("debug_step_metadata") = out_debug_step_metadata,
    Rcpp::Named("traversal_size_support") = Rcpp::List::create(
      Rcpp::Named("has_cladogenesis") = has_cladogenesis
    ),
    Rcpp::Named("terminal_records") = Rcpp::DataFrame::create(
      Rcpp::Named("tip_id") = terminal_tip_ids,
      Rcpp::Named("final_phylo_edge_id") = terminal_phylo_edge_ids,
      Rcpp::Named("final_scenario_edge_id") = terminal_scenario_edge_ids,
      Rcpp::Named("stable_scenario_id") = terminal_scenario_edge_ids,
      Rcpp::Named("final_path_id") = terminal_path_ids,
      Rcpp::Named("terminal_time") = terminal_times
    ),
    Rcpp::Named("metadata") = Rcpp::List::create(
      Rcpp::Named("root_node") = root_node,
      Rcpp::Named("tree_height") = tree_height,
      Rcpp::Named("is_ultrametric") = is_ultrametric,
      Rcpp::Named("time_tolerance") = time_tolerance,
      Rcpp::Named("invariant_validation") = "every_capture_and_commit",
      Rcpp::Named("validated_capture_count") =
        online_counters.validated_capture_count(),
      Rcpp::Named("validated_commit_count") =
        online_counters.validated_commit_count(),
      Rcpp::Named("validated_reconciliation_count") =
        online_counters.validated_reconciliation_count(),
      Rcpp::Named("active_path_reconciliation") =
        "gross_ledgers_equal_post_batch_frontier",
      Rcpp::Named("public_slice_capture") = "one_per_validated_batch",
      Rcpp::Named("public_slice_count") = nslice,
      Rcpp::Named("batch_transaction_mode") =
        "prevalidated_source_removals",
      Rcpp::Named("causal_phase_order") =
        "incoming_cladogenesis_outgoing_terminal",
      Rcpp::Named("debug_failed") = debug_failed,
      Rcpp::Named("debug_attempted_batch_id") =
        debug_attempted_batch_id,
      Rcpp::Named("debug_attempted_time") = debug_attempted_time,
      Rcpp::Named("debug_failure_message") = debug_failure_message,
      Rcpp::Named("stage_timings") = Rcpp::List::create(
        Rcpp::Named("setup_schedule") = std::chrono::duration<double>(
          traversal_start - setup_start
        ).count(),
        Rcpp::Named("traversal_event_loop") = std::chrono::duration<double>(
          traversal_end - traversal_start
        ).count(),
        Rcpp::Named("batch_counter_capture") =
          batch_counter_capture_seconds,
        Rcpp::Named("batch_commit") = batch_commit_seconds,
        Rcpp::Named("scenario_batch_frontier") =
          scenario_batch_frontier_seconds,
        Rcpp::Named("scenario_batch_history_copy") =
          scenario_batch_history_copy_seconds,
        Rcpp::Named("scenario_batch_event_phases") =
          scenario_batch_event_phase_seconds,
        Rcpp::Named("scenario_batch_topology_commit") =
          scenario_batch_topology_commit_seconds,
        Rcpp::Named("scenario_batch_history_restore") =
          scenario_batch_history_restore_seconds,
        Rcpp::Named("scenario_batch_commit_core") =
          scenario_batch_commit_core_seconds,
        Rcpp::Named("scenario_batch_classification") =
          scenario_batch_classification_seconds,
        Rcpp::Named("cpp_output_materialization") =
          std::chrono::duration<double>(
            materialization_end - traversal_end
          ).count()
      )
    )
  );
}
