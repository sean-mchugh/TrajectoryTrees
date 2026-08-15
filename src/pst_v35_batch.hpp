#ifndef PST_V35_BATCH_HPP
#define PST_V35_BATCH_HPP

#include "pst_v35_types.hpp"

#include <map>
#include <set>
#include <vector>

/** Causal ordering inside one integer-time biological transaction. */
enum V35BatchPhase {
  V35_BATCH_INCOMING_ANAGENETIC = 1,
  V35_BATCH_CLADOGENESIS = 2,
  V35_BATCH_OUTGOING_ZERO_OFFSET = 3,
  V35_BATCH_TERMINAL = 4
};

/**
 * One scheduled mapped-edge boundary inside a transactional same-time batch.
 *
 * `scenario_edge_id` follows causal replacements. The batch-origin id remains
 * frozen for events already pending at the same time, allowing source capacity
 * validation and post-batch scenario resolution without event-order dependence.
 */
struct V35BatchEvent {
  int edge_id;
  int step;
  int scenario_edge_id;
  int batch_origin_scenario_edge_id;
  double time;
  int priority_tip;
  int phase;
};

/**
 * Planned source withdrawals and phase progress for one same-time batch.
 *
 * Capacities are frozen when a source scenario is first encountered. Every
 * event plans one lineage removal before its handler creates destinations.
 * The transaction is private to traversal and is discarded after validation.
 */
struct V35BatchTransaction {
  long long time_key;
  double time;
  std::map<int, int> source_capacity;
  std::map<int, int> source_removals;
  int event_count;
  int last_phase;
};

/** Frozen endpoint ownership and final survivor intent for one biological tip. */
struct V35ScenarioBatchTipIntent {
  int tip_id;
  int source_scenario_id;
  int endpoint_state_id;
  int endpoint_path_id;
  int endpoint_phylo_edge_id;
  bool terminal;
  int survivor_state_id;
  int survivor_path_id;
  int survivor_phylo_edge_id;
};

/** Scenario-dependent transition whose destination ID is resolved at commit. */
struct V35PendingScenarioTransition {
  double time;
  int from_state_id;
  int to_state_id;
  int phylo_edge_id;
  int source_path_id;
  int destination_path_id;
  int source_scenario_id;
  int representative_tip_id;
  bool phylo_edge_is_terminal;
};

/** Complete per-tip scenario intent for one integer-time transaction. */
struct V35ScenarioBatchIntent {
  long long time_key;
  double time;
  std::vector<V35ScenarioBatchTipIntent> tips;
  std::vector<V35PendingScenarioTransition> transitions;
  std::set<int> affected_source_scenario_ids;
};

/** Convert a positive time tolerance to the shared integer-key scale. */
double pst_v35_batch_time_scale(double time_tolerance);

/** Convert one absolute event time to a stable caller-scaled integer batch key. */
long long pst_v35_batch_time_key(double time, double time_scale);

/** Return the public numeric time represented by the stable batch key. */
double pst_v35_batch_round_time(double time, double time_scale);

/** Return the deterministic smallest descendant-tip tie-breaker. */
int pst_v35_batch_min_tip(const std::vector<int>& tips);

/** Select the next queue index by time, tip tie-breaker, then edge id. */
int pst_v35_batch_next_event_index(
  const std::vector<V35BatchEvent>& queue,
  double time_scale);

/**
 * Reassign mutable scenario ownership after a causal split while preserving
 * frozen batch-origin ownership for same-time scoring.
 */
void pst_v35_batch_reassign_scenario(
  std::vector<V35BatchEvent>& queue,
  int from_edge_id,
  int to_edge_id);

/** Initialize an empty source-removal transaction for one integer time key. */
V35BatchTransaction pst_v35_batch_begin(double time, double time_scale);

/**
 * Plan one lineage withdrawal before destination mutation and validate that
 * simultaneous events cannot overdraw the frozen source scenario.
 */
void pst_v35_batch_plan_source_removal(
  V35BatchTransaction& transaction,
  int scenario_edge_id,
  int source_capacity,
  int phase);

/** Validate all planned withdrawals and causal phase progression before commit. */
void pst_v35_batch_validate(const V35BatchTransaction& transaction);

/** Freeze all active endpoint members before same-key causal phases mutate them. */
V35ScenarioBatchIntent pst_v35_scenario_batch_begin(
  long long time_key,
  double time,
  const std::vector<int>& current_scenario_edge_id,
  const std::vector<int>& current_state_id,
  const std::vector<int>& current_path_id,
  const std::vector<int>& current_phylo_edge_id,
  const std::vector<int>& previously_terminal_tip_ids);

/** Return the unique mutable intent row for one biological tip. */
V35ScenarioBatchTipIntent& pst_v35_scenario_batch_tip(
  V35ScenarioBatchIntent& intent,
  int tip_id);

/** Update one surviving tip's final state/path without changing scenario identity. */
void pst_v35_scenario_batch_update_survivor(
  V35ScenarioBatchIntent& intent,
  int tip_id,
  int state_id,
  int path_id);

/** Mark one endpoint member terminal while retaining its frozen source ID. */
void pst_v35_scenario_batch_mark_terminal(
  V35ScenarioBatchIntent& intent,
  int tip_id);

/** Mark one frozen source scenario as biologically affected by this batch. */
void pst_v35_scenario_batch_mark_source_affected(
  V35ScenarioBatchIntent& intent,
  int scenario_edge_id);

#endif
