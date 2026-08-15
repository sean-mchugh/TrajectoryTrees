#ifndef PST_V35_STATE_HPP
#define PST_V35_STATE_HPP

#include "pst_v35_batch.hpp"
#include "pst_v35_types.hpp"

#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

/**
 * One biological state change observed while the traversal is applying an
 * event batch.
 *
 * Ownership: the online counter state owns these records until final TT
 * serialization. All ids remain 1-based R-facing ids. A record represents one
 * changed phylogeny edge; parallel changes therefore remain distinct CTT
 * events even when their source and destination paths are identical.
 */
struct V35OnlineChangeRecord {
  double time;
  int from_state_id;
  int to_state_id;
  int phylo_edge_id;
  int source_path_id;
  int path_id;
  int scenario_edge_before;
  int scenario_edge_after;
  bool phylo_edge_is_terminal;
};

/**
 * A committed through-time row captured from the live traversal frontier.
 *
 * The row stores counter values rather than a copy of the public scenario
 * matrices. Membership pairs exist only while the batch is pending. Commit
 * consumes them into path, scenario, and transition-run registries and releases
 * them before the compact public TT row is retained.
 */
struct V35OnlineCounterRow {
  long long time_key;
  double time;
  int support_phase;
  std::map<int, double> lineage_by_state;
  std::map<int, double> lineage_by_path;
  std::map<int, double> scenario_by_state;
  std::map<int, double> scenario_by_path;
  std::map<int, double> post_event_scenario_by_state;
  std::map<int, double> post_event_scenario_by_path;
  std::vector<double> cumulative_lineage_by_path;
  std::vector<double> cumulative_scenario_by_path;
  std::vector<double> cumulative_ctt_by_path;
  std::vector<double> cumulative_ldif_by_path;
  std::vector<double> cumulative_puniq_by_path;
  std::map<int, double> sequence_by_state;
  std::map<int, double> sequence_by_path;
  std::map<int, double> interval_lineage_by_state;
  std::map<int, double> interval_lineage_by_path;
  std::vector<double> exposure_by_path;
  std::vector<int> endpoint_tip_membership;
  std::vector<int> post_event_survivor_membership;
  std::vector<std::pair<int, int> > tip_state_memberships;
  std::vector<std::pair<int, int> > tip_edge_node_memberships;
  std::vector<std::pair<int, int> > tip_path_memberships;
  std::vector<std::pair<int, int> > scenario_path_memberships;
};

/**
 * Own traversal-time PST through-time counters and scored event buffers.
 *
 * Inputs: state lookup and root semantics are immutable construction metadata;
 * capture methods receive read-only views of the current frontier arrays.
 * Outputs: `serialize_trans()` returns the C++ TT-core payload consumed by R
 * formatting.
 * Invariants: one row per rounded time key, nonnegative counts, terminal tips
 * remain active at their exact endpoint, and transition records are never
 * reconstructed from public matrices. Failure: invalid ids or inconsistent
 * frontier arrays raise an R error before a malformed TT payload is returned.
 */
class V35LiveTTState {
 public:
  V35LiveTTState(
    const Rcpp::DataFrame& state_lookup,
    int ntips,
    int root_anchor_state_id,
    bool root_is_synthetic,
    double time_scale);

  /** Capture the initialized root frontier as the first pending TT row. */
  void initialize_root(
    double time,
    int event_phase,
    const std::vector<int>& current_state_id,
    const std::vector<int>& current_edge_node_id,
    const std::vector<int>& current_phylo_edge_id,
    const std::vector<int>& current_scenario_edge_id,
    const std::vector<int>& current_path_id);

  /** Begin one traversal-owned atomic biological batch. */
  V35BatchTransaction begin_batch(double time) const;

  void mark_terminal(int tip_id, double terminal_time);

  /** Register one canonical path before a batch can expose it publicly. */
  void register_path(
    int path_id,
    int parent_path_id,
    int added_state_id);

  void record_ctt_event(
    double time,
    int from_state_id,
    int to_state_id,
    int phylo_edge_id,
    int source_path_id,
    int destination_path_id,
    int scenario_edge_before,
    int scenario_edge_after,
    bool phylo_edge_is_terminal);

  /** Record one gross lineage withdrawal for atomic batch accounting. */
  void record_lineage_loss(double time, int path_id, int count = 1);

  /** Record one gross lineage entry for cumulative lineage accounting. */
  void record_lineage_gain(double time, int path_id, int count = 1);

  /** Record one gross scenario withdrawal for atomic batch accounting. */
  void record_scenario_loss(double time, int path_id, int count = 1);

  /** Record one gross scenario entry for cumulative scenario accounting. */
  void record_scenario_gain(double time, int path_id, int count = 1);

  /** Score one resolved post-batch leaving-edge formation as an LDIF event. */
  void record_ldif_leave_edge(
    double time,
    int source_path_id,
    int destination_path_id,
    int source_state_id,
    int destination_state_id,
    int parent_scenario_id,
    int child_scenario_id);

  /** Score the first committed public appearance of one canonical path. */
  void record_puniq_appearance(double time, int path_id);

  void capture_committed_row(
    double time,
    int event_phase,
    const std::vector<int>& current_state_id,
    const std::vector<int>& current_edge_node_id,
    const std::vector<int>& current_phylo_edge_id,
    const std::vector<int>& current_scenario_edge_id,
    const std::vector<int>& current_path_id);

  /** Atomically commit one complete pending biological/accounting row. */
  void commit_batch();

  /** Close final support state, commit its row, and validate empty ledgers. */
  void flush_pending_batch();

  Rcpp::List serialize_trans(
    const std::vector<std::vector<int> >& path_components);

  /** Serialize transition runs already committed with biological batches. */
  Rcpp::List transition_support() const;

  /** Serialize path-size counters already committed with biological batches. */
  Rcpp::List size_support(
    const std::vector<std::vector<int> >& path_components,
    bool has_cladogenesis) const;

  /** Return the number of provisional frontiers that passed ownership checks. */
  int validated_capture_count() const;

  /** Return the number of atomic rows that passed commit checks. */
  int validated_commit_count() const;

  /** Return the number of active path rows reconciled against their frontier. */
  int validated_reconciliation_count() const;

  /** Return elapsed seconds spent inside atomic counter commit operations. */
  double batch_commit_seconds() const;

 private:
  /** Reset reusable direct-index storage for one reconciliation family. */
  void begin_reconciliation_scratch();

  /** Load one sparse path map into the active reconciliation generation. */
  void load_reconciliation_values(
    const std::map<int, double>& values,
    std::vector<double>& destination);

  /** Validate endpoint-inclusive LTT against separate gross path ledgers. */
  void require_lineage_reconciliation(
    const V35OnlineCounterRow& previous,
    const V35OnlineCounterRow& current,
    const std::map<int, double>& losses,
    const std::map<int, double>& gains);

  /** Validate STT against separate gross scenario path ledgers. */
  void require_scenario_reconciliation(
    const std::map<int, double>& previous,
    const std::map<int, double>& current,
    const std::map<int, double>& losses,
    const std::map<int, double>& gains);

  /** Reset reusable semantic prefix/source grouping for one complete batch. */
  void begin_semantic_batch_scratch();

  /** Record one resolved current prefix and return true on its first carrier. */
  bool record_semantic_prefix(
    int semantic_prefix_id,
    int path_id,
    int raw_scenario_id);

  /** Record one unique source-prefix/destination-path child in sorted order. */
  void record_semantic_child(
    int source_prefix_id,
    int destination_path_id,
    int child_prefix_id);

  std::vector<int> state_ids_;
  std::vector<std::string> state_labels_by_id_;
  int ntips_;
  int root_anchor_state_id_;
  bool root_is_synthetic_;
  double time_scale_;
  std::vector<double> terminal_time_by_tip_;
  std::vector<long long> terminal_time_key_by_tip_;
  std::vector<V35OnlineCounterRow> rows_;
  bool has_pending_row_;
  V35OnlineCounterRow pending_row_;
  std::vector<V35OnlineChangeRecord> changes_;
  std::vector<V35OnlineChangeRecord> lineage_differentiating_changes_;
  std::map<long long, std::map<int, double> > lineage_losses_by_time_;
  std::map<long long, std::map<int, double> > lineage_gains_by_time_;
  std::map<long long, std::map<int, double> > scenario_losses_by_time_;
  std::map<long long, std::map<int, double> > scenario_gains_by_time_;
  std::map<long long, std::map<int, double> > ctt_gains_by_time_;
  std::map<long long, std::map<int, double> > ldif_gains_by_time_;
  std::map<long long, std::map<int, double> > puniq_gains_by_time_;
  std::size_t next_change_to_score_;
  std::vector<int> path_parent_by_id_;
  std::vector<int> path_added_state_by_id_;
  std::set<int> appeared_paths_;
  std::vector<double> path_appearance_time_;
  std::vector<int> path_appearance_id_;
  std::vector<int> path_appearance_parent_id_;
  std::vector<int> path_appearance_from_state_id_;
  std::vector<int> path_appearance_to_state_id_;
  std::vector<int> semantic_prefix_by_tip_;
  std::map<std::pair<int, int>, int> semantic_prefix_by_step_;
  int next_semantic_prefix_id_;
  std::vector<int> active_path_by_semantic_prefix_;
  std::vector<int> active_semantic_prefix_ids_;
  std::vector<int> semantic_prefix_generation_by_id_;
  std::vector<int> semantic_prefix_path_by_id_;
  std::vector<int> semantic_prefix_primary_raw_scenario_by_id_;
  std::vector<int> semantic_touched_prefix_ids_;
  std::vector<int> semantic_source_generation_by_id_;
  std::vector<int> semantic_source_first_child_by_id_;
  std::vector<int> semantic_source_child_count_by_id_;
  std::vector<int> semantic_touched_source_prefix_ids_;
  std::vector<int> semantic_child_path_ids_;
  std::vector<int> semantic_child_prefix_ids_;
  std::vector<int> semantic_child_next_indices_;
  int semantic_batch_generation_;
  std::vector<double> running_cumulative_lineage_by_path_;
  std::vector<double> running_cumulative_scenario_by_path_;
  std::vector<double> running_cumulative_ctt_by_path_;
  std::vector<double> running_cumulative_ldif_by_path_;
  std::vector<double> running_cumulative_puniq_by_path_;
  std::vector<double> running_exposure_by_path_;
  std::set<std::tuple<long long, int, int, int, int, int> > scored_ldif_groups_;
  std::vector<std::vector<int> > transition_run_states_by_tip_;
  std::vector<std::vector<int> > transition_run_edge_nodes_by_tip_;
  std::vector<std::vector<int> > transition_run_path_ids_by_tip_;
  std::vector<int> terminal_size_counted_by_tip_;
  std::map<int, int> terminal_lineage_count_by_path_;
  std::set<std::pair<int, int> > terminal_scenario_pairs_;
  std::vector<int> active_phylo_generation_;
  std::vector<int> active_scenario_generation_;
  std::vector<int> active_path_generation_;
  std::vector<int> phylo_owner_state_by_id_;
  std::vector<int> phylo_owner_path_by_id_;
  std::vector<int> phylo_owner_scenario_by_id_;
  std::vector<int> scenario_lineage_generation_;
  std::vector<int> scenario_lineage_count_by_id_;
  std::vector<int> interval_phylo_generation_;
  std::vector<int> reconciliation_generation_by_path_;
  std::vector<int> reconciliation_touched_path_ids_;
  std::vector<double> reconciliation_previous_active_by_path_;
  std::vector<double> reconciliation_previous_interval_by_path_;
  std::vector<double> reconciliation_current_active_by_path_;
  std::vector<double> reconciliation_current_interval_by_path_;
  std::vector<double> reconciliation_losses_by_path_;
  std::vector<double> reconciliation_gains_by_path_;
  std::vector<double> reconciliation_derived_by_path_;
  int reconciliation_generation_;
  int capture_generation_;
  int validated_capture_count_;
  int validated_commit_count_;
  int validated_reconciliation_count_;
  double batch_commit_seconds_;
};

#endif
