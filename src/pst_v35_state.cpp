#include "pst_v35_state.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <sstream>
#include <tuple>

namespace {

// Floating-point accounting comparisons use a narrow numeric guard. This is
// deliberately independent of the public event-time resolution because it
// protects exact counter reconciliation rather than grouping biological time.
constexpr double V35_ONLINE_NUMERIC_TOLERANCE = 1e-10;

/**
 * Convert a biological event time to the stable integer key shared by batch
 * capture and final TT rows.
 *
 * Input/output: accepts an absolute root-relative time and the caller-selected
 * scale, then returns their shared rounded key.
 */
long long v35_online_time_key(double time, double time_scale) {
  return pst_v35_batch_time_key(time, time_scale);
}

/** Return the public numeric time represented by one online key. */
double v35_online_round_time(long long time_key, double time_scale) {
  return static_cast<double>(time_key) / time_scale;
}

/**
 * Return true when an R-facing integer id is available for counter indexing.
 *
 * NA ids occur only before a tip joins the active frontier. They must not enter
 * active lineage, scenario, or path counters.
 */
bool v35_online_has_id(int value) {
  return value != NA_INTEGER && value > 0;
}

/**
 * Escape one biological state for the public `|`-delimited path grammar.
 *
 * Percent is escaped first conceptually so literal escape-like state text
 * remains distinguishable from delimiter escapes. All other bytes are copied
 * unchanged; labels remain deterministic UTF-8 byte strings.
 */
std::string v35_online_escape_path_state(const std::string& state) {
  std::string escaped;
  escaped.reserve(state.size());
  // Encode every byte that can collide with the public path grammar. After
  // each iteration, `escaped` decodes to the exact processed state prefix.
  for (char value : state) {
    // Literal percent must not be mistaken for an escape sequence.
    if (value == '%') {
      escaped += "%25";
    } else if (value == '|') {
      // The path delimiter is represented as its percent escape inside a state.
      escaped += "%7C";
    } else {
      // Ordinary state bytes have no path-grammar role and remain unchanged.
      escaped.push_back(value);
    }
  }
  return escaped;
}

/**
 * Build a path label from its canonical sequence of state ids.
 *
 * The leading separator matches the existing V30/V35 public path contract.
 * Unknown ids are rejected because silently emitting a partial path would make
 * TT columns disagree with PST path labels.
 */
std::string v35_online_path_label(
    const std::vector<int>& components,
    const std::vector<std::string>& state_labels_by_id) {
  std::ostringstream label;
  // Walk the ordered transition sequence and append each biological state;
  // after every iteration the label represents the exact path prefix seen so far.
  for (int state_id : components) {
    // A path component must resolve to a public state label before TT
    // materialization; otherwise the path registry and state lookup diverged.
    if (state_id <= 0 || state_id >= static_cast<int>(state_labels_by_id.size())) {
      Rcpp::stop("V35 online counter saw path component outside the state lookup");
    }
    label << "|" << v35_online_escape_path_state(
      state_labels_by_id[state_id]
    );
  }
  return label.str();
}

/**
 * Materialize one grouped online counter family as a legacy-shaped TT matrix.
 *
 * Rows are committed batch times; columns are the supplied semantic levels,
 * followed by `total` and `summary_time_vec`. Missing groups are explicit zeroes.
 */
Rcpp::NumericMatrix v35_online_group_matrix(
    const std::vector<V35OnlineCounterRow>& rows,
    const std::vector<int>& levels,
    const std::vector<std::string>& labels,
    const std::vector<std::map<int, double> V35OnlineCounterRow::*>& selectors) {
  // Exactly one selector identifies the requested counter family. Keeping this
  // validation here prevents a caller from accidentally combining TT families.
  if (selectors.size() != 1) {
    Rcpp::stop("V35 online matrix materialization requires one counter selector");
  }
  Rcpp::NumericMatrix out(rows.size(), levels.size() + 2);
  Rcpp::CharacterVector names(levels.size() + 2);
  const std::map<int, double> V35OnlineCounterRow::* selector = selectors.front();

  // Assign semantic state/path labels once; their order remains identical for
  // every row and therefore aligns all TT families on a common column contract.
  for (int col = 0; col < static_cast<int>(levels.size()); ++col) {
    names[col] = labels[col];
  }
  names[levels.size()] = "total";
  names[levels.size() + 1] = "summary_time_vec";

  // Serialize committed rows without consulting scenario matrices or public
  // trees; each pass writes one complete, nonnegative batch snapshot.
  for (int row = 0; row < static_cast<int>(rows.size()); ++row) {
    const std::map<int, double>& values = rows[row].*selector;
    double total = 0.0;
    // Fill every declared semantic level, using zero when it was absent from
    // the active frontier; after each pass `total` equals the written prefix sum.
    for (int col = 0; col < static_cast<int>(levels.size()); ++col) {
      std::map<int, double>::const_iterator found = values.find(levels[col]);
      // An absent state or path has zero active support at this committed time.
      if (found == values.end()) {
        out(row, col) = 0.0;
      } else {
        // A present counter is copied exactly; negative online state is a
        // traversal invariant violation and must not be serialized.
        if (found->second < 0.0) {
          Rcpp::stop("V35 online counter became negative before TT serialization");
        }
        out(row, col) = found->second;
      }
      total += out(row, col);
    }
    out(row, levels.size()) = total;
    out(row, levels.size() + 1) = rows[row].time;
  }
  Rcpp::colnames(out) = names;
  return out;
}

/**
 * Materialize one direct-index path counter family as a public TT matrix.
 *
 * Each committed vector is 1-based and may be shorter than the final path
 * registry because paths created later are implicit zeroes in earlier rows.
 * Values are copied directly; this function performs no biological replay.
 */
Rcpp::NumericMatrix v35_online_path_vector_matrix(
    const std::vector<V35OnlineCounterRow>& rows,
    const std::vector<int>& path_ids,
    const std::vector<std::string>& labels,
    const std::vector<double> V35OnlineCounterRow::* selector) {
  Rcpp::NumericMatrix out(rows.size(), path_ids.size() + 2);
  Rcpp::CharacterVector names(path_ids.size() + 2);
  // Assign the final canonical path vocabulary once for every committed row.
  for (int col = 0; col < static_cast<int>(path_ids.size()); ++col) {
    names[col] = labels[col];
  }
  names[path_ids.size()] = "total";
  names[path_ids.size() + 1] = "summary_time_vec";

  // Copy each live-committed path snapshot into its public dense row. Paths
  // absent when a row committed are trailing zeroes by registry construction.
  for (int row = 0; row < static_cast<int>(rows.size()); ++row) {
    const std::vector<double>& values = rows[row].*selector;
    double total = 0.0;
    // Direct path-id indexing avoids one ordered-map lookup per output cell.
    for (int col = 0; col < static_cast<int>(path_ids.size()); ++col) {
      int path_id = path_ids[col];
      double value = path_id < static_cast<int>(values.size()) ?
        values[path_id] : 0.0;
      if (!R_finite(value) || value < 0.0) {
        Rcpp::stop("V35 direct-index path counter became invalid");
      }
      out(row, col) = value;
      total += value;
    }
    out(row, path_ids.size()) = total;
    out(row, path_ids.size() + 1) = rows[row].time;
  }
  Rcpp::colnames(out) = names;
  return out;
}

/**
 * Convert online biological change records into the shared event-table shape.
 *
 * The table is deliberately lossless for CTT: only duplicate representations
 * of the same edge event are removed by the traversal before this boundary.
 */
Rcpp::DataFrame v35_online_change_frame(
    const std::vector<V35OnlineChangeRecord>& records) {
  std::vector<double> time;
  std::vector<int> from_state;
  std::vector<int> to_state;
  std::vector<int> phylo_edge;
  std::vector<int> source_path;
  std::vector<int> destination_path;
  std::vector<int> scenario_before;
  std::vector<int> scenario_after;

  // Preserve traversal event order for diagnostics while copying each record
  // into columnar R storage; no biological classification is recomputed here.
  for (const V35OnlineChangeRecord& record : records) {
    time.push_back(record.time);
    from_state.push_back(record.from_state_id);
    to_state.push_back(record.to_state_id);
    phylo_edge.push_back(record.phylo_edge_id);
    source_path.push_back(record.source_path_id);
    destination_path.push_back(record.path_id);
    scenario_before.push_back(record.scenario_edge_before);
    scenario_after.push_back(record.scenario_edge_after);
  }
  return Rcpp::DataFrame::create(
    Rcpp::Named("time") = time,
    Rcpp::Named("from_state_id") = from_state,
    Rcpp::Named("to_state_id") = to_state,
    Rcpp::Named("phylo_edge_id") = phylo_edge,
    Rcpp::Named("source_path_id") = source_path,
    Rcpp::Named("path_id") = destination_path,
    Rcpp::Named("scenario_edge_before") = scenario_before,
    Rcpp::Named("scenario_edge_after") = scenario_after,
    Rcpp::Named("parent_scenario_id") = scenario_before,
    Rcpp::Named("child_scenario_id") = scenario_after,
    Rcpp::Named("stringsAsFactors") = false
  );
}

/**
 * Format positive path counters as a named integer vector.
 *
 * The input uses canonical 1-based path ids as vector indexes. Zero support is
 * omitted to preserve the legacy sparse lookup contract used by size-map
 * attachment; no counter is recalculated or combined here.
 */
Rcpp::IntegerVector v35_online_named_positive_counts(
    const std::vector<int>& counts) {
  int positive_count = 0;
  // Count public entries before allocating the exact R vector length. Index
  // zero is an internal sentinel and is never emitted.
  for (int path_id = 1; path_id < static_cast<int>(counts.size()); ++path_id) {
    // Only positive support participates in sparse path-size lookups.
    if (counts[path_id] > 0) {
      ++positive_count;
    }
  }
  Rcpp::IntegerVector out(positive_count);
  Rcpp::CharacterVector names(positive_count);
  int output_index = 0;
  // Emit positive counters in canonical path-id order so names and values are
  // deterministic across event iteration orders.
  for (int path_id = 1; path_id < static_cast<int>(counts.size()); ++path_id) {
    // Zero support remains absent rather than becoming a misleading map entry.
    if (counts[path_id] <= 0) {
      continue;
    }
    out[output_index] = counts[path_id];
    names[output_index] = std::to_string(path_id);
    ++output_index;
  }
  out.attr("names") = names;
  return out;
}

/**
 * Sum one committed counter map while enforcing finite nonnegative values.
 *
 * The caller supplies the biological family name for precise failures. This
 * helper reads only the current transaction and returns its exact subtotal.
 */
double v35_online_validated_sum(
    const std::map<int, double>& counters,
    const std::string& family) {
  double total = 0.0;
  // Inspect every semantic state/path counter in this family. After each pass,
  // `total` equals the finite nonnegative subtotal of the processed entries.
  for (const std::pair<const int, double>& item : counters) {
    // Positive ids are required because zero is an internal sentinel and
    // negative/NA ids cannot own biological counter support.
    if (item.first < 1) {
      Rcpp::stop("V35 " + family + " contains a nonpositive semantic id");
    }
    // Negative, infinite, or NaN counters indicate an invalid batch delta and
    // must fail before the frontier or row can become public.
    if (!R_finite(item.second) || item.second < 0.0) {
      Rcpp::stop("V35 " + family + " contains invalid support");
    }
    total += item.second;
  }
  return total;
}

/** Sum one 1-based direct-index counter while validating every public slot. */
double v35_online_validated_sum(
    const std::vector<double>& counters,
    const std::string& family) {
  double total = 0.0;
  // Index zero is the internal path sentinel; all public path values follow it.
  for (int path_id = 1; path_id < static_cast<int>(counters.size()); ++path_id) {
    double value = counters[path_id];
    if (!R_finite(value) || value < 0.0) {
      Rcpp::stop("V35 " + family + " contains invalid support");
    }
    total += value;
  }
  return total;
}

/** Compare two accounting subtotals under the traversal numeric tolerance. */
void v35_online_require_equal_total(
    double left,
    double right,
    const std::string& family) {
  // State/path partitions of one biological family must describe the same
  // active frontier; a mismatch identifies incomplete or duplicated ownership.
  if (std::fabs(left - right) > V35_ONLINE_NUMERIC_TOLERANCE) {
    Rcpp::stop("V35 " + family + " state/path totals diverged");
  }
}

/** Return one sparse counter value, treating an absent path as zero. */
double v35_online_path_value(
    const std::map<int, double>& counters,
    int path_id) {
  std::map<int, double>::const_iterator found = counters.find(path_id);
  return found == counters.end() ? 0.0 : found->second;
}

/** Return one direct-index counter value, treating a future path as zero. */
double v35_online_path_value(
    const std::vector<double>& counters,
    int path_id) {
  return path_id > 0 && path_id < static_cast<int>(counters.size()) ?
    counters[path_id] : 0.0;
}

/** Grow one 1-based path vector before a live counter mutation. */
void v35_online_ensure_path(
    std::vector<double>& counters,
    int path_id) {
  if (path_id < 1) {
    Rcpp::stop("V35 direct-index counter received a nonpositive path id");
  }
  if (path_id >= static_cast<int>(counters.size())) {
    counters.resize(path_id + 1, 0.0);
  }
}

/** Copy one sparse active root row into its direct-index cumulative owner. */
void v35_online_initialize_running_paths(
    std::vector<double>& running,
    const std::map<int, double>& initial,
    const std::string& family) {
  running.clear();
  running.resize(1, 0.0);
  // Root occupancy establishes the initial cumulative value on each path.
  for (const std::pair<const int, double>& item : initial) {
    if (!R_finite(item.second) || item.second < 0.0) {
      Rcpp::stop("V35 " + family + " root contains invalid support");
    }
    v35_online_ensure_path(running, item.first);
    running[item.first] = item.second;
  }
}

/** Apply sparse gross gains to one traversal-owned cumulative path vector. */
void v35_online_apply_path_gains(
    std::vector<double>& running,
    const std::map<int, double>& gains,
    const std::string& family) {
  // Each sparse gain mutates only its destination path in the live owner.
  for (const std::pair<const int, double>& gain : gains) {
    if (!R_finite(gain.second) || gain.second < 0.0) {
      Rcpp::stop("V35 " + family + " gain contains invalid support");
    }
    v35_online_ensure_path(running, gain.first);
    running[gain.first] += gain.second;
  }
}

}  // namespace

/**
 * Initialize online counters before the root frontier is created.
 *
 * State labels are copied because the R lookup does not outlive the bridge
 * call. Terminal times begin at infinity, meaning every lineage is active until
 * its terminal event is committed.
 */
V35LiveTTState::V35LiveTTState(
    const Rcpp::DataFrame& state_lookup,
    int ntips,
    int root_anchor_state_id,
    bool root_is_synthetic,
    double time_scale)
    : ntips_(ntips),
      root_anchor_state_id_(root_anchor_state_id),
      root_is_synthetic_(root_is_synthetic),
      time_scale_(time_scale),
      terminal_time_by_tip_(ntips + 1, std::numeric_limits<double>::infinity()),
      terminal_time_key_by_tip_(
        ntips + 1,
        std::numeric_limits<long long>::max()
      ),
      has_pending_row_(false),
      next_change_to_score_(0),
      semantic_prefix_by_tip_(ntips + 1, 0),
      next_semantic_prefix_id_(1),
      semantic_batch_generation_(0),
      transition_run_states_by_tip_(ntips + 1),
      transition_run_edge_nodes_by_tip_(ntips + 1),
      transition_run_path_ids_by_tip_(ntips + 1),
      terminal_size_counted_by_tip_(ntips + 1, 0),
      active_phylo_generation_(ntips * 2 + 1, 0),
      active_scenario_generation_(ntips * 2 + 1, 0),
      active_path_generation_(ntips + 2, 0),
      phylo_owner_state_by_id_(ntips * 2 + 1, 0),
      phylo_owner_path_by_id_(ntips * 2 + 1, 0),
      phylo_owner_scenario_by_id_(ntips * 2 + 1, 0),
      scenario_lineage_generation_(ntips * 2 + 1, 0),
      scenario_lineage_count_by_id_(ntips * 2 + 1, 0),
      interval_phylo_generation_(ntips * 2 + 1, 0),
      reconciliation_generation_(0),
      capture_generation_(0),
      validated_capture_count_(0),
      validated_commit_count_(0),
      validated_reconciliation_count_(0),
      batch_commit_seconds_(0.0) {
  if (!std::isfinite(time_scale_) || time_scale_ <= 0.0) {
    Rcpp::stop("V35 live TT time scale must be finite and positive");
  }
  Rcpp::IntegerVector ids = state_lookup["state_id"];
  Rcpp::CharacterVector labels = state_lookup["state"];
  int max_id = 0;
  // Find the largest R-facing state id so label lookup can use direct indexing.
  for (int state_id : ids) {
    max_id = std::max(max_id, state_id);
  }
  state_labels_by_id_.resize(max_id + 1);
  // Copy lookup rows once; each iteration establishes one immutable id/label pair.
  for (int row = 0; row < ids.size(); ++row) {
    state_ids_.push_back(ids[row]);
    state_labels_by_id_[ids[row]] = Rcpp::as<std::string>(labels[row]);
  }
}

/** Begin one reusable direct-index reconciliation generation. */
void V35LiveTTState::begin_reconciliation_scratch() {
  reconciliation_touched_path_ids_.clear();
  // Generation rollover is formally handled even though one build cannot
  // realistically commit enough rows to exhaust a signed integer.
  if (reconciliation_generation_ == std::numeric_limits<int>::max()) {
    std::fill(
      reconciliation_generation_by_path_.begin(),
      reconciliation_generation_by_path_.end(),
      0
    );
    reconciliation_generation_ = 1;
  } else {
    ++reconciliation_generation_;
  }
}

/** Load one sparse path family into reusable direct-index reconciliation data. */
void V35LiveTTState::load_reconciliation_values(
    const std::map<int, double>& values,
    std::vector<double>& destination) {
  // Visit every populated semantic path once, initialize its scratch slot for
  // this generation, and write the requested family value directly by path id.
  for (const std::pair<const int, double>& item : values) {
    int path_id = item.first;
    if (path_id < 1) {
      Rcpp::stop("V35 reconciliation received a nonpositive path id");
    }
    int required_size = path_id + 1;
    // All parallel vectors grow together so any touched path has a complete,
    // zero-initializable scratch record regardless of which family found it.
    if (required_size >
        static_cast<int>(reconciliation_generation_by_path_.size())) {
      reconciliation_generation_by_path_.resize(required_size, 0);
      reconciliation_previous_active_by_path_.resize(required_size, 0.0);
      reconciliation_previous_interval_by_path_.resize(required_size, 0.0);
      reconciliation_current_active_by_path_.resize(required_size, 0.0);
      reconciliation_current_interval_by_path_.resize(required_size, 0.0);
      reconciliation_losses_by_path_.resize(required_size, 0.0);
      reconciliation_gains_by_path_.resize(required_size, 0.0);
      reconciliation_derived_by_path_.resize(required_size, 0.0);
    }
    // First touch in this generation clears stale values from every family and
    // adds the path once to the compact validation domain.
    if (reconciliation_generation_by_path_[path_id] !=
        reconciliation_generation_) {
      reconciliation_generation_by_path_[path_id] = reconciliation_generation_;
      reconciliation_touched_path_ids_.push_back(path_id);
      reconciliation_previous_active_by_path_[path_id] = 0.0;
      reconciliation_previous_interval_by_path_[path_id] = 0.0;
      reconciliation_current_active_by_path_[path_id] = 0.0;
      reconciliation_current_interval_by_path_[path_id] = 0.0;
      reconciliation_losses_by_path_[path_id] = 0.0;
      reconciliation_gains_by_path_[path_id] = 0.0;
      reconciliation_derived_by_path_[path_id] = 0.0;
    }
    destination[path_id] = item.second;
  }
}

/** Validate one endpoint-inclusive lineage row against gross path ledgers. */
void V35LiveTTState::require_lineage_reconciliation(
    const V35OnlineCounterRow& previous,
    const V35OnlineCounterRow& current,
    const std::map<int, double>& losses,
    const std::map<int, double>& gains) {
  begin_reconciliation_scratch();
  load_reconciliation_values(
    previous.lineage_by_path,
    reconciliation_previous_active_by_path_
  );
  load_reconciliation_values(
    previous.interval_lineage_by_path,
    reconciliation_previous_interval_by_path_
  );
  load_reconciliation_values(
    current.lineage_by_path,
    reconciliation_current_active_by_path_
  );
  load_reconciliation_values(
    current.interval_lineage_by_path,
    reconciliation_current_interval_by_path_
  );
  load_reconciliation_values(losses, reconciliation_losses_by_path_);
  load_reconciliation_values(gains, reconciliation_gains_by_path_);
  std::sort(
    reconciliation_touched_path_ids_.begin(),
    reconciliation_touched_path_ids_.end()
  );

  // Validate the preceding endpoint across the complete direct-index domain
  // before checking the current endpoint, preserving the former error order.
  for (int path_id : reconciliation_touched_path_ids_) {
    double previous_deferred =
      reconciliation_previous_active_by_path_[path_id] -
      reconciliation_previous_interval_by_path_[path_id];
    if (previous_deferred < -V35_ONLINE_NUMERIC_TOLERANCE) {
      Rcpp::stop(
        "V35 previous endpoint interval lineage exceeds endpoint support"
      );
    }
  }
  // Current endpoint inclusion may defer a terminal loss until the next row,
  // but its interval occupancy can never exceed its visible endpoint support.
  for (int path_id : reconciliation_touched_path_ids_) {
    double current_deferred =
      reconciliation_current_active_by_path_[path_id] -
      reconciliation_current_interval_by_path_[path_id];
    if (current_deferred < -V35_ONLINE_NUMERIC_TOLERANCE) {
      Rcpp::stop(
        "V35 current endpoint interval lineage exceeds endpoint support"
      );
    }
  }
  // Derive every ledger-owned active value and validate capacity before any
  // frontier equality check can report a later-stage reconciliation mismatch.
  for (int path_id : reconciliation_touched_path_ids_) {
    double previous_deferred =
      reconciliation_previous_active_by_path_[path_id] -
      reconciliation_previous_interval_by_path_[path_id];
    double current_deferred =
      reconciliation_current_active_by_path_[path_id] -
      reconciliation_current_interval_by_path_[path_id];
    double prior_available =
      reconciliation_previous_active_by_path_[path_id] - previous_deferred;
    double loss = reconciliation_losses_by_path_[path_id];
    double gain = reconciliation_gains_by_path_[path_id];
    if (prior_available < -V35_ONLINE_NUMERIC_TOLERANCE ||
        loss > prior_available + gain + V35_ONLINE_NUMERIC_TOLERANCE) {
      Rcpp::stop("V35 gross lineage losses exceed available path capacity");
    }
    if (current_deferred > loss + V35_ONLINE_NUMERIC_TOLERANCE) {
      Rcpp::stop("V35 terminal lineage endpoint is missing its gross loss");
    }
    double value = prior_available + gain - loss + current_deferred;
    if (value < -V35_ONLINE_NUMERIC_TOLERANCE) {
      Rcpp::stop("V35 reconciled lineage path became negative");
    }
    reconciliation_derived_by_path_[path_id] =
      value > V35_ONLINE_NUMERIC_TOLERANCE ? value : 0.0;
  }
  // Compare the complete independently captured frontier after all ledger
  // arithmetic has validated, retaining the exact path-specific failure text.
  for (int path_id : reconciliation_touched_path_ids_) {
    if (std::fabs(
          reconciliation_derived_by_path_[path_id] -
          reconciliation_current_active_by_path_[path_id]
        ) > V35_ONLINE_NUMERIC_TOLERANCE) {
      Rcpp::stop(
        "V35 lineage gross ledgers do not reconcile post-batch path " +
        std::to_string(path_id)
      );
    }
  }
}

/** Validate one active scenario row against gross scenario path ledgers. */
void V35LiveTTState::require_scenario_reconciliation(
    const std::map<int, double>& previous,
    const std::map<int, double>& current,
    const std::map<int, double>& losses,
    const std::map<int, double>& gains) {
  begin_reconciliation_scratch();
  load_reconciliation_values(
    previous,
    reconciliation_previous_active_by_path_
  );
  load_reconciliation_values(
    current,
    reconciliation_current_active_by_path_
  );
  load_reconciliation_values(losses, reconciliation_losses_by_path_);
  load_reconciliation_values(gains, reconciliation_gains_by_path_);
  std::sort(
    reconciliation_touched_path_ids_.begin(),
    reconciliation_touched_path_ids_.end()
  );

  // Validate all scenario ledger arithmetic before comparing any captured cell.
  for (int path_id : reconciliation_touched_path_ids_) {
    double prior = reconciliation_previous_active_by_path_[path_id];
    double loss = reconciliation_losses_by_path_[path_id];
    double gain = reconciliation_gains_by_path_[path_id];
    if (loss > prior + gain + V35_ONLINE_NUMERIC_TOLERANCE) {
      Rcpp::stop("V35 gross scenario losses exceed available path capacity");
    }
    double value = prior + gain - loss;
    if (value < -V35_ONLINE_NUMERIC_TOLERANCE) {
      Rcpp::stop("V35 reconciled scenario path became negative");
    }
    reconciliation_derived_by_path_[path_id] =
      value > V35_ONLINE_NUMERIC_TOLERANCE ? value : 0.0;
  }
  // Captured semantic STT remains authoritative only after it equals the
  // independently derived gross-ledger result on every touched canonical path.
  for (int path_id : reconciliation_touched_path_ids_) {
    if (std::fabs(
          reconciliation_derived_by_path_[path_id] -
          reconciliation_current_active_by_path_[path_id]
        ) > V35_ONLINE_NUMERIC_TOLERANCE) {
      Rcpp::stop(
        "V35 scenario gross ledgers do not reconcile post-batch path " +
        std::to_string(path_id)
      );
    }
  }
}

/** Reset reusable prefix/source scratch for one fully resolved batch. */
void V35LiveTTState::begin_semantic_batch_scratch() {
  semantic_touched_prefix_ids_.clear();
  semantic_touched_source_prefix_ids_.clear();
  semantic_child_path_ids_.clear();
  semantic_child_prefix_ids_.clear();
  semantic_child_next_indices_.clear();
  // A generation rollover clears only marker arrays; value arrays are reset on
  // first touch and therefore never require a full per-batch fill.
  if (semantic_batch_generation_ == std::numeric_limits<int>::max()) {
    std::fill(
      semantic_prefix_generation_by_id_.begin(),
      semantic_prefix_generation_by_id_.end(),
      0
    );
    std::fill(
      semantic_source_generation_by_id_.begin(),
      semantic_source_generation_by_id_.end(),
      0
    );
    semantic_batch_generation_ = 1;
  } else {
    ++semantic_batch_generation_;
  }
}

/** Record one current semantic prefix in direct-index batch storage. */
bool V35LiveTTState::record_semantic_prefix(
    int semantic_prefix_id,
    int path_id,
    int raw_scenario_id) {
  if (semantic_prefix_id < 1 || path_id < 1 || raw_scenario_id < 1) {
    Rcpp::stop("V35 semantic prefix scratch received an invalid id");
  }
  int required_size = semantic_prefix_id + 1;
  // Keep every prefix-indexed array parallel when recursive history creates a
  // new dense id above the current capacity.
  if (required_size >
      static_cast<int>(semantic_prefix_generation_by_id_.size())) {
    semantic_prefix_generation_by_id_.resize(required_size, 0);
    semantic_prefix_path_by_id_.resize(required_size, 0);
    semantic_prefix_primary_raw_scenario_by_id_.resize(required_size, 0);
    active_path_by_semantic_prefix_.resize(required_size, 0);
  }
  bool first_carrier =
    semantic_prefix_generation_by_id_[semantic_prefix_id] !=
      semantic_batch_generation_;
  // The first active carrier establishes this batch's canonical path and raw
  // representative; later equivalent carriers can only reduce the raw id.
  if (first_carrier) {
    semantic_prefix_generation_by_id_[semantic_prefix_id] =
      semantic_batch_generation_;
    semantic_prefix_path_by_id_[semantic_prefix_id] = path_id;
    semantic_prefix_primary_raw_scenario_by_id_[semantic_prefix_id] =
      raw_scenario_id;
    semantic_touched_prefix_ids_.push_back(semantic_prefix_id);
  } else {
    // One recursive semantic child cannot own two canonical paths.
    if (semantic_prefix_path_by_id_[semantic_prefix_id] != path_id) {
      Rcpp::stop("V35 semantic scenario child owns conflicting paths");
    }
    semantic_prefix_primary_raw_scenario_by_id_[semantic_prefix_id] = std::min(
      semantic_prefix_primary_raw_scenario_by_id_[semantic_prefix_id],
      raw_scenario_id
    );
  }
  return first_carrier;
}

/** Record one distinct child under a source prefix in destination-path order. */
void V35LiveTTState::record_semantic_child(
    int source_prefix_id,
    int destination_path_id,
    int child_prefix_id) {
  if (source_prefix_id < 1 || destination_path_id < 1 || child_prefix_id < 1) {
    Rcpp::stop("V35 semantic child scratch received an invalid id");
  }
  int required_size = source_prefix_id + 1;
  // Source-indexed ownership arrays grow together and are initialized lazily by
  // generation so inactive historical prefix ids incur no clearing cost.
  if (required_size >
      static_cast<int>(semantic_source_generation_by_id_.size())) {
    semantic_source_generation_by_id_.resize(required_size, 0);
    semantic_source_first_child_by_id_.resize(required_size, -1);
    semantic_source_child_count_by_id_.resize(required_size, 0);
  }
  if (semantic_source_generation_by_id_[source_prefix_id] !=
      semantic_batch_generation_) {
    semantic_source_generation_by_id_[source_prefix_id] =
      semantic_batch_generation_;
    semantic_source_first_child_by_id_[source_prefix_id] = -1;
    semantic_source_child_count_by_id_[source_prefix_id] = 0;
    semantic_touched_source_prefix_ids_.push_back(source_prefix_id);
  }

  int previous_index = -1;
  int child_index = semantic_source_first_child_by_id_[source_prefix_id];
  // Follow the source's short sorted child chain to deduplicate one child per
  // destination path while preserving the old ordered-map iteration order.
  while (child_index >= 0 &&
         semantic_child_path_ids_[child_index] < destination_path_id) {
    previous_index = child_index;
    child_index = semantic_child_next_indices_[child_index];
  }
  // Repeated carriers of one source/path must resolve to the identical child.
  if (child_index >= 0 &&
      semantic_child_path_ids_[child_index] == destination_path_id) {
    if (semantic_child_prefix_ids_[child_index] != child_prefix_id) {
      Rcpp::stop("V35 batch created duplicate semantic children on one path");
    }
    return;
  }

  int new_index = static_cast<int>(semantic_child_path_ids_.size());
  semantic_child_path_ids_.push_back(destination_path_id);
  semantic_child_prefix_ids_.push_back(child_prefix_id);
  semantic_child_next_indices_.push_back(child_index);
  // Insert the flat record at its sorted linked-list position for this source.
  if (previous_index < 0) {
    semantic_source_first_child_by_id_[source_prefix_id] = new_index;
  } else {
    semantic_child_next_indices_[previous_index] = new_index;
  }
  ++semantic_source_child_count_by_id_[source_prefix_id];
}

/**
 * Capture the traversal's initialized root frontier exactly once.
 *
 * Root initialization delegates to the same validated frontier capture used by
 * later batches. Existing pending or committed state is rejected so callers
 * cannot reclassify a biological batch as a second root condition.
 */
void V35LiveTTState::initialize_root(
    double time,
    int event_phase,
    const std::vector<int>& current_state_id,
    const std::vector<int>& current_edge_node_id,
    const std::vector<int>& current_phylo_edge_id,
    const std::vector<int>& current_scenario_edge_id,
    const std::vector<int>& current_path_id) {
  // The root is the sole initial frontier and must precede every committed or
  // provisional batch row.
  if (has_pending_row_ || !rows_.empty()) {
    Rcpp::stop("V35 live TT root was initialized more than once");
  }
  capture_committed_row(
    time,
    event_phase,
    current_state_id,
    current_edge_node_id,
    current_phylo_edge_id,
    current_scenario_edge_id,
    current_path_id
  );
}

/** Begin one atomic batch through the shared stable-time transaction helper. */
V35BatchTransaction V35LiveTTState::begin_batch(double time) const {
  return pst_v35_batch_begin(time, time_scale_);
}

/**
 * Mark a lineage terminal while preserving endpoint inclusion.
 *
 * The stored terminal time is consulted with a `<=` comparison during row
 * capture, so the lineage remains represented at this exact time and disappears
 * only from later intervals.
 */
void V35LiveTTState::mark_terminal(int tip_id, double terminal_time) {
  // Terminal ids must name one biological tip; accepting an invalid id would
  // corrupt both LTT endpoint counts and path-through membership.
  if (tip_id < 1 || tip_id > ntips_) {
    Rcpp::stop("V35 online counter received an invalid terminal tip id");
  }
  terminal_time_by_tip_[tip_id] = terminal_time;
  terminal_time_key_by_tip_[tip_id] = v35_online_time_key(
    terminal_time, time_scale_
  );
}

/**
 * Register one canonical path at the moment traversal creates it.
 *
 * Atomic PUNIQ scoring needs the parent path and added state before a committed
 * frontier can expose the path. Re-registration is accepted only when the
 * semantics are identical; conflicting ownership fails before batch mutation.
 */
void V35LiveTTState::register_path(
    int path_id,
    int parent_path_id,
    int added_state_id) {
  // Public path ids are positive; zero remains an internal sentinel.
  if (path_id < 1) {
    Rcpp::stop("V35 online path registry received a nonpositive path id");
  }
  // Keep parent and added-state arrays parallel whenever the registry grows.
  if (path_id >= static_cast<int>(path_parent_by_id_.size())) {
    path_parent_by_id_.resize(path_id + 1, NA_INTEGER);
    path_added_state_by_id_.resize(path_id + 1, NA_INTEGER);
  }
  bool registered = path_id == 1 ||
    path_parent_by_id_[path_id] != NA_INTEGER ||
    path_added_state_by_id_[path_id] != NA_INTEGER;
  // Existing ids must resolve to the exact same transition semantics.
  if (registered) {
    // Re-registering an id with a different parent/state would merge distinct
    // transition histories and corrupt every path-keyed counter.
    if (path_parent_by_id_[path_id] != parent_path_id ||
        path_added_state_by_id_[path_id] != added_state_id) {
      Rcpp::stop("V35 online path id has conflicting registrations");
    }
    return;
  }
  path_parent_by_id_[path_id] = parent_path_id;
  path_added_state_by_id_[path_id] = added_state_id;
}

/**
 * Record one real state-changing phylogeny-edge event during traversal.
 *
 * State continuations are rejected here rather than becoming zero-effect CTT
 * records. Scenario equality is retained because CTT counts the biological
 * change even when a singleton scenario edge continues in place; LDIF applies
 * its stricter scenario-change rule during batch finalization.
 */
void V35LiveTTState::record_ctt_event(
    double time,
    int from_state_id,
    int to_state_id,
    int phylo_edge_id,
    int source_path_id,
    int destination_path_id,
    int scenario_edge_before,
    int scenario_edge_after,
    bool phylo_edge_is_terminal) {
  // Equal states are continuations, not transition events, and must not alter
  // CTT, LDIF, PUNIQ, or the path registry.
  if (from_state_id == to_state_id) {
    return;
  }
  changes_.push_back(V35OnlineChangeRecord{
    v35_online_round_time(v35_online_time_key(time, time_scale_), time_scale_),
    from_state_id,
    to_state_id,
    phylo_edge_id,
    source_path_id,
    destination_path_id,
    scenario_edge_before,
    scenario_edge_after,
    phylo_edge_is_terminal
  });
  ctt_gains_by_time_[v35_online_time_key(time, time_scale_)][destination_path_id] += 1.0;
}

/** Record one validated gross count under its atomic batch time and path. */
static void v35_record_gross_path_change(
    std::map<long long, std::map<int, double> >& changes_by_time,
    double time,
    double time_scale,
    int path_id,
    int count,
    const std::string& family) {
  // Gross accounting accepts only positive biological multiplicities on a
  // registered public path. Zero would hide a missing event classification.
  if (path_id < 1 || count < 1 || !R_finite(time)) {
    Rcpp::stop("V35 invalid gross " + family + " accounting event");
  }
  changes_by_time[v35_online_time_key(time, time_scale)][path_id] += count;
}

/** Record one lineage leaving its current path during this batch. */
void V35LiveTTState::record_lineage_loss(
    double time,
    int path_id,
    int count) {
  v35_record_gross_path_change(
    lineage_losses_by_time_, time, time_scale_, path_id, count, "lineage loss"
  );
}

/** Record one lineage entering a path during this batch. */
void V35LiveTTState::record_lineage_gain(
    double time,
    int path_id,
    int count) {
  v35_record_gross_path_change(
    lineage_gains_by_time_, time, time_scale_, path_id, count, "lineage gain"
  );
}

/** Record one scenario leaving its current path during this batch. */
void V35LiveTTState::record_scenario_loss(
    double time,
    int path_id,
    int count) {
  v35_record_gross_path_change(
    scenario_losses_by_time_, time, time_scale_, path_id, count, "scenario loss"
  );
}

/** Record one scenario entering a path during this batch. */
void V35LiveTTState::record_scenario_gain(
    double time,
    int path_id,
    int count) {
  v35_record_gross_path_change(
    scenario_gains_by_time_, time, time_scale_, path_id, count, "scenario gain"
  );
}

/**
 * Score LDIF after one complete batch resolves a leaving scenario edge.
 *
 * The public scenario-edge id is part of the deduplication key. Repeated
 * traversal records for one public edge remain one event, while two distinct
 * public Leave edges score separately even when every biological label and
 * timestamp is identical.
 */
void V35LiveTTState::record_ldif_leave_edge(
    double time,
    int source_path_id,
    int destination_path_id,
    int source_state_id,
    int destination_state_id,
    int parent_scenario_id,
    int child_scenario_id) {
  long long time_key = v35_online_time_key(time, time_scale_);
  std::tuple<long long, int, int, int, int, int> key(
    time_key,
    source_path_id,
    destination_path_id,
    source_state_id,
    destination_state_id,
    child_scenario_id
  );
  // Duplicate records for the same public scenario edge do not score twice.
  if (!scored_ldif_groups_.insert(key).second) {
    return;
  }
  if (source_path_id < 1 || destination_path_id < 1 ||
      source_state_id < 1 || destination_state_id < 1 ||
      parent_scenario_id < 1 || child_scenario_id < 1 ||
      parent_scenario_id == child_scenario_id) {
    Rcpp::stop("V35 leaving-edge LDIF record contains an invalid id");
  }
  lineage_differentiating_changes_.push_back(V35OnlineChangeRecord{
    v35_online_round_time(time_key, time_scale_),
    source_state_id,
    destination_state_id,
    NA_INTEGER,
    source_path_id,
    destination_path_id,
    parent_scenario_id,
    child_scenario_id,
    false
  });
  ldif_gains_by_time_[time_key][destination_path_id] += 1.0;
}

/**
 * Score one canonical path's first committed public appearance.
 *
 * The operation is called only from atomic batch commit after all simultaneous
 * events are resolved. It preserves root and synthetic-root suppression while
 * keeping re-entry idempotent through the existing appeared-path registry.
 */
void V35LiveTTState::record_puniq_appearance(double time, int path_id) {
  // Every exposed path must have been registered when traversal created it.
  if (path_id < 1 ||
      path_id >= static_cast<int>(path_parent_by_id_.size()) ||
      (path_parent_by_id_[path_id] == NA_INTEGER && path_id != 1)) {
    Rcpp::stop("V35 batch exposed an unregistered path");
  }
  bool first_appearance = appeared_paths_.insert(path_id).second;
  // Previously observed paths never rescore after disappearance/reappearance.
  if (!first_appearance) {
    return;
  }
  // The root condition and synthetic-root daughters are initial conditions,
  // not transition-created public paths.
  if (path_id == 1 ||
      (root_is_synthetic_ && path_parent_by_id_[path_id] == 1)) {
    return;
  }
  int parent_path_id = path_parent_by_id_[path_id];
  int from_state_id = root_anchor_state_id_;
  // Non-root parents contribute their added terminal state as the departure.
  if (parent_path_id > 0 &&
      parent_path_id < static_cast<int>(path_added_state_by_id_.size()) &&
      path_added_state_by_id_[parent_path_id] != NA_INTEGER) {
    from_state_id = path_added_state_by_id_[parent_path_id];
  }
  long long time_key = v35_online_time_key(time, time_scale_);
  path_appearance_time_.push_back(
    v35_online_round_time(time_key, time_scale_)
  );
  path_appearance_id_.push_back(path_id);
  path_appearance_parent_id_.push_back(parent_path_id);
  path_appearance_from_state_id_.push_back(from_state_id);
  path_appearance_to_state_id_.push_back(path_added_state_by_id_[path_id]);
  puniq_gains_by_time_[time_key][path_id] += 1.0;
}

/**
 * Capture the current frontier as the committed result of one event time key.
 *
 * Repeated calls at the same rounded time replace the earlier provisional row.
 * This is the online transaction boundary: zero-duration intermediate paths and
 * order-dependent partial same-time states cannot become public TT rows.
 */
void V35LiveTTState::capture_committed_row(
    double time,
    int event_phase,
    const std::vector<int>& current_state_id,
    const std::vector<int>& current_edge_node_id,
    const std::vector<int>& current_phylo_edge_id,
    const std::vector<int>& current_scenario_edge_id,
    const std::vector<int>& current_path_id) {
  // Every frontier vector is 1-based with one slot per biological tip. Reject a
  // malformed traversal view before any counter family is updated.
  if (current_state_id.size() != static_cast<std::size_t>(ntips_ + 1) ||
      current_edge_node_id.size() != static_cast<std::size_t>(ntips_ + 1) ||
      current_phylo_edge_id.size() != static_cast<std::size_t>(ntips_ + 1) ||
      current_scenario_edge_id.size() != static_cast<std::size_t>(ntips_ + 1) ||
      current_path_id.size() != static_cast<std::size_t>(ntips_ + 1)) {
    Rcpp::stop("V35 online counter frontier vectors are not tip-aligned");
  }

  V35OnlineCounterRow row;
  row.time_key = v35_online_time_key(time, time_scale_);
  row.time = v35_online_round_time(row.time_key, time_scale_);
  // Advance reusable uniqueness generations instead of allocating ownership
  // sets for every batch. Overflow is practically unreachable, but resetting
  // every generation-indexed family keeps repeated builds formally correct.
  if (capture_generation_ == std::numeric_limits<int>::max()) {
    std::fill(active_phylo_generation_.begin(), active_phylo_generation_.end(), 0);
    std::fill(active_scenario_generation_.begin(), active_scenario_generation_.end(), 0);
    std::fill(active_path_generation_.begin(), active_path_generation_.end(), 0);
    std::fill(scenario_lineage_generation_.begin(), scenario_lineage_generation_.end(), 0);
    std::fill(interval_phylo_generation_.begin(), interval_phylo_generation_.end(), 0);
    capture_generation_ = 1;
  } else {
    // Ordinary captures advance the generation so stale ownership marks are ignored without clearing arrays.
    ++capture_generation_;
  }
  row.tip_state_memberships.reserve(static_cast<std::size_t>(ntips_));
  row.tip_edge_node_memberships.reserve(static_cast<std::size_t>(ntips_));
  row.tip_path_memberships.reserve(static_cast<std::size_t>(ntips_));
  row.scenario_path_memberships.reserve(static_cast<std::size_t>(ntips_));
  row.endpoint_tip_membership.assign(ntips_ + 1, 0);
  row.post_event_survivor_membership.assign(ntips_ + 1, 0);
  int unique_phylo_lineages = 0;
  int unique_scenarios = 0;
  std::vector<int> active_scenario_ids;
  active_scenario_ids.reserve(static_cast<std::size_t>(ntips_));

  // Visit the live tip frontier to update counter sets. This reads current
  // traversal state, not a completed scenario matrix; after each iteration the
  // sets contain each active phylogeny/scenario/path identity at most once.
  for (int tip = 1; tip <= ntips_; ++tip) {
    // Endpoint inclusion follows the same integer transaction key as batching.
    // A raw floating tolerance would retain a prior-key tip in an adjacent key.
    if (row.time_key > terminal_time_key_by_tip_[tip]) {
      continue;
    }
    int state_id = current_state_id[tip];
    int edge_node_id = current_edge_node_id[tip];
    int phylo_edge_id = current_phylo_edge_id[tip];
    int scenario_edge_id = current_scenario_edge_id[tip];
    int path_id = current_path_id[tip];
    // Tips that have not yet entered the root frontier carry NA ids and do not
    // represent an active lineage or sequence.
    if (!v35_online_has_id(state_id) ||
        !v35_online_has_id(phylo_edge_id) ||
        !v35_online_has_id(edge_node_id) ||
        !v35_online_has_id(scenario_edge_id) ||
        !v35_online_has_id(path_id)) {
      continue;
    }
    // Grow reusable generation vectors only when traversal creates an id above
    // the conservative topology/path estimate. Existing generations remain
    // valid because resize preserves all previously written slots.
    if (phylo_edge_id >= static_cast<int>(active_phylo_generation_.size())) {
      active_phylo_generation_.resize(phylo_edge_id + 1, 0);
      phylo_owner_state_by_id_.resize(phylo_edge_id + 1, 0);
      phylo_owner_path_by_id_.resize(phylo_edge_id + 1, 0);
      phylo_owner_scenario_by_id_.resize(phylo_edge_id + 1, 0);
      interval_phylo_generation_.resize(phylo_edge_id + 1, 0);
    }
    // Scenario ids can grow beyond the topology estimate when event splits
    // create online scenario children.
    if (scenario_edge_id >= static_cast<int>(active_scenario_generation_.size())) {
      active_scenario_generation_.resize(scenario_edge_id + 1, 0);
      scenario_lineage_generation_.resize(scenario_edge_id + 1, 0);
      scenario_lineage_count_by_id_.resize(scenario_edge_id + 1, 0);
    }
    // New canonical paths likewise require a generation slot before sequence
    // occupancy can be deduplicated.
    if (path_id >= static_cast<int>(active_path_generation_.size())) {
      active_path_generation_.resize(path_id + 1, 0);
    }
    // One active phylogeny edge owns one state, path, and scenario. Its first
    // descendant stores that tuple and increments LTT; later descendants must
    // match it exactly.
    if (active_phylo_generation_[phylo_edge_id] != capture_generation_) {
      active_phylo_generation_[phylo_edge_id] = capture_generation_;
      phylo_owner_state_by_id_[phylo_edge_id] = state_id;
      phylo_owner_path_by_id_[phylo_edge_id] = path_id;
      phylo_owner_scenario_by_id_[phylo_edge_id] = scenario_edge_id;
      ++unique_phylo_lineages;
      // Initialize this scenario's lineage subtotal once per generation, then
      // add the newly encountered active phylogeny lineage.
      if (scenario_lineage_generation_[scenario_edge_id] != capture_generation_) {
        scenario_lineage_generation_[scenario_edge_id] = capture_generation_;
        scenario_lineage_count_by_id_[scenario_edge_id] = 0;
        active_scenario_ids.push_back(scenario_edge_id);
      }
      ++scenario_lineage_count_by_id_[scenario_edge_id];
      row.lineage_by_state[state_id] += 1.0;
      row.lineage_by_path[path_id] += 1.0;
    } else {
      // Conflicting descendants expose a partially applied event batch and
      // cannot be allowed into any online counter family.
      if (phylo_owner_state_by_id_[phylo_edge_id] != state_id ||
          phylo_owner_path_by_id_[phylo_edge_id] != path_id ||
          phylo_owner_scenario_by_id_[phylo_edge_id] != scenario_edge_id) {
        Rcpp::stop("V35 active phylogeny lineage has conflicting ownership");
      }
    }
    // Endpoint-inclusive LTT and future-interval exposure are distinct. A
    // terminal lineage is visible now but must not occupy the interval after
    // this row; every surviving phylogeny edge contributes once.
    if (row.time_key < terminal_time_key_by_tip_[tip] &&
        interval_phylo_generation_[phylo_edge_id] != capture_generation_) {
      interval_phylo_generation_[phylo_edge_id] = capture_generation_;
      row.interval_lineage_by_state[state_id] += 1.0;
      row.interval_lineage_by_path[path_id] += 1.0;
    }
    // One active scenario edge likewise owns one state and one path. Its first
    // descendant in this generation increments both raw STT groupings.
    if (active_scenario_generation_[scenario_edge_id] != capture_generation_) {
      active_scenario_generation_[scenario_edge_id] = capture_generation_;
      ++unique_scenarios;
      row.scenario_by_state[state_id] += 1.0;
      row.scenario_by_path[path_id] += 1.0;
    }
    // A represented path has one terminal state, so its first occurrence
    // increments both state-grouped and global transition-sequence counters.
    if (active_path_generation_[path_id] != capture_generation_) {
      active_path_generation_[path_id] = capture_generation_;
      row.sequence_by_state[state_id] += 1.0;
      row.sequence_by_path[path_id] = 1.0;
    }
    row.tip_state_memberships.push_back(std::make_pair(tip, state_id));
    row.tip_edge_node_memberships.push_back(std::make_pair(tip, edge_node_id));
    row.tip_path_memberships.push_back(std::make_pair(tip, path_id));
    row.endpoint_tip_membership[tip] = 1;
    row.post_event_survivor_membership[tip] =
      row.time_key < terminal_time_key_by_tip_[tip] ? 1 : 0;
    // Prefix normalization needs the raw scenario id parallel to every active
    // tip membership, even when several tips share one scenario/path counter.
    row.scenario_path_memberships.push_back(
      std::make_pair(scenario_edge_id, path_id)
    );
  }

  // Validate all raw frontier partitions before this capture can replace a
  // provisional same-time row. These checks run for intermediate causal phases
  // as well as final post-batch views.
  double lineage_state_total = v35_online_validated_sum(
    row.lineage_by_state,
    "lineage-by-state capture"
  );
  double lineage_path_total = v35_online_validated_sum(
    row.lineage_by_path,
    "lineage-by-path capture"
  );
  v35_online_require_equal_total(
    lineage_state_total,
    lineage_path_total,
    "lineage capture"
  );
  // The unique phylogeny-owner registry is an independent count of active
  // lineages and must agree with both grouped LTT partitions.
  if (std::fabs(
        lineage_state_total - static_cast<double>(unique_phylo_lineages)
      ) > V35_ONLINE_NUMERIC_TOLERANCE) {
    Rcpp::stop("V35 captured lineage total disagrees with frontier ownership");
  }
  int scenario_lineage_total = 0;
  // Sum direct-indexed scenario lineage subtotals. Every active scenario must
  // own positive support, and their sum must recover global LTT.
  for (int scenario_edge_id : active_scenario_ids) {
    int scenario_size = scenario_lineage_count_by_id_[scenario_edge_id];
    // A zero scenario subtotal would mean scenario identity was counted without
    // any active phylogeny lineage carrier.
    if (scenario_size < 1) {
      Rcpp::stop("V35 active scenario has nonpositive lineage membership");
    }
    scenario_lineage_total += scenario_size;
  }
  // Scenario-local sizes partition the same active lineage set as LTT.
  if (std::fabs(
        lineage_state_total - static_cast<double>(scenario_lineage_total)
      ) > V35_ONLINE_NUMERIC_TOLERANCE) {
    Rcpp::stop("V35 scenario sizes do not sum to active lineages");
  }
  double raw_scenario_state_total = v35_online_validated_sum(
    row.scenario_by_state,
    "scenario-by-state capture"
  );
  double raw_scenario_path_total = v35_online_validated_sum(
    row.scenario_by_path,
    "scenario-by-path capture"
  );
  v35_online_require_equal_total(
    raw_scenario_state_total,
    raw_scenario_path_total,
    "scenario capture"
  );
  // Raw scenario identity count must agree with both grouped scenario totals.
  if (std::fabs(
        raw_scenario_state_total -
          static_cast<double>(unique_scenarios)
      ) > V35_ONLINE_NUMERIC_TOLERANCE) {
    Rcpp::stop("V35 captured scenario total disagrees with frontier ownership");
  }
  double sequence_state_total = v35_online_validated_sum(
    row.sequence_by_state,
    "sequence-by-state capture"
  );
  double sequence_path_total = v35_online_validated_sum(
    row.sequence_by_path,
    "sequence-by-path capture"
  );
  v35_online_require_equal_total(
    sequence_state_total,
    sequence_path_total,
    "sequence capture"
  );
  double interval_state_total = v35_online_validated_sum(
    row.interval_lineage_by_state,
    "future interval lineage-by-state capture"
  );
  double interval_path_total = v35_online_validated_sum(
    row.interval_lineage_by_path,
    "future interval lineage-by-path capture"
  );
  v35_online_require_equal_total(
    interval_state_total,
    interval_path_total,
    "future interval lineage capture"
  );
  // Future occupancy is a subset of endpoint-inclusive LTT and can never add
  // lineages that are absent from the public row.
  if (interval_state_total > lineage_state_total + V35_ONLINE_NUMERIC_TOLERANCE) {
    Rcpp::stop("V35 future interval occupancy exceeds endpoint LTT");
  }
  // Every represented canonical path contributes exactly one active sequence,
  // regardless of its lineage or scenario multiplicity.
  for (const std::pair<const int, double>& item : row.sequence_by_path) {
    // Any value other than one means the path registry failed to deduplicate
    // multiple active carriers of the same transition sequence.
    if (std::fabs(item.second - 1.0) > V35_ONLINE_NUMERIC_TOLERANCE) {
      Rcpp::stop("V35 active path contributes more than one sequence");
    }
  }
  ++validated_capture_count_;

  int normalized_phase = event_phase == 21 ? 2 : event_phase;
  row.support_phase = normalized_phase;

  // Advancing time proves that every causal phase at the previous key has been
  // observed. Commit its complete post-batch frontier before accepting the new
  // transaction; same-time intermediate views never become public state.
  if (has_pending_row_ && pending_row_.time_key < row.time_key) {
    commit_batch();
  }
  // Chronological time is an invariant for exposure and post-batch LDIF.
  if (has_pending_row_ && pending_row_.time_key > row.time_key) {
    Rcpp::stop("V35 online counter received event times out of order");
  }
  pending_row_ = std::move(row);
  has_pending_row_ = true;
}

/**
 * Atomically commit the complete post-batch biological/accounting frontier.
 *
 * This method performs every semantic operation that formerly replayed the row
 * records in `serialize_trans()`: exposure integration, semantic scenario-prefix
 * assignment, PUNIQ first appearance, post-batch LDIF, and through support.
 * The committed row is compacted before storage, so finalization only formats
 * already-scored counters and records.
 */
void V35LiveTTState::commit_batch() {
  // A missing pending row is a no-op; the final flush uses this after closing
  // support state even when traversal ended without another biological row.
  if (!has_pending_row_) {
    return;
  }
  std::chrono::steady_clock::time_point commit_started =
    std::chrono::steady_clock::now();
  V35OnlineCounterRow row = std::move(pending_row_);
  has_pending_row_ = false;

  // A committed public row must advance beyond the preceding batch key. Same-
  // time provisional captures replace one another and therefore never reach
  // this branch as separate committed rows.
  if (!rows_.empty() && row.time_key <= rows_.back().time_key) {
    Rcpp::stop("V35 committed batch times are not strictly increasing");
  }

  // Parallel membership arrays describe the same active-tip frontier. A
  // mismatch would mix state, path, and raw scenario ownership across tips.
  if (row.tip_state_memberships.size() != row.tip_path_memberships.size() ||
      row.tip_edge_node_memberships.size() != row.tip_path_memberships.size() ||
      row.scenario_path_memberships.size() != row.tip_path_memberships.size() ||
      row.endpoint_tip_membership.size() !=
        static_cast<std::size_t>(ntips_ + 1) ||
      row.post_event_survivor_membership.size() !=
        static_cast<std::size_t>(ntips_ + 1)) {
    Rcpp::stop("V35 pending batch memberships are not tip-aligned");
  }

  // The previous committed LTT row owns the interval ending at this batch.
  // Integrate exposure before storing the current cumulative state; terminal
  // lineages are included at their endpoint but not in future intervals.
  if (!rows_.empty()) {
    const V35OnlineCounterRow& previous = rows_.back();
    double duration = row.time - previous.time;
    // Chronological committed keys guarantee nonnegative exposure intervals.
    if (duration < -V35_ONLINE_NUMERIC_TOLERANCE) {
      Rcpp::stop("V35 batch exposure duration became negative");
    }
    // Integrate lineages surviving beyond the prior endpoint under their
    // canonical paths. State exposure is derived from this path buffer.
    for (const std::pair<const int, double>& item :
         previous.interval_lineage_by_path) {
      v35_online_ensure_path(running_exposure_by_path_, item.first);
      running_exposure_by_path_[item.first] += item.second * duration;
    }
  }
  row.exposure_by_path = running_exposure_by_path_;

  // Rebuild only this pending row's semantic scenario counters from the live
  // per-tip prefix registry. No earlier row is revisited or reinterpreted.
  std::vector<std::pair<int, int> > raw_scenario_memberships =
    std::move(row.scenario_path_memberships);
  row.scenario_by_state.clear();
  row.scenario_by_path.clear();
  row.post_event_scenario_by_state.clear();
  row.post_event_scenario_by_path.clear();
  begin_semantic_batch_scratch();
  // Update every endpoint-inclusive transition run, then assign only tips that
  // survive beyond this batch to semantic scenario children. A terminal tip
  // retains its transition endpoint without allocating a scenario identity.
  int surviving_tip_count = 0;
  for (int membership_id = 0;
       membership_id < static_cast<int>(row.tip_path_memberships.size());
       ++membership_id) {
    int tip_id = row.tip_path_memberships[membership_id].first;
    int path_id = row.tip_path_memberships[membership_id].second;
    int state_id = row.tip_state_memberships[membership_id].second;
    int edge_node_id = row.tip_edge_node_memberships[membership_id].second;
    int raw_scenario_id = raw_scenario_memberships[membership_id].first;
    int source_semantic_prefix_id = semantic_prefix_by_tip_[tip_id];
    // Every live membership must resolve to valid biological and registry ids
    // before it can mutate semantic scenario state.
    if (tip_id < 1 || tip_id > ntips_ ||
        path_id < 1 ||
        path_id >= static_cast<int>(path_parent_by_id_.size()) ||
        state_id < 1 ||
        state_id >= static_cast<int>(state_labels_by_id_.size()) ||
        !v35_online_has_id(edge_node_id) ||
        !v35_online_has_id(raw_scenario_id)) {
      Rcpp::stop("V35 pending semantic scenario contains an invalid id");
    }
    std::vector<int>& run_states = transition_run_states_by_tip_[tip_id];
    std::vector<int>& run_edge_nodes =
      transition_run_edge_nodes_by_tip_[tip_id];
    std::vector<int>& run_path_ids =
      transition_run_path_ids_by_tip_[tip_id];
    // Initialize the first committed biological path directly.
    if (run_states.empty()) {
      run_states.push_back(state_id);
      run_edge_nodes.push_back(edge_node_id);
      run_path_ids.push_back(path_id);
    } else if (run_path_ids.back() != path_id) {
      // Several biological transitions can quantize to one atomic time key and
      // return to the state visible in the previous public row. The canonical
      // path registry retains every such event. Walk from the final path back
      // to the previously committed ancestor, then append those transitions in
      // causal order at this committed edge node.
      std::vector<int> advanced_paths;
      int cursor_path_id = path_id;
      while (cursor_path_id != run_path_ids.back()) {
        if (cursor_path_id < 1 ||
            cursor_path_id >= static_cast<int>(path_parent_by_id_.size()) ||
            cursor_path_id >=
              static_cast<int>(path_added_state_by_id_.size())) {
          Rcpp::stop("V35 transition run advanced outside the path registry");
        }
        advanced_paths.push_back(cursor_path_id);
        cursor_path_id = path_parent_by_id_[cursor_path_id];
        if (cursor_path_id == NA_INTEGER || cursor_path_id < 1) {
          Rcpp::stop(
            "V35 transition run changed to a non-descendant canonical path"
          );
        }
      }
      std::reverse(advanced_paths.begin(), advanced_paths.end());
      for (int advanced_path_id : advanced_paths) {
        int advanced_state_id =
          path_added_state_by_id_[advanced_path_id];
        if (advanced_state_id < 1 ||
            advanced_state_id >=
              static_cast<int>(state_labels_by_id_.size())) {
          Rcpp::stop("V35 transition path lacks its biological state");
        }
        run_states.push_back(advanced_state_id);
        run_edge_nodes.push_back(edge_node_id);
        run_path_ids.push_back(advanced_path_id);
      }
      if (run_states.back() != state_id || run_path_ids.back() != path_id) {
        Rcpp::stop(
          "V35 transition path registry disagrees with the committed frontier"
        );
      }
    } else {
      // A true state/path continuation extends the current run to the latest
      // committed edge node without creating a transition.
      if (run_edge_nodes.empty() || run_path_ids.empty()) {
        Rcpp::stop("V35 transition run lost its edge-node ownership");
      }
      if (run_states.back() != state_id) {
        std::ostringstream diagnostic;
        diagnostic
          << "V35 transition run changed state without advancing its path"
          << " [tip_id=" << tip_id
          << ", state_id=" << state_id
          << ", previous_path_id=" << run_path_ids.back()
          << ", current_path_id=" << path_id
          << ", time_key=" << row.time_key
          << ", time=" << row.time
          << "]";
        Rcpp::stop(diagnostic.str());
      }
      run_edge_nodes.back() = edge_node_id;
    }
    if (!row.post_event_survivor_membership[tip_id]) {
      // An endpoint may quantize to the root key and replace the provisional
      // root row before it was committed. Give that endpoint the same
      // deterministic initial semantic identity without publishing it as a
      // post-event survivor.
      if (source_semantic_prefix_id < 1) {
        std::pair<int, int> endpoint_step(0, path_id);
        std::map<std::pair<int, int>, int>::iterator known_endpoint =
          semantic_prefix_by_step_.find(endpoint_step);
        int endpoint_prefix_id = 0;
        if (known_endpoint == semantic_prefix_by_step_.end()) {
          endpoint_prefix_id = next_semantic_prefix_id_++;
          semantic_prefix_by_step_[endpoint_step] = endpoint_prefix_id;
        } else {
          endpoint_prefix_id = known_endpoint->second;
        }
        semantic_prefix_by_tip_[tip_id] = endpoint_prefix_id;
      }
      continue;
    }
    ++surviving_tip_count;
    std::pair<int, int> step(source_semantic_prefix_id, path_id);
    std::map<std::pair<int, int>, int>::iterator known =
      semantic_prefix_by_step_.find(step);
    int semantic_prefix_id = 0;
    // The first occurrence of a prior-prefix/path child creates one canonical
    // semantic scenario shared by every equivalent simultaneous lineage.
    if (known == semantic_prefix_by_step_.end()) {
      semantic_prefix_id = next_semantic_prefix_id_++;
      semantic_prefix_by_step_[step] = semantic_prefix_id;
    } else {
      // Equivalent histories reuse their canonical semantic scenario id.
      semantic_prefix_id = known->second;
    }
    semantic_prefix_by_tip_[tip_id] = semantic_prefix_id;
    bool first_carrier = record_semantic_prefix(
      semantic_prefix_id,
      path_id,
      raw_scenario_id
    );
    // The first carrier of a canonical prefix contributes one active scenario;
    // equivalent simultaneous lineages do not multiply STT occupancy.
    if (first_carrier) {
      row.post_event_scenario_by_state[state_id] += 1.0;
      row.post_event_scenario_by_path[path_id] += 1.0;
    }
    // Root initialization has no source semantic scenario. Every later active
    // member records its resolved child by the pre-batch scenario and final
    // path; repeated members of one child do not multiply scenario gains.
    if (source_semantic_prefix_id > 0) {
      if (source_semantic_prefix_id >=
            static_cast<int>(active_path_by_semantic_prefix_.size()) ||
          active_path_by_semantic_prefix_[source_semantic_prefix_id] < 1) {
        Rcpp::stop("V35 batch child has no active pre-batch scenario");
      }
      record_semantic_child(
        source_semantic_prefix_id,
        path_id,
        semantic_prefix_id
      );
    }
  }

  // Rebuild endpoint-inclusive historical STT without changing survivor-only
  // semantic ownership. A terminal tip retains its prior prefix but contributes
  // its exact endpoint path/state to this public row.
  std::set<std::pair<int, int> > endpoint_scenario_paths;
  for (int membership_id = 0;
       membership_id < static_cast<int>(row.tip_path_memberships.size());
       ++membership_id) {
    int tip_id = row.tip_path_memberships[membership_id].first;
    int path_id = row.tip_path_memberships[membership_id].second;
    int state_id = row.tip_state_memberships[membership_id].second;
    int semantic_prefix_id = semantic_prefix_by_tip_[tip_id];
    if (semantic_prefix_id < 1) {
      std::ostringstream diagnostic;
      diagnostic
        << "V35 endpoint scenario has no semantic prefix"
        << " [tip_id=" << tip_id
        << ", path_id=" << path_id
        << ", time=" << row.time
        << ", time_key=" << row.time_key
        << ", post_event_survivor="
        << (row.post_event_survivor_membership[tip_id] ? 1 : 0)
        << "]";
      Rcpp::stop(diagnostic.str());
    }
    if (endpoint_scenario_paths.insert(
          std::make_pair(semantic_prefix_id, path_id)
        ).second) {
      row.scenario_by_state[state_id] += 1.0;
      row.scenario_by_path[path_id] += 1.0;
    }
  }

  // Resolve gross scenario accounting and LDIF only after every simultaneous
  // event has produced the complete child set of each pre-batch scenario.
  // This prevents a sequential lineage update from becoming a provisional
  // Stay/Leave decision while retaining distinct gains from distinct sources.
  for (int source_prefix_id : active_semantic_prefix_ids_) {
    int source_path_id = active_path_by_semantic_prefix_[source_prefix_id];
    bool has_children =
      source_prefix_id <
        static_cast<int>(semantic_source_generation_by_id_.size()) &&
      semantic_source_generation_by_id_[source_prefix_id] ==
        semantic_batch_generation_ &&
      semantic_source_child_count_by_id_[source_prefix_id] > 0;
    // A scenario absent from the post-event survivor frontier contributes one
    // gross loss and no cumulative gain.
    if (!has_children) {
      record_scenario_loss(row.time, source_path_id);
      continue;
    }
    int first_child_index =
      semantic_source_first_child_by_id_[source_prefix_id];
    int child_count = semantic_source_child_count_by_id_[source_prefix_id];
    // One resolved child is scenario continuation. A changed path is one gross
    // step, regardless of how many member lineages moved during the batch.
    if (child_count == 1) {
      int destination_path_id = semantic_child_path_ids_[first_child_index];
      if (destination_path_id != source_path_id) {
        record_scenario_loss(row.time, source_path_id);
        record_scenario_gain(row.time, destination_path_id);
      }
      continue;
    }

    bool retains_source_path = false;
    // The source-local linked list is sorted by path, preserving deterministic
    // Leave ordering without allocating an ordered map for this batch.
    for (int child_index = first_child_index;
         child_index >= 0;
         child_index = semantic_child_next_indices_[child_index]) {
      if (semantic_child_path_ids_[child_index] == source_path_id) {
        retains_source_path = true;
        break;
      }
    }
    // Without a staying child, the source scenario is withdrawn once even
    // when it differentiates into several leaving destinations.
    if (!retains_source_path) {
      record_scenario_loss(row.time, source_path_id);
    }
    // Each distinct changed-path survivor child is one scenario gain and one
    // resolved Leave. The source-path child, when present, is the unscored Stay.
    for (int child_index = first_child_index;
         child_index >= 0;
         child_index = semantic_child_next_indices_[child_index]) {
      int destination_path_id = semantic_child_path_ids_[child_index];
      int child_prefix_id = semantic_child_prefix_ids_[child_index];
      if (destination_path_id == source_path_id) {
        continue;
      }
      record_scenario_gain(row.time, destination_path_id);
      if (surviving_tip_count < 2) {
        continue;
      }
      int source_state_id = source_path_id == 1 ?
        root_anchor_state_id_ : path_added_state_by_id_[source_path_id];
      int destination_state_id = destination_path_id == 1 ?
        root_anchor_state_id_ : path_added_state_by_id_[destination_path_id];
      bool has_raw_child =
        child_prefix_id <
          static_cast<int>(semantic_prefix_generation_by_id_.size()) &&
        semantic_prefix_generation_by_id_[child_prefix_id] ==
          semantic_batch_generation_ &&
        semantic_prefix_primary_raw_scenario_by_id_[child_prefix_id] > 0;
      if (source_state_id < 1 || destination_state_id < 1 ||
          !has_raw_child) {
        Rcpp::stop("V35 resolved Leave has incomplete path or edge identity");
      }
      int parent_scenario_id =
        semantic_prefix_primary_raw_scenario_by_id_[source_prefix_id];
      if (parent_scenario_id < 1) {
        Rcpp::stop("V35 resolved Leave lacks an exact parent scenario ID");
      }
      record_ldif_leave_edge(
        row.time,
        source_path_id,
        destination_path_id,
        source_state_id,
        destination_state_id,
        parent_scenario_id,
        semantic_prefix_primary_raw_scenario_by_id_[child_prefix_id]
      );
    }
  }
  // Replace active semantic ownership only after all pre-batch sources have
  // been scored. Sorting preserves the former ordered-map source iteration.
  for (int prior_prefix_id : active_semantic_prefix_ids_) {
    active_path_by_semantic_prefix_[prior_prefix_id] = 0;
  }
  std::sort(
    semantic_touched_prefix_ids_.begin(),
    semantic_touched_prefix_ids_.end()
  );
  for (int current_prefix_id : semantic_touched_prefix_ids_) {
    active_path_by_semantic_prefix_[current_prefix_id] =
      semantic_prefix_path_by_id_[current_prefix_id];
  }
  active_semantic_prefix_ids_ = semantic_touched_prefix_ids_;
  // Visit each active tip/path pair once to retain exact terminal contributions.
  // Internal through sizes come exclusively from the final CLTT/CSTT rows.
  for (const std::pair<int, int>& membership : row.tip_path_memberships) {
    int tip_id = membership.first;
    int path_id = membership.second;
    int semantic_prefix_id = semantic_prefix_by_tip_[tip_id];
    bool survives_event = row.post_event_survivor_membership[tip_id] != 0;
    bool has_primary_raw_scenario =
      semantic_prefix_id > 0 &&
      semantic_prefix_id <
        static_cast<int>(semantic_prefix_generation_by_id_.size()) &&
      semantic_prefix_generation_by_id_[semantic_prefix_id] ==
        semantic_batch_generation_ &&
      semantic_prefix_primary_raw_scenario_by_id_[semantic_prefix_id] > 0;
    // Every survivor received a current prefix and representative in the first
    // pass. Terminal-only members intentionally retain their prior prefix.
    if (survives_event && !has_primary_raw_scenario) {
      Rcpp::stop("V35 semantic scenario has no raw representative");
    }
    bool reaches_terminal_endpoint =
      row.time_key == terminal_time_key_by_tip_[tip_id];
    // A terminal lineage contributes once at its exact endpoint. The semantic
    // prefix distinguishes independent terminal scenarios sharing one path.
    if (reaches_terminal_endpoint && semantic_prefix_id < 1) {
      Rcpp::stop("V35 terminal endpoint has no occupied semantic scenario");
    }
    if (reaches_terminal_endpoint && !terminal_size_counted_by_tip_[tip_id]) {
      ++terminal_lineage_count_by_path_[path_id];
      terminal_scenario_pairs_.insert(
        std::make_pair(path_id, semantic_prefix_id)
      );
      terminal_size_counted_by_tip_[tip_id] = 1;
    }
  }

  // Score first public path appearances from the complete post-batch frontier.
  // Paths absent from this committed slice, including zero-duration
  // intermediates, never enter the appeared-path registry.
  for (const std::pair<const int, double>& item : row.sequence_by_path) {
    record_puniq_appearance(row.time, item.first);
  }

  // Advance CTT record ownership through this committed batch. LDIF has already
  // been resolved above from complete pre/post-batch semantic scenario groups.
  while (next_change_to_score_ < changes_.size()) {
    const V35OnlineChangeRecord& change = changes_[next_change_to_score_];
    long long change_key = v35_online_time_key(change.time, time_scale_);
    // A future change belongs to a later pending transaction.
    if (change_key > row.time_key) {
      break;
    }
    // A past unscored change proves its batch was committed incompletely.
    if (change_key < row.time_key) {
      Rcpp::stop("V35 batch commit skipped a change record");
    }
    ++next_change_to_score_;
  }

  const std::map<int, double>& lineage_losses =
    lineage_losses_by_time_[row.time_key];
  const std::map<int, double>& lineage_gains =
    lineage_gains_by_time_[row.time_key];
  const std::map<int, double>& scenario_losses =
    scenario_losses_by_time_[row.time_key];
  const std::map<int, double>& scenario_gains =
    scenario_gains_by_time_[row.time_key];
  v35_online_validated_sum(lineage_losses, "gross lineage losses");
  v35_online_validated_sum(lineage_gains, "gross lineage gains");
  v35_online_validated_sum(scenario_losses, "gross scenario losses");
  v35_online_validated_sum(scenario_gains, "gross scenario gains");

  // Root occupancy is the sole initialization oracle. Every later active row
  // is derived from the previous committed row and separate gross ledgers, then
  // checked against the independently captured complete post-batch frontier.
  if (!rows_.empty()) {
    require_lineage_reconciliation(
      rows_.back(),
      row,
      lineage_losses,
      lineage_gains
    );
    require_scenario_reconciliation(
      rows_.back().post_event_scenario_by_path,
      row.post_event_scenario_by_path,
      scenario_losses,
      scenario_gains
    );
  }
  ++validated_reconciliation_count_;

  // Commit cumulative path matrices from traversal-recorded gross gains. The
  // first row is initialization and therefore exactly equals active occupancy;
  // later losses never subtract from either cumulative family.
  if (rows_.empty()) {
    v35_online_initialize_running_paths(
      running_cumulative_lineage_by_path_,
      row.lineage_by_path,
      "cumulative lineage"
    );
    v35_online_initialize_running_paths(
      running_cumulative_scenario_by_path_,
      row.post_event_scenario_by_path,
      "cumulative scenario"
    );
    // Initialization is occupancy, never a transition-event score.
    running_cumulative_ctt_by_path_.assign(1, 0.0);
    running_cumulative_ldif_by_path_.assign(1, 0.0);
    running_cumulative_puniq_by_path_.assign(1, 0.0);
  } else {
    // Every gross lineage entry increments its destination path once.
    v35_online_apply_path_gains(
      running_cumulative_lineage_by_path_,
      lineage_gains,
      "cumulative lineage"
    );
    // Scenario entry includes both leaving-edge formation and an in-place edge
    // step onto a new path; neither can be recovered from the active net row.
    v35_online_apply_path_gains(
      running_cumulative_scenario_by_path_,
      scenario_gains,
      "cumulative scenario"
    );
    v35_online_apply_path_gains(
      running_cumulative_ctt_by_path_,
      ctt_gains_by_time_[row.time_key],
      "cumulative CTT"
    );
    v35_online_apply_path_gains(
      running_cumulative_ldif_by_path_,
      ldif_gains_by_time_[row.time_key],
      "cumulative LDIF"
    );
    v35_online_apply_path_gains(
      running_cumulative_puniq_by_path_,
      puniq_gains_by_time_[row.time_key],
      "cumulative PUNIQ"
    );
  }
  // Snapshot the six direct-index running owners only after every atomic gain
  // has been applied. These values are the canonical live rows serialized later.
  row.cumulative_lineage_by_path = running_cumulative_lineage_by_path_;
  row.cumulative_scenario_by_path = running_cumulative_scenario_by_path_;
  row.cumulative_ctt_by_path = running_cumulative_ctt_by_path_;
  row.cumulative_ldif_by_path = running_cumulative_ldif_by_path_;
  row.cumulative_puniq_by_path = running_cumulative_puniq_by_path_;
  v35_online_validated_sum(
    row.cumulative_lineage_by_path,
    "committed cumulative lineage-by-path"
  );
  v35_online_validated_sum(
    row.cumulative_scenario_by_path,
    "committed cumulative scenario-by-path"
  );

  // Loss ledgers remain separate through commit even though cumulative views
  // intentionally ignore them. Validating then releasing all four sparse maps
  // prevents a batch from leaking accounting events into the next time key.
  v35_online_validated_sum(
    lineage_losses_by_time_[row.time_key],
    "committed gross lineage losses"
  );
  v35_online_validated_sum(
    scenario_losses_by_time_[row.time_key],
    "committed gross scenario losses"
  );
  lineage_losses_by_time_.erase(row.time_key);
  lineage_gains_by_time_.erase(row.time_key);
  scenario_losses_by_time_.erase(row.time_key);
  scenario_gains_by_time_.erase(row.time_key);
  ctt_gains_by_time_.erase(row.time_key);
  ldif_gains_by_time_.erase(row.time_key);
  puniq_gains_by_time_.erase(row.time_key);

  // Revalidate the finalized semantic row after prefix assignment and event
  // scoring. This catches errors introduced between raw capture and atomic
  // commit rather than relying only on the provisional frontier checks.
  double committed_lineage_state_total = v35_online_validated_sum(
    row.lineage_by_state,
    "committed lineage-by-state"
  );
  double committed_lineage_path_total = v35_online_validated_sum(
    row.lineage_by_path,
    "committed lineage-by-path"
  );
  v35_online_require_equal_total(
    committed_lineage_state_total,
    committed_lineage_path_total,
    "committed lineage"
  );
  double committed_scenario_state_total = v35_online_validated_sum(
    row.scenario_by_state,
    "committed scenario-by-state"
  );
  double committed_scenario_path_total = v35_online_validated_sum(
    row.scenario_by_path,
    "committed scenario-by-path"
  );
  v35_online_require_equal_total(
    committed_scenario_state_total,
    committed_scenario_path_total,
    "committed scenario"
  );
  double committed_post_event_scenario_state_total = v35_online_validated_sum(
    row.post_event_scenario_by_state,
    "committed post-event scenario-by-state"
  );
  double committed_post_event_scenario_path_total = v35_online_validated_sum(
    row.post_event_scenario_by_path,
    "committed post-event scenario-by-path"
  );
  v35_online_require_equal_total(
    committed_post_event_scenario_state_total,
    committed_post_event_scenario_path_total,
    "committed post-event scenario"
  );
  double committed_sequence_state_total = v35_online_validated_sum(
    row.sequence_by_state,
    "committed sequence-by-state"
  );
  double committed_sequence_path_total = v35_online_validated_sum(
    row.sequence_by_path,
    "committed sequence-by-path"
  );
  v35_online_require_equal_total(
    committed_sequence_state_total,
    committed_sequence_path_total,
    "committed sequence"
  );
  // Semantic scenario and path totals cannot exceed active lineage support;
  // each scenario/path must own at least one represented lineage.
  if (committed_scenario_state_total >
        committed_lineage_state_total + V35_ONLINE_NUMERIC_TOLERANCE ||
      committed_sequence_state_total >
        committed_lineage_state_total + V35_ONLINE_NUMERIC_TOLERANCE) {
    Rcpp::stop("V35 committed scenario or sequence support exceeds lineages");
  }
  v35_online_validated_sum(
    row.exposure_by_path,
    "committed path exposure"
  );
  ++validated_commit_count_;

  // Per-tip memberships have now updated every live registry and support set.
  // Release them before storing the compact committed row; retaining them would
  // recreate the post-traversal reconstruction pass this architecture removes.
  std::vector<std::pair<int, int> >().swap(row.tip_state_memberships);
  std::vector<std::pair<int, int> >().swap(row.tip_edge_node_memberships);
  std::vector<std::pair<int, int> >().swap(row.tip_path_memberships);
  std::vector<int>().swap(row.endpoint_tip_membership);
  std::vector<int>().swap(row.post_event_survivor_membership);
  std::vector<std::pair<int, int> >().swap(raw_scenario_memberships);
  // State maps were needed only for same-batch cross-checks. Public state
  // surfaces derive from canonical path matrices, so release duplicate maps
  // before retaining this committed row for serialization.
  row.lineage_by_state.clear();
  row.scenario_by_state.clear();
  row.sequence_by_state.clear();
  row.interval_lineage_by_state.clear();
  rows_.push_back(std::move(row));
  batch_commit_seconds_ += std::chrono::duration<double>(
    std::chrono::steady_clock::now() - commit_started
  ).count();
}

/**
 * Close the final pending support phase and biological batch.
 *
 * Traversal calls this after terminal processing and before any PST/TT/SS
 * serialization. Every recorded state change must then belong to a committed
 * batch; otherwise an interrupted transaction is reported.
 */
void V35LiveTTState::flush_pending_batch() {
  commit_batch();
  // All real changes must have been scored when their batch committed.
  if (next_change_to_score_ != changes_.size()) {
    Rcpp::stop("V35 traversal ended with uncommitted change records");
  }
  if (!lineage_losses_by_time_.empty() ||
      !lineage_gains_by_time_.empty() ||
      !scenario_losses_by_time_.empty() ||
      !scenario_gains_by_time_.empty() ||
      !ctt_gains_by_time_.empty() ||
      !ldif_gains_by_time_.empty() ||
      !puniq_gains_by_time_.empty()) {
    Rcpp::stop("V35 traversal ended with uncommitted gross accounting events");
  }
}

/**
 * Serialize per-tip transition runs committed by atomic biological batches.
 *
 * Inputs are the immutable state lookup and the live run registries owned by
 * this counter state. Output matches the legacy transition-tree support shape.
 * Ownership remains read-only: shorter terminal histories are represented by
 * their committed final biological path and require no synthetic state.
 * Failure is reported if a tip has no committed run or if parallel run vectors
 * diverge, because either condition means traversal lost transition topology.
 */
Rcpp::List V35LiveTTState::transition_support() const {
  if (rows_.empty()) {
    Rcpp::stop("V35 transition support requested before cumulative TT commit");
  }
  const std::vector<double>& cumulative_lineage_by_path =
    rows_.back().cumulative_lineage_by_path;

  std::map<int, int> state_by_path;
  std::map<int, std::set<int> > support_tips_by_path;
  std::map<int, std::set<int> > terminal_tips_by_path;
  // Fold each tip's already-committed run registry into canonical path records.
  // Every iteration reads one aligned run sequence and writes path-keyed sets;
  // no completed matrix or public tree is created or scanned.
  for (int tip_id = 1; tip_id <= ntips_; ++tip_id) {
    const std::vector<int>& states = transition_run_states_by_tip_[tip_id];
    const std::vector<int>& edge_nodes =
      transition_run_edge_nodes_by_tip_[tip_id];
    const std::vector<int>& path_ids =
      transition_run_path_ids_by_tip_[tip_id];
    // Every initialized tip must own one state, edge node, and path per run.
    if (states.empty() ||
        states.size() != edge_nodes.size() ||
        states.size() != path_ids.size()) {
      Rcpp::stop("V35 transition support contains incomplete tip runs");
    }
    // Register every run directly under its canonical path. The initial run
    // uses one synthetic support node, preserving the established root size of
    // one; later runs retain their committed phylogeny edge-node identities.
    for (int run_id = 0; run_id < static_cast<int>(states.size()); ++run_id) {
      int path_id = path_ids[run_id];
      int state_id = states[run_id];
      // A path may never change biological state across represented tips.
      std::map<int, int>::const_iterator known_state = state_by_path.find(path_id);
      // Conflicting ownership means committed tip runs disagree, so fail before
      // publishing one path edge with an ambiguous terminal state.
      if (known_state != state_by_path.end() && known_state->second != state_id) {
        Rcpp::stop("V35 transition path owns conflicting biological states");
      }
      state_by_path[path_id] = state_id;
      support_tips_by_path[path_id].insert(tip_id);
    }
    int final_path_id = path_ids.back();
    // A tip marked terminal by traversal contributes to the terminal history
    // group for its final canonical path. Partial debug trajectories instead
    // leave this set empty and use active support tips as provisional labels.
    if (R_finite(terminal_time_by_tip_[tip_id])) {
      terminal_tips_by_path[final_path_id].insert(tip_id);
    }
  }

  std::vector<int> emitted_paths;
  // Preserve canonical path-id order while excluding registry entries that
  // never appeared in a committed transition run.
  for (const std::pair<const int, int>& item : state_by_path) {
    emitted_paths.push_back(item.first);
  }
  Rcpp::IntegerVector edge_parent_path_id(emitted_paths.size());
  Rcpp::IntegerVector edge_path_id(emitted_paths.size());
  Rcpp::IntegerVector edge_state_id(emitted_paths.size());
  Rcpp::IntegerVector edge_lineage_size(emitted_paths.size());
  Rcpp::List edge_tip_ids(emitted_paths.size());
  // Serialize one compact transition-tree edge record per committed path.
  // Numeric support and represented tips were fixed at batch commit above.
  for (int record_id = 0;
       record_id < static_cast<int>(emitted_paths.size());
       ++record_id) {
    int path_id = emitted_paths[record_id];
    // Every emitted path must have a registered parent slot, including an NA
    // initial parent; otherwise traversal and counter registries diverged.
    if (path_id < 1 ||
        path_id >= static_cast<int>(path_parent_by_id_.size())) {
      Rcpp::stop("V35 transition support contains an unknown path id");
    }
    edge_parent_path_id[record_id] = path_parent_by_id_[path_id];
    edge_path_id[record_id] = path_id;
    edge_state_id[record_id] = state_by_path[path_id];
    edge_lineage_size[record_id] = static_cast<int>(std::round(
      v35_online_path_value(cumulative_lineage_by_path, path_id)
    ));
    if (edge_lineage_size[record_id] < 1) {
      Rcpp::stop("V35 transition path has no final cumulative lineage support");
    }
    edge_tip_ids[record_id] = Rcpp::wrap(
      std::vector<int>(
        support_tips_by_path[path_id].begin(),
        support_tips_by_path[path_id].end()
      )
    );
  }

  Rcpp::IntegerVector terminal_path_id(terminal_tips_by_path.size());
  Rcpp::List terminal_tip_ids(terminal_tips_by_path.size());
  int terminal_record_id = 0;
  // Serialize one terminal-history group per final path. Grouping is already
  // canonical because each tip stored its live final path at batch commit.
  for (const std::pair<const int, std::set<int> >& item :
       terminal_tips_by_path) {
    terminal_path_id[terminal_record_id] = item.first;
    terminal_tip_ids[terminal_record_id] = Rcpp::wrap(
      std::vector<int>(item.second.begin(), item.second.end())
    );
    ++terminal_record_id;
  }

  Rcpp::IntegerVector path_parent_by_id(path_parent_by_id_.size() - 1);
  // Preserve the complete canonical parent registry so a formatter can bypass
  // zero-duration paths that never appeared in a committed run without
  // reconstructing ancestry from labels or tip histories.
  for (int path_id = 1;
       path_id < static_cast<int>(path_parent_by_id_.size());
       ++path_id) {
    path_parent_by_id[path_id - 1] = path_parent_by_id_[path_id];
  }

  Rcpp::CharacterVector state_labels_by_id(state_ids_.size());
  Rcpp::CharacterVector state_label_names(state_ids_.size());
  // Publish biological labels in stable state-id order.
  for (int index = 0; index < static_cast<int>(state_ids_.size()); ++index) {
    int state_id = state_ids_[index];
    state_labels_by_id[index] = state_labels_by_id_[state_id];
    state_label_names[index] = std::to_string(state_id);
  }
  state_labels_by_id.attr("names") = state_label_names;

  std::set<std::string> emitted_state_labels;
  // Build the mapped-edge label domain from states actually represented by
  // compact path records.
  for (const std::pair<const int, int>& item : state_by_path) {
    int state_id = item.second;
    // Every emitted biological state must resolve through the constructor lookup.
    if (state_id < 1 ||
        state_id >= static_cast<int>(state_labels_by_id_.size()) ||
        state_labels_by_id_[state_id].empty()) {
      Rcpp::stop("V35 transition run contains an unknown state id");
    }
    emitted_state_labels.insert(state_labels_by_id_[state_id]);
  }
  Rcpp::CharacterVector expanded_state_space(
    emitted_state_labels.size() * static_cast<std::size_t>(ntips_)
  );
  int expanded_index = 0;
  // Expand each emitted state across every possible lineage support so public
  // state-size maps retain the existing deterministic column domain.
  for (const std::string& state_label : emitted_state_labels) {
    // Scenario sizes are positive and cannot exceed the number of phylogeny tips.
    for (int size = 1; size <= ntips_; ++size) {
      expanded_state_space[expanded_index++] =
        state_label + "_" + std::to_string(size);
    }
  }

  Rcpp::List support = Rcpp::List::create(
    Rcpp::Named("edge_parent_path_id") = edge_parent_path_id,
    Rcpp::Named("edge_path_id") = edge_path_id,
    Rcpp::Named("edge_state_id") = edge_state_id,
    Rcpp::Named("edge_lineage_size") = edge_lineage_size,
    Rcpp::Named("edge_tip_ids") = edge_tip_ids,
    Rcpp::Named("terminal_path_id") = terminal_path_id,
    Rcpp::Named("terminal_tip_ids") = terminal_tip_ids,
    Rcpp::Named("path_parent_by_id") = path_parent_by_id,
    Rcpp::Named("state_labels_by_id") = state_labels_by_id,
    Rcpp::Named("expanded_state_space") = expanded_state_space,
    Rcpp::Named("construction_phase") = "atomic_batch_commit"
  );
  std::ostringstream digest;
  digest << "v35-transition-support:";
  const auto append_integer_vector = [&digest](
      const Rcpp::IntegerVector& values) {
    digest << values.size() << ":";
    for (int value : values) {
      digest << value << ",";
    }
    digest << ";";
  };
  append_integer_vector(edge_parent_path_id);
  append_integer_vector(edge_path_id);
  append_integer_vector(edge_state_id);
  append_integer_vector(edge_lineage_size);
  for (int record_id = 0; record_id < edge_tip_ids.size(); ++record_id) {
    append_integer_vector(Rcpp::as<Rcpp::IntegerVector>(
      edge_tip_ids[record_id]
    ));
  }
  append_integer_vector(terminal_path_id);
  for (int record_id = 0; record_id < terminal_tip_ids.size(); ++record_id) {
    append_integer_vector(Rcpp::as<Rcpp::IntegerVector>(
      terminal_tip_ids[record_id]
    ));
  }
  append_integer_vector(path_parent_by_id);
  digest << state_labels_by_id.size() << ":";
  for (int label_id = 0; label_id < state_labels_by_id.size(); ++label_id) {
    std::string label = Rcpp::as<std::string>(state_labels_by_id[label_id]);
    digest << label.size() << ":" << label << ";";
  }
  support.attr("pst_v35_transition_support_digest") = digest.str();
  return support;
}

/**
 * Serialize path-size support committed by atomic biological batches.
 *
 * Path components provide labels only. Internal numerical compatibility fields
 * are direct copies of final CLTT/CSTT; terminal fields come from the exact
 * endpoint ledgers. No historical occupancy or transition matrix is rescored.
 */
Rcpp::List V35LiveTTState::size_support(
    const std::vector<std::vector<int> >& path_components,
    bool has_cladogenesis) const {
  (void)has_cladogenesis;
  int path_count = static_cast<int>(path_components.size()) - 1;
  Rcpp::CharacterVector labels_by_id(path_count);
  // Format one label per canonical path. This loop reads the path registry and
  // writes labels only; all numeric support remains traversal-owned.
  for (int path_id = 1; path_id <= path_count; ++path_id) {
    labels_by_id[path_id - 1] = v35_online_path_label(
      path_components[path_id],
      state_labels_by_id_
    );
  }

  std::vector<int> path_tail_counts(path_count + 1, 0);
  std::vector<int> path_through_lineage_counts(path_count + 1, 0);
  std::vector<int> path_through_scenario_counts(path_count + 1, 0);
  std::vector<int> terminal_lineage_counts(path_count + 1, 0);
  std::vector<int> terminal_scenario_counts(path_count + 1, 0);
  std::vector<int> cumulative_lineage_counts(path_count + 1, 0);
  std::vector<int> cumulative_scenario_counts(path_count + 1, 0);
  if (rows_.empty()) {
    Rcpp::stop("V35 size support requested before cumulative TT commit");
  }
  const V35OnlineCounterRow& final_cumulative_row = rows_.back();
  for (int path_id = 1; path_id <= path_count; ++path_id) {
    cumulative_lineage_counts[path_id] = static_cast<int>(std::round(
      v35_online_path_value(
        final_cumulative_row.cumulative_lineage_by_path,
        path_id
      )
    ));
    cumulative_scenario_counts[path_id] = static_cast<int>(std::round(
      v35_online_path_value(
        final_cumulative_row.cumulative_scenario_by_path,
        path_id
      )
    ));
    // All compatibility through/tail values use the same canonical cumulative
    // source as public transition-tree internal sizes and SS by-path values.
    path_tail_counts[path_id] = cumulative_lineage_counts[path_id];
    path_through_lineage_counts[path_id] = cumulative_lineage_counts[path_id];
    path_through_scenario_counts[path_id] = cumulative_scenario_counts[path_id];
  }
  // Copy endpoint lineage counts already committed at each tip's terminal time.
  for (const std::pair<const int, int>& item : terminal_lineage_count_by_path_) {
    int path_id = item.first;
    // Terminal counters share the same canonical path registry as TT rows.
    if (path_id < 1 || path_id > path_count) {
      Rcpp::stop("V35 terminal size support contains an unknown path id");
    }
    terminal_lineage_counts[path_id] = item.second;
  }
  // Aggregate unique `(path, semantic terminal scenario)` pairs. Each pair was
  // inserted once when its biological tip reached its exact endpoint.
  for (const std::pair<int, int>& terminal_pair : terminal_scenario_pairs_) {
    int path_id = terminal_pair.first;
    // A terminal scenario without a registered path is an ownership failure.
    if (path_id < 1 || path_id > path_count) {
      Rcpp::stop("V35 terminal scenario support contains an unknown path id");
    }
    ++terminal_scenario_counts[path_id];
  }

  // No production consumer uses the former independently scored transition-
  // scenario support matrix. Retain only an empty compatibility field.
  Rcpp::NumericMatrix transition_counts(0, 0);

  Rcpp::IntegerVector empty_root_overrides(0);
  empty_root_overrides.attr("names") = Rcpp::CharacterVector(0);
  return Rcpp::List::create(
    Rcpp::Named("labels_by_id") = labels_by_id,
    Rcpp::Named("path_tail_lineage_counts_by_id") =
      v35_online_named_positive_counts(path_tail_counts),
    Rcpp::Named("path_through_lineage_counts_by_id") =
      v35_online_named_positive_counts(path_through_lineage_counts),
    Rcpp::Named("path_through_scenario_counts_by_id") =
      v35_online_named_positive_counts(path_through_scenario_counts),
    Rcpp::Named("path_terminal_lineage_counts_by_id") =
      v35_online_named_positive_counts(terminal_lineage_counts),
    Rcpp::Named("path_terminal_scenario_counts_by_id") =
      v35_online_named_positive_counts(terminal_scenario_counts),
    Rcpp::Named("path_transition_counts_by_id") = transition_counts,
    Rcpp::Named("path_cumulative_lineage_counts_by_id") =
      v35_online_named_positive_counts(cumulative_lineage_counts),
    Rcpp::Named("path_cumulative_scenario_counts_by_id") =
      v35_online_named_positive_counts(cumulative_scenario_counts),
    Rcpp::Named("root_lineage_size_overrides") = empty_root_overrides,
    Rcpp::Named("construction_phase") = "atomic_batch_commit"
  );
}

/**
 * Serialize all online counter families after traversal finishes.
 *
 * Finalization formats counters, event buffers, and path support that were
 * already completed by atomic batch commit. It performs no biological scoring,
 * identity reconstruction, dense-matrix scan, or public-tree interpretation.
 */
Rcpp::List V35LiveTTState::serialize_trans(
    const std::vector<std::vector<int> >& path_components) {
  // Traversal must emit at least the initialized root row; an empty TT object
  // would violate both public shape and exposure integration contracts.
  if (rows_.empty()) {
    Rcpp::stop("V35 online counter finalized without a root row");
  }

  std::vector<int> path_ids;
  std::vector<std::string> path_labels;
  // Materialize the canonical path registry in id order. Path zero is an
  // internal sentinel and is deliberately omitted from all public matrices.
  for (int path_id = 1; path_id < static_cast<int>(path_components.size()); ++path_id) {
    path_ids.push_back(path_id);
    path_labels.push_back(v35_online_path_label(path_components[path_id], state_labels_by_id_));
  }
  // Batch commit has already assigned semantic scenarios, scored events, and
  // integrated exposure. Serialization reads the compact committed rows only.
  const std::vector<V35OnlineCounterRow>& committed_rows = rows_;

  Rcpp::NumericMatrix lineage_path = v35_online_group_matrix(
    committed_rows, path_ids, path_labels, {&V35OnlineCounterRow::lineage_by_path});
  Rcpp::NumericMatrix cumulative_lineage_path = v35_online_path_vector_matrix(
    committed_rows,
    path_ids,
    path_labels,
    &V35OnlineCounterRow::cumulative_lineage_by_path
  );
  Rcpp::NumericMatrix scenario_path = v35_online_group_matrix(
    committed_rows, path_ids, path_labels, {&V35OnlineCounterRow::scenario_by_path});
  Rcpp::NumericMatrix cumulative_scenario_path = v35_online_path_vector_matrix(
    committed_rows,
    path_ids,
    path_labels,
    &V35OnlineCounterRow::cumulative_scenario_by_path
  );
  Rcpp::NumericMatrix sequence_path = v35_online_group_matrix(
    committed_rows, path_ids, path_labels, {&V35OnlineCounterRow::sequence_by_path});
  Rcpp::NumericMatrix ctt_path = v35_online_path_vector_matrix(
    committed_rows,
    path_ids,
    path_labels,
    &V35OnlineCounterRow::cumulative_ctt_by_path
  );
  Rcpp::NumericMatrix ldif_path = v35_online_path_vector_matrix(
    committed_rows,
    path_ids,
    path_labels,
    &V35OnlineCounterRow::cumulative_ldif_by_path
  );
  Rcpp::NumericMatrix puniq_path = v35_online_path_vector_matrix(
    committed_rows,
    path_ids,
    path_labels,
    &V35OnlineCounterRow::cumulative_puniq_by_path
  );

  // Reorder already-scored PUNIQ records by canonical path id for stable public
  // diagnostics. No first-appearance decision is made during serialization.
  std::map<int, int> appearance_index_by_path;
  // Index each batch-committed path record by its unique canonical id.
  for (int record_id = 0;
       record_id < static_cast<int>(path_appearance_id_.size());
       ++record_id) {
    appearance_index_by_path[path_appearance_id_[record_id]] = record_id;
  }
  std::vector<double> path_time;
  std::vector<int> unique_path;
  std::vector<int> parent_path;
  std::vector<int> path_to_state;
  std::vector<int> path_from_state;
  // Copy committed PUNIQ buffers in path-id order without re-evaluating them.
  for (const std::pair<const int, int>& item : appearance_index_by_path) {
    int record_id = item.second;
    path_time.push_back(path_appearance_time_[record_id]);
    unique_path.push_back(path_appearance_id_[record_id]);
    parent_path.push_back(path_appearance_parent_id_[record_id]);
    path_to_state.push_back(path_appearance_to_state_id_[record_id]);
    path_from_state.push_back(path_appearance_from_state_id_[record_id]);
  }
  Rcpp::DataFrame path_records = Rcpp::DataFrame::create(
    Rcpp::Named("time") = path_time,
    Rcpp::Named("path_id") = unique_path,
    Rcpp::Named("parent_path_id") = parent_path,
    Rcpp::Named("to_state_id") = path_to_state,
    Rcpp::Named("from_state_id") = path_from_state,
    Rcpp::Named("stringsAsFactors") = false
  );

  // Batch commit already classified LDIF records against complete post-batch
  // source support. Serialization copies that authoritative event buffer.

  std::vector<int> through_path_id;
  std::vector<std::string> through_path_label;
  std::vector<int> lineage_through;
  std::vector<int> scenario_through;
  // Retain the compatibility metadata table without an independent ledger: its
  // values are direct final CLTT/CSTT path cells.
  const V35OnlineCounterRow& final_row = committed_rows.back();
  for (int path_id : path_ids) {
    double lineage_value = v35_online_path_value(
      final_row.cumulative_lineage_by_path,
      path_id
    );
    double scenario_value = v35_online_path_value(
      final_row.cumulative_scenario_by_path,
      path_id
    );
    if (lineage_value <= V35_ONLINE_NUMERIC_TOLERANCE &&
        scenario_value <= V35_ONLINE_NUMERIC_TOLERANCE) {
      continue;
    }
    through_path_id.push_back(path_id);
    through_path_label.push_back(path_labels[path_id - 1]);
    lineage_through.push_back(static_cast<int>(std::round(lineage_value)));
    scenario_through.push_back(static_cast<int>(std::round(scenario_value)));
  }
  Rcpp::DataFrame path_through = Rcpp::DataFrame::create(
    Rcpp::Named("path_id") = through_path_id,
    Rcpp::Named("path") = through_path_label,
    Rcpp::Named("lin_through") = lineage_through,
    Rcpp::Named("scn_through") = scenario_through,
    Rcpp::Named("stringsAsFactors") = false
  );

  Rcpp::NumericVector time(committed_rows.size());
  // Copy the one-row-per-batch time grid used by every online TT family.
  for (int row = 0; row < static_cast<int>(committed_rows.size()); ++row) {
    time[row] = committed_rows[row].time;
  }

  Rcpp::DataFrame change_frame = v35_online_change_frame(changes_);
  Rcpp::DataFrame lineage_differentiating_frame =
    v35_online_change_frame(lineage_differentiating_changes_);
  std::set<std::string> raw_scenario_change_levels;
  std::set<std::string> all_change_levels;
  // Collect the complete transition vocabulary and separately identify raw
  // scenario-edge changes. The latter lead the saved LDIF column order; other
  // CTT event types remain visible afterward as stable zero-valued columns.
  for (const V35OnlineChangeRecord& change : changes_) {
    // Transition labels require both state ids to resolve in the immutable lookup.
    if (change.from_state_id < 1 ||
        change.from_state_id >= static_cast<int>(state_labels_by_id_.size()) ||
        change.to_state_id < 1 ||
        change.to_state_id >= static_cast<int>(state_labels_by_id_.size())) {
      Rcpp::stop("V35 online change record contains an invalid state id");
    }
    std::string event_label =
      state_labels_by_id_[change.from_state_id] + "->" +
      state_labels_by_id_[change.to_state_id];
    all_change_levels.insert(event_label);
    // Raw scenario changes establish the primary LDIF event-level order even
    // when post-batch source depletion suppresses their numeric score.
    if (v35_online_has_id(change.scenario_edge_before) &&
        v35_online_has_id(change.scenario_edge_after) &&
        change.scenario_edge_before != change.scenario_edge_after) {
      raw_scenario_change_levels.insert(event_label);
    }
  }
  std::vector<std::string> complete_event_levels(
    raw_scenario_change_levels.begin(),
    raw_scenario_change_levels.end()
  );
  // Append event types not represented by a raw scenario-edge change. This
  // includes virtual terminal-scenario differentiation and CTT-only events.
  for (const std::string& event_label : all_change_levels) {
    // Primary levels are already present and must not be duplicated.
    if (raw_scenario_change_levels.find(event_label) != raw_scenario_change_levels.end()) {
      continue;
    }
    complete_event_levels.push_back(event_label);
  }
  lineage_differentiating_frame.attr("event_levels") = Rcpp::wrap(
    complete_event_levels
  );

  // Exposure maps were integrated when each atomic batch committed. Formatting
  // now writes the canonical path buffer directly without replaying LTT rows.
  Rcpp::NumericMatrix exposure_path = v35_online_path_vector_matrix(
    committed_rows,
    path_ids,
    path_labels,
    &V35OnlineCounterRow::exposure_by_path
  );
  Rcpp::IntegerVector final_semantic_prefix_by_tip(ntips_);
  // Preserve each tip's final recursive `(prior prefix, path)` identity. Equal
  // ids prove equal complete committed path histories without replaying dense
  // state/path matrices during scenario-tree canonicalization.
  for (int tip_id = 1; tip_id <= ntips_; ++tip_id) {
    final_semantic_prefix_by_tip[tip_id - 1] =
      semantic_prefix_by_tip_[tip_id];
  }

  // The core retains only canonical path matrices and scored records. Public
  // state/total views are derived from each path matrix in `tt_output.cpp`.
  Rcpp::List out(19);
  out[0] = time;
  out[1] = lineage_path;
  out[2] = cumulative_lineage_path;
  out[3] = scenario_path;
  out[4] = cumulative_scenario_path;
  out[5] = sequence_path;
  out[6] = ctt_path;
  out[7] = ldif_path;
  out[8] = puniq_path;
  out[9] = change_frame;
  out[10] = lineage_differentiating_frame;
  out[11] = path_records;
  out[12] = path_through;
  out[13] = exposure_path;
  out[14] = final_semantic_prefix_by_tip;
  out[15] = "online_traversal_counters";
  out[16] = "atomic_same_time_gross_gains";
  out[17] = "atomic_batch_commit";
  out[18] = "canonical_trans_dependent_views";
  std::vector<std::string> field_names{
    "time_vec", "lin_path", "cumulative_lin_path",
    "scn_path", "cumulative_scn_path", "seq_path",
    "ctt_path", "ldif_path", "puniq_path",
    "change_records", "lineage_differentiating_records", "path_records",
    "path_through_counts", "tlen_path", "final_semantic_prefix_by_tip",
    "construction_mode", "batch_semantics", "scoring_phase",
    "finalization_mode"
  };
  out.attr("names") = Rcpp::wrap(field_names);
  return out;
}

/** Return how many provisional frontiers passed ownership validation. */
int V35LiveTTState::validated_capture_count() const {
  return validated_capture_count_;
}

/** Return how many atomic rows passed finalized counter validation. */
int V35LiveTTState::validated_commit_count() const {
  return validated_commit_count_;
}

/** Return how many committed path rows passed ledger/frontier reconciliation. */
int V35LiveTTState::validated_reconciliation_count() const {
  return validated_reconciliation_count_;
}

/** Return elapsed time owned by full-row atomic commits. */
double V35LiveTTState::batch_commit_seconds() const {
  return batch_commit_seconds_;
}
