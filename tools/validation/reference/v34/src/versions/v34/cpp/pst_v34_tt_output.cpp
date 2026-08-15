#include "pst_v34_tt_output.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

static std::vector<std::string> v34_labels_from_lookup(Rcpp::DataFrame lookup, const char* id_col, const char* label_col);

/**
 * Check integer NA sentinels from R before using ids in C++ indexing.
 *
 * Input and output: accepts one R integer and returns whether it is usable as
 * an id. It owns no state and cannot fail.
 */
static bool v34_not_na_int(int value) {
  return value != NA_INTEGER;
}

/**
 * Attach totals and time to a dense TT value matrix.
 *
 * `values` is row-major `(time, semantic level)` data already scored by
 * traversal. This function adds only public labels, totals, and the shared time
 * column; mismatched dimensions are a caller invariant.
 */
static Rcpp::NumericMatrix v34_make_tt_matrix_cpp(
    const std::vector<double>& values,
    int ntime,
    const std::vector<std::string>& labels,
    const std::vector<double>& totals,
    const Rcpp::NumericVector& time_vec) {
  int nlevels = static_cast<int>(labels.size());
  Rcpp::NumericMatrix out(ntime, nlevels + 2);
  Rcpp::CharacterVector colnames(nlevels + 2);
  // Serialize one semantic event level at a time. Each inner pass copies that
  // level across the complete batch axis without changing cumulative values.
  for (int col = 0; col < nlevels; ++col) {
    colnames[col] = labels[col];
    // Copy every committed batch value for this event level; matrix rows remain
    // aligned with `time_vec` after each iteration.
    for (int row = 0; row < ntime; ++row) {
      out(row, col) = values[row * nlevels + col];
    }
  }
  colnames[nlevels] = "total";
  colnames[nlevels + 1] = "summary_time_vec";
  // Append the precomputed total and time coordinate for each committed batch.
  for (int row = 0; row < ntime; ++row) {
    out(row, nlevels) = totals[row];
    out(row, nlevels + 1) = time_vec[row];
  }
  out.attr("dimnames") = Rcpp::List::create(R_NilValue, colnames);
  return out;
}

/**
 * Create a total-only TT matrix on the shared time axis.
 *
 * Totals were scored before this formatting boundary. The function pairs them
 * with the common batch time axis and performs no event classification.
 */
static Rcpp::NumericMatrix v34_make_total_matrix_cpp(
    const std::vector<double>& totals,
    const Rcpp::NumericVector& time_vec) {
  int ntime = time_vec.size();
  Rcpp::NumericMatrix out(ntime, 2);
  Rcpp::CharacterVector colnames = Rcpp::CharacterVector::create("total", "summary_time_vec");
  // Preserve one total/time pair per committed batch in chronological order.
  for (int row = 0; row < ntime; ++row) {
    out(row, 0) = totals[row];
    out(row, 1) = time_vec[row];
  }
  out.attr("dimnames") = Rcpp::List::create(R_NilValue, colnames);
  return out;
}

/**
 * Translate integer ids to public labels through an R lookup table.
 *
 * The lookup owns 1-based traversal ids. The returned direct-index vector keeps
 * slot zero as a sentinel and copies labels without changing their order.
 */
static std::vector<std::string> v34_labels_from_lookup(Rcpp::DataFrame lookup, const char* id_col, const char* label_col) {
  Rcpp::IntegerVector ids = lookup[id_col];
  Rcpp::CharacterVector labels = lookup[label_col];
  int max_id = 0;
  // Determine the direct-index vector size from every nonmissing public id.
  for (int id : ids) {
    // Missing lookup rows carry no label and cannot affect storage size.
    if (v34_not_na_int(id)) {
      max_id = std::max(max_id, id);
    }
  }
  std::vector<std::string> out(max_id + 1);
  // Install each lookup label at its canonical id; after every iteration all
  // processed ids resolve in O(1) time during event formatting.
  for (int i = 0; i < ids.size(); ++i) {
    int id = ids[i];
    // Ignore missing or negative sentinels because they are not public ids.
    if (v34_not_na_int(id) && id >= 0) {
      out[id] = Rcpp::as<std::string>(labels[i]);
    }
  }
  return out;
}

/**
 * Translate integer ids to public labels through an R lookup table.
 *
 * Unlike direct id lookup, this preserves data-frame row order for public TT
 * columns. It performs no sorting or deduplication.
 */
static std::vector<std::string> v34_label_levels_from_lookup(Rcpp::DataFrame lookup, const char* label_col) {
  Rcpp::CharacterVector labels = lookup[label_col];
  std::vector<std::string> out(labels.size());
  // Copy each declared level once so all dependent views share column order.
  for (int i = 0; i < labels.size(); ++i) {
    out[i] = Rcpp::as<std::string>(labels[i]);
  }
  return out;
}


}  // namespace

/**
 * Derive TT total and state views from one canonical path matrix.
 *
 * `path_matrix` is the sole numeric source. `path_lookup` maps each escaped
 * public path label to its terminal state id, and `state_lookup` defines output
 * column order. The function performs deterministic aggregation only; it owns
 * no biological counters and cannot change path values.
 */
Rcpp::List pst_v34_path_counter_views_cpp(
    Rcpp::NumericMatrix path_matrix,
    Rcpp::DataFrame path_lookup,
    Rcpp::DataFrame state_lookup) {
  Rcpp::CharacterVector path_columns = Rcpp::colnames(path_matrix);
  int total_column = -1;
  int time_column = -1;
  // Locate the two nonsemantic columns once before aggregating path values.
  for (int column = 0; column < path_columns.size(); ++column) {
    std::string label = Rcpp::as<std::string>(path_columns[column]);
    // The total column is copied verbatim into dependent public views.
    if (label == "total") {
      total_column = column;
    } else if (label == "summary_time_vec") {
      // The shared batch time column remains unchanged across all views.
      time_column = column;
    }
  }
  // A canonical path matrix always carries both formatting columns.
  if (total_column < 0 || time_column < 0) {
    Rcpp::stop("V34 canonical path matrix omits total or time");
  }

  Rcpp::IntegerVector state_ids = state_lookup["state_id"];
  Rcpp::CharacterVector state_labels = state_lookup["state"];
  std::vector<std::string> state_levels(state_labels.size());
  std::map<int, int> state_index_by_id;
  // Preserve public state order while indexing traversal ids for aggregation.
  for (int state_index = 0; state_index < state_labels.size(); ++state_index) {
    state_levels[state_index] = Rcpp::as<std::string>(state_labels[state_index]);
    state_index_by_id[state_ids[state_index]] = state_index;
  }

  Rcpp::CharacterVector lookup_labels = path_lookup["label"];
  Rcpp::IntegerVector terminal_state_ids = path_lookup["terminal_state_id"];
  std::map<std::string, int> terminal_state_by_path;
  // Index each canonical path label by its traversal-owned terminal state id.
  for (int path_index = 0; path_index < lookup_labels.size(); ++path_index) {
    // Synthetic placeholders have no state and cannot receive numeric support.
    if (terminal_state_ids[path_index] == NA_INTEGER) {
      continue;
    }
    terminal_state_by_path[
      Rcpp::as<std::string>(lookup_labels[path_index])
    ] = terminal_state_ids[path_index];
  }

  std::vector<int> state_index_by_column(path_columns.size(), -1);
  // Path-column loop:
  //   Resolve every semantic path label to its public terminal-state column
  //   once. This lookup depends only on matrix columns, not time rows, and must
  //   not be repeated across the potentially large `(time x path)` surface.
  for (int column = 0; column < path_columns.size(); ++column) {
    // Total and time are dependent formatting columns and have no terminal state.
    if (column == total_column || column == time_column) {
      continue;
    }
    std::string path_label = Rcpp::as<std::string>(path_columns[column]);
    std::map<std::string, int>::const_iterator terminal_state =
      terminal_state_by_path.find(path_label);
    // Every numeric path column must resolve through the canonical lookup.
    if (terminal_state == terminal_state_by_path.end()) {
      Rcpp::stop("V34 path counter column is absent from path lookup");
    }
    std::map<int, int>::const_iterator state_index =
      state_index_by_id.find(terminal_state->second);
    // A path terminal state outside the declared state domain is malformed.
    if (state_index == state_index_by_id.end()) {
      Rcpp::stop("V34 path counter terminal state is absent from state lookup");
    }
    state_index_by_column[column] = state_index->second;
  }

  int ntime = path_matrix.nrow();
  int nstate = state_levels.size();
  std::vector<double> state_values(
    static_cast<std::size_t>(ntime) * nstate,
    0.0
  );
  std::vector<double> totals(ntime, 0.0);
  Rcpp::NumericVector time_vec(ntime);
  // Recalculate totals from semantic path columns and aggregate those same
  // values by terminal state. The serialized `total` column is deliberately
  // not an independent numerical source.
  for (int row = 0; row < ntime; ++row) {
    totals[row] = 0.0;
    time_vec[row] = path_matrix(row, time_column);
    // Path-column loop:
    //   Aggregate numeric support through the precomputed column-to-state map.
    //   After each row, state cells represent the same canonical path values
    //   without repeating label conversion or associative lookup work.
    for (int column = 0; column < path_columns.size(); ++column) {
      int state_index = state_index_by_column[column];
      // A negative precomputed index marks total/time rather than a path value.
      if (state_index < 0) {
        continue;
      }
      double value = path_matrix(row, column);
      totals[row] += value;
      state_values[static_cast<std::size_t>(row) * nstate + state_index] +=
        value;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("tot") = v34_make_total_matrix_cpp(totals, time_vec),
    Rcpp::Named("state") = v34_make_tt_matrix_cpp(
      state_values,
      ntime,
      state_levels,
      totals,
      time_vec
    ),
    Rcpp::Named("trans") = path_matrix
  );
}

/**
 * Build cumulative transition-event TT views on the supplied shared time axis.
 *
 * `canonical_trans` contains traversal-owned CTT, LDIF, or PUNIQ rows. Event
 * records contribute only the ordered legacy event vocabulary; no biological
 * value is recounted from them. All numerical outputs are deterministic views
 * of the supplied canonical path matrix.
 */
Rcpp::List pst_v34_transition_views_cpp(
    Rcpp::DataFrame records,
    Rcpp::NumericMatrix canonical_trans,
    Rcpp::List summary,
    Rcpp::NumericVector time_vec,
    double time_tolerance) {
  Rcpp::DataFrame state_lookup = Rcpp::as<Rcpp::DataFrame>(summary["state_lookup"]);
  Rcpp::DataFrame path_lookup = Rcpp::as<Rcpp::DataFrame>(summary["path_lookup"]);
  std::vector<std::string> state_by_id = v34_labels_from_lookup(state_lookup, "state_id", "state");
  std::vector<std::string> state_names = v34_label_levels_from_lookup(state_lookup, "state");
  Rcpp::IntegerVector lookup_path_ids = path_lookup["path_id"];
  Rcpp::IntegerVector lookup_parent_path_ids = path_lookup["parent_path_id"];
  Rcpp::IntegerVector lookup_terminal_state_ids = path_lookup["terminal_state_id"];
  Rcpp::CharacterVector lookup_path_labels = path_lookup["label"];
  int ntime = time_vec.size();

  std::vector<std::string> event_levels;
  std::unordered_set<std::string> event_level_seen;
  SEXP event_level_attr = records.attr("event_levels");
  // Preserve an explicitly supplied event-level order. LDIF uses this to keep
  // raw scenario-changing events first while retaining zero-valued CTT events.
  if (!Rf_isNull(event_level_attr)) {
    Rcpp::CharacterVector attr_levels(event_level_attr);
    // Copy each declared level once and in order; later records may append only
    // event types absent from this explicit vocabulary.
    for (int i = 0; i < attr_levels.size(); ++i) {
      // Missing labels do not represent a public event column.
      if (attr_levels[i] == NA_STRING) {
        continue;
      }
      std::string level = Rcpp::as<std::string>(attr_levels[i]);
      // Duplicate labels would create ambiguous TT columns, so retain only the first.
      if (event_level_seen.insert(level).second) {
        event_levels.push_back(level);
      }
    }
  }

  Rcpp::IntegerVector from_state_id = records["from_state_id"];
  Rcpp::IntegerVector to_state_id = records["to_state_id"];
  int nrecord = records.nrows();
  std::vector<std::string> event_labels(nrecord);
  // Translate every committed record into the labels required by each
  // dependent view. No record is reclassified as CTT, LDIF, or PUNIQ here.
  for (int i = 0; i < nrecord; ++i) {
    std::string from = from_state_id[i] >= 0 && from_state_id[i] < static_cast<int>(state_by_id.size()) ?
      state_by_id[from_state_id[i]] : std::string("");
    std::string to = to_state_id[i] >= 0 && to_state_id[i] < static_cast<int>(state_by_id.size()) ?
      state_by_id[to_state_id[i]] : std::string("");
    event_labels[i] = from + "->" + to;
    // Append record-derived event types only when they were not declared by
    // the ordered attribute or an earlier record.
    if (event_level_seen.insert(event_labels[i]).second) {
      event_levels.push_back(event_labels[i]);
    }
  }
  // Families without an explicit ordering attribute retain the historical
  // lexicographic event-column contract used by CTT and PUNIQ references.
  if (Rf_isNull(event_level_attr)) {
    std::sort(event_levels.begin(), event_levels.end());
  }

  if (canonical_trans.nrow() != ntime) {
    Rcpp::stop("V34 transition trans rows do not match the shared time axis");
  }
  Rcpp::List canonical_views = pst_v34_path_counter_views_cpp(
    canonical_trans,
    path_lookup,
    state_lookup
  );
  Rcpp::NumericMatrix total_matrix = canonical_views["tot"];
  Rcpp::NumericMatrix state_matrix = canonical_views["state"];
  Rcpp::NumericMatrix path_matrix = canonical_views["trans"];
  std::vector<double> cumulative_totals(ntime, 0.0);
  // The total and time views are deterministic projections of canonical trans.
  for (int row = 0; row < ntime; ++row) {
    cumulative_totals[row] = total_matrix(row, 0);
    if (std::fabs(total_matrix(row, 1) - time_vec[row]) > time_tolerance) {
      Rcpp::stop("V34 transition trans time does not match the shared time axis");
    }
  }

  std::map<int, int> terminal_state_by_path_id;
  std::map<std::string, int> path_id_by_label;
  // Index canonical path ancestry and terminal-state identity once. These
  // lookups are the only semantic inputs used to derive non-path views.
  for (int lookup_row = 0; lookup_row < lookup_path_ids.size(); ++lookup_row) {
    int canonical_path_id = lookup_path_ids[lookup_row];
    terminal_state_by_path_id[canonical_path_id] =
      lookup_terminal_state_ids[lookup_row];
    path_id_by_label[
      Rcpp::as<std::string>(lookup_path_labels[lookup_row])
    ] = canonical_path_id;
  }
  std::map<int, int> state_output_index;
  Rcpp::IntegerVector output_state_ids = state_lookup["state_id"];
  // Preserve the declared public state-column order during path aggregation.
  for (int state_index = 0; state_index < output_state_ids.size(); ++state_index) {
    state_output_index[output_state_ids[state_index]] = state_index;
  }
  std::map<std::string, int> event_output_index;
  // Cache the established event-level order, including explicit zero columns.
  for (int event_index = 0;
       event_index < static_cast<int>(event_levels.size());
       ++event_index) {
    event_output_index[event_levels[event_index]] = event_index;
  }

  int nstate = state_names.size();
  int nevent = event_levels.size();
  std::vector<double> leaving_values(
    static_cast<std::size_t>(ntime) * nstate,
    0.0
  );
  std::vector<double> event_values(
    static_cast<std::size_t>(ntime) * nevent,
    0.0
  );
  Rcpp::CharacterVector path_columns = Rcpp::colnames(path_matrix);
  // Aggregate canonical destination-path counts into legacy event and leaving
  // views. State and arrival are already the same terminal-state aggregation.
  for (int path_column = 0; path_column < path_columns.size(); ++path_column) {
    std::string path_label = Rcpp::as<std::string>(path_columns[path_column]);
    // Total and time are dependent formatting columns, not transition paths.
    if (path_label == "total" || path_label == "summary_time_vec") {
      continue;
    }
    std::map<std::string, int>::const_iterator path_id_entry =
      path_id_by_label.find(path_label);
    // Every canonical path column must resolve through the path lookup.
    if (path_id_entry == path_id_by_label.end()) {
      Rcpp::stop("V34 transition path column is absent from path lookup");
    }
    int destination_path_id = path_id_entry->second;
    int destination_state_id = terminal_state_by_path_id[destination_path_id];
    int lookup_row = destination_path_id - 1;
    // Path ids are contiguous and row-aligned by construction; an invalid row
    // indicates that public lookup serialization lost canonical ordering.
    if (lookup_row < 0 || lookup_row >= lookup_parent_path_ids.size()) {
      Rcpp::stop("V34 transition path id is outside parent lookup");
    }
    int source_path_id = lookup_parent_path_ids[lookup_row];
    int source_state_id = source_path_id == NA_INTEGER ? NA_INTEGER :
      terminal_state_by_path_id[source_path_id];
    std::map<int, int>::const_iterator destination_state_index =
      state_output_index.find(destination_state_id);
    // A scored destination path must terminate in one declared biological state.
    if (destination_state_index == state_output_index.end()) {
      Rcpp::stop("V34 transition destination state is absent from lookup");
    }
    std::string from = source_state_id == NA_INTEGER ? std::string("") :
      state_by_id[source_state_id];
    std::string to = state_by_id[destination_state_id];
    std::string event_label = from + "->" + to;
    std::map<std::string, int>::const_iterator event_index =
      event_output_index.find(event_label);
    std::map<int, int>::const_iterator source_state_index =
      state_output_index.find(source_state_id);
    // Copy this cumulative path count into every deterministic dependent view.
    for (int row = 0; row < ntime; ++row) {
      double value = path_matrix(row, path_column);
      // Cross-state transitions contribute leaving support to their source;
      // same-state path arrivals intentionally do not.
      if (source_state_id != NA_INTEGER &&
          source_state_id != destination_state_id &&
          source_state_index != state_output_index.end()) {
        leaving_values[static_cast<std::size_t>(row) * nstate +
                       source_state_index->second] += value;
      }
      // Explicit event vocabularies may retain zero columns not represented by
      // a path; represented path events add to their one canonical column.
      if (event_index != event_output_index.end()) {
        event_values[static_cast<std::size_t>(row) * nevent +
                     event_index->second] += value;
      }
    }
  }

  Rcpp::NumericMatrix event_matrix = v34_make_tt_matrix_cpp(
    event_values, ntime, event_levels, cumulative_totals, time_vec
  );
  Rcpp::NumericMatrix arr_matrix = Rcpp::clone(state_matrix);
  Rcpp::NumericMatrix leave_matrix = v34_make_tt_matrix_cpp(
    leaving_values, ntime, state_names, cumulative_totals, time_vec
  );

  return Rcpp::List::create(
    Rcpp::Named("tot") = total_matrix,
    Rcpp::Named("state") = state_matrix,
    Rcpp::Named("event") = event_matrix,
    Rcpp::Named("arr") = arr_matrix,
    Rcpp::Named("leave") = leave_matrix,
    Rcpp::Named("trans") = path_matrix
  );
}
