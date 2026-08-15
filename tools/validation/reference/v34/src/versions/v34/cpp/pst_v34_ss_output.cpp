#include "pst_v34_ss_output.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct V34TreeView {
  Rcpp::IntegerMatrix edge;
  Rcpp::NumericVector edge_length;
  Rcpp::List maps;
  Rcpp::NumericMatrix mapped_edge;
  Rcpp::CharacterVector tip_label;
  std::vector<std::string> state_levels;
};

struct V34TreeSummaryParts {
  Rcpp::List tipstats;
  Rcpp::List transtats;
  Rcpp::NumericVector treelength_tree;
  Rcpp::List treelength_bystate;
};

struct V34PathSizeSummary {
  std::map<std::string, double> lineage_through;
  std::map<std::string, double> scenario_through;
  std::map<std::string, double> lineage_stopped;
  std::map<std::string, double> scenario_stopped;
};

/**
 * Convert one R character scalar to its C++ label representation.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static std::string v34_as_string(SEXP value) {
  return Rcpp::as<std::string>(value);
}

/**
 * Extract mapped-segment names in their original edge-local order.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static std::vector<std::string> v34_character_names(SEXP x) {
  Rcpp::CharacterVector names = Rf_getAttrib(x, R_NamesSymbol);
  std::vector<std::string> out;
  out.reserve(names.size());
  // Copy map labels in segment order so later SS counts remain paired with their durations.
  for (int i = 0; i < names.size(); ++i) {
    out.push_back(v34_as_string(names[i]));
  }
  return out;
}

/**
 * Read matrix column labels used as the canonical SS state domain.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static std::vector<std::string> v34_dimnames_colnames(SEXP x) {
  Rcpp::List dimnames = Rf_getAttrib(x, R_DimNamesSymbol);
  // A matrix without column names defines no state/path domain for dependent SS formatting.
  if (dimnames.size() < 2 || Rf_isNull(dimnames[1])) {
    return std::vector<std::string>();
  }
  Rcpp::CharacterVector names = dimnames[1];
  std::vector<std::string> out;
  out.reserve(names.size());
  // Preserve the declared matrix column order while moving labels into owned C++ storage.
  for (int i = 0; i < names.size(); ++i) {
    out.push_back(v34_as_string(names[i]));
  }
  return out;
}

/**
 * Sort and deduplicate public labels for deterministic SS row ordering.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static std::vector<std::string> v34_sorted_unique(std::vector<std::string> values) {
  std::sort(values.begin(), values.end());
  values.erase(std::unique(values.begin(), values.end()), values.end());
  return values;
}

/**
 * Project an R tree list into a compact C++ view used by summary-statistics helpers.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static V34TreeView v34_tree_view(Rcpp::List tree) {
  V34TreeView out;
  out.edge = Rcpp::as<Rcpp::IntegerMatrix>(tree["edge"]);
  out.edge_length = Rcpp::as<Rcpp::NumericVector>(tree["edge.length"]);
  out.maps = Rcpp::as<Rcpp::List>(tree["maps"]);
  out.mapped_edge = Rcpp::as<Rcpp::NumericMatrix>(tree["mapped.edge"]);
  out.tip_label = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]);
  out.state_levels = v34_sorted_unique(v34_dimnames_colnames(out.mapped_edge));
  return out;
}

/**
 * Parse or summarize mapped state labels for tree and TT statistics.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static std::string v34_first_map_state(SEXP map_sexp) {
  std::vector<std::string> names = v34_character_names(map_sexp);
  return names.empty() ? std::string() : names.front();
}

/**
 * Parse or summarize mapped state labels for tree and TT statistics.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static std::string v34_last_map_state(SEXP map_sexp) {
  std::vector<std::string> names = v34_character_names(map_sexp);
  return names.empty() ? std::string() : names.back();
}

/**
 * Split an encoded path or label string into stable components.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static std::vector<std::string> v34_split_path(const std::string& path) {
  std::vector<std::string> out;
  std::string current;
  // Split each canonical path delimiter while retaining empty root-prefix components.
  for (char ch : path) {
    // A delimiter closes the current path step, including the leading empty root marker.
    if (ch == '|') {
      out.push_back(current);
      current.clear();
    } else {
      // Non-delimiter bytes belong to the current encoded state step.
      current.push_back(ch);
    }
  }
  out.push_back(current);
  return out;
}

/**
 * Return the terminal state component of a canonical transition path.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static std::string v34_last_step(const std::string& path) {
  std::vector<std::string> pieces = v34_split_path(path);
  return pieces.empty() ? std::string() : pieces.back();
}

/**
 * Count state components in a canonical transition path.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static int v34_count_steps(const std::string& path) {
  return static_cast<int>(v34_split_path(path).size());
}

/**
 * Resolve state labels to stable rank/index positions for summary tables.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static int v34_state_rank(
    const std::string& state,
    const std::vector<std::string>& state_levels) {
  // Find the state's public order rank; unknown states sort after the declared domain.
  for (int i = 0; i < static_cast<int>(state_levels.size()); ++i) {
    // Exact label equality identifies the state's stable summary-table position.
    if (state_levels[i] == state) {
      return i;
    }
  }
  return static_cast<int>(state_levels.size());
}

/**
 * Create a named integer table with stable columns.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static Rcpp::IntegerVector v34_named_table(
    const std::vector<std::string>& names,
    const std::vector<int>& values) {
  Rcpp::IntegerVector out(values.size());
  Rcpp::CharacterVector out_names(values.size());
  // Copy aligned names/counts into the legacy-compatible named `table` vector.
  for (int i = 0; i < static_cast<int>(values.size()); ++i) {
    out[i] = values[i];
    out_names[i] = names[i];
  }
  out.attr("names") = out_names;
  out.attr("class") = "table";
  return out;
}

/**
 * Create a named numeric vector with stable output labels.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static Rcpp::NumericVector v34_named_numeric(
    const std::vector<std::string>& names,
    const std::vector<double>& values) {
  Rcpp::NumericVector out(values.size());
  Rcpp::CharacterVector out_names(values.size());
  // Copy aligned labels and values into one deterministic numeric SS vector.
  for (int i = 0; i < static_cast<int>(values.size()); ++i) {
    out[i] = values[i];
    out_names[i] = names[i];
  }
  out.attr("names") = out_names;
  return out;
}

/**
 * Read or build a matrix-shaped value while preserving V34 row/time alignment.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static Rcpp::NumericMatrix v34_zero_matrix(const std::vector<std::string>& states) {
  int n = static_cast<int>(states.size());
  Rcpp::NumericMatrix out(n, n);
  Rcpp::CharacterVector dim_names(n);
  // Use the same ordered state labels on transition-matrix rows and columns.
  for (int i = 0; i < n; ++i) {
    dim_names[i] = states[i];
  }
  out.attr("dimnames") = Rcpp::List::create(dim_names, Rcpp::clone(dim_names));
  return out;
}

/**
 * Resolve state labels to stable rank/index positions for summary tables.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static int v34_state_index(
    const std::vector<std::string>& states,
    const std::string& state) {
  // Resolve one state label against the exact matrix domain used by this SS surface.
  for (int i = 0; i < static_cast<int>(states.size()); ++i) {
    // Exact equality returns the corresponding row/column index.
    if (states[i] == state) {
      return i;
    }
  }
  return -1;
}

/**
 * Build state-transition count matrices from a tree view.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static Rcpp::NumericMatrix v34_transition_matrix(
    const V34TreeView& tree,
    const std::vector<std::string>& states,
    bool include_edge_boundaries) {
  Rcpp::NumericMatrix out = v34_zero_matrix(states);

  // Count state changes inside each mapped edge; adjacent equal labels are continuations.
  for (int edge_i = 0; edge_i < tree.maps.size(); ++edge_i) {
    std::vector<std::string> labels = v34_character_names(tree.maps[edge_i]);
    // Compare consecutive map segments in biological time order.
    for (int j = 1; j < static_cast<int>(labels.size()); ++j) {
      // Equal adjacent states do not define a transition event.
      if (labels[j - 1] == labels[j]) {
        continue;
      }
      int from = v34_state_index(states, labels[j - 1]);
      int to = v34_state_index(states, labels[j]);
      // Count only transitions whose endpoints belong to the declared state domain.
      if (from >= 0 && to >= 0) {
        out(from, to) += 1.0;
      }
    }
  }

  // Scenario and transition trees may encode a state change across an edge boundary;
  // phylogeny maps do not request this extra topology-level comparison.
  if (include_edge_boundaries) {
    std::map<int, int> edge_by_child;
    // Index each incoming edge by child node for parent-terminal-state lookup.
    for (int i = 0; i < tree.edge.nrow(); ++i) {
      edge_by_child[tree.edge(i, 1)] = i;
    }
    // Compare every nonroot edge's first state with its parent edge's terminal state.
    for (int i = 0; i < tree.edge.nrow(); ++i) {
      int parent_node = tree.edge(i, 0);
      auto parent_edge = edge_by_child.find(parent_node);
      // Root edges have no incoming parent state and therefore no boundary transition.
      if (parent_edge == edge_by_child.end()) {
        continue;
      }
      std::string from_state = v34_last_map_state(tree.maps[parent_edge->second]);
      std::string to_state = v34_first_map_state(tree.maps[i]);
      // Empty labels and state continuations do not contribute a biological transition.
      if (from_state.empty() || to_state.empty() || from_state == to_state) {
        continue;
      }
      int from = v34_state_index(states, from_state);
      int to = v34_state_index(states, to_state);
      // Count only boundary changes represented in the requested state domain.
      if (from >= 0 && to >= 0) {
        out(from, to) += 1.0;
      }
    }
  }

  return out;
}

/**
 * Parse or summarize mapped state labels for tree and TT statistics.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static Rcpp::IntegerVector v34_ntips_bystate(const V34TreeView& tree) {
  std::map<std::string, int> counts;
  int ntips = tree.tip_label.size();
  // Count the terminal mapped state of every edge ending at a public tree tip.
  for (int i = 0; i < tree.edge.nrow(); ++i) {
    // Only terminal child nodes contribute tip-state diversity.
    if (tree.edge(i, 1) <= ntips) {
      std::string state = v34_last_map_state(tree.maps[i]);
      counts[state] += 1;
    }
  }

  std::vector<std::string> names;
  std::vector<int> values;
  // Serialize lexical state/count pairs into the expected named table.
  for (const auto& item : counts) {
    names.push_back(item.first);
    values.push_back(item.second);
  }
  return v34_named_table(names, values);
}

/**
 * Parse, count, or summarize path-history labels for transition/path statistics.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static std::vector<int> v34_node_path(const V34TreeView& tree, int tip_node) {
  std::vector<int> reversed;
  int node = tip_node;
  reversed.push_back(node);
  // Follow incoming edges from the tip to the root, collecting node ids for later reversal.
  while (true) {
    int parent = NA_INTEGER;
    // Locate the unique incoming edge for the current node.
    for (int i = 0; i < tree.edge.nrow(); ++i) {
      // A matching child endpoint identifies the next ancestor.
      if (tree.edge(i, 1) == node) {
        parent = tree.edge(i, 0);
        break;
      }
    }
    // Absence of an incoming edge means the root has been reached.
    if (parent == NA_INTEGER) {
      break;
    }
    reversed.push_back(parent);
    node = parent;
  }
  std::reverse(reversed.begin(), reversed.end());
  return reversed;
}

/**
 * Count represented tips by transition or state labels.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static Rcpp::IntegerVector v34_ntips_bytrans(
    const V34TreeView& tree,
    const std::vector<std::string>& state_levels) {
  std::map<int, int> edge_by_child;
  // Index each tree edge by child node for root-to-tip map reconstruction.
  for (int i = 0; i < tree.edge.nrow(); ++i) {
    edge_by_child[tree.edge(i, 1)] = i;
  }

  std::map<std::string, int> counts;
  int ntips = tree.tip_label.size();
  // Reconstruct and count the complete state-transition path for every public tip.
  for (int tip = 1; tip <= ntips; ++tip) {
    std::vector<int> path = v34_node_path(tree, tip);
    std::vector<std::string> states;
    // Visit root-to-tip nodes and append their incoming edge maps chronologically.
    for (int node : path) {
      auto edge_it = edge_by_child.find(node);
      // The root node has no incoming edge and contributes no mapped state segment.
      if (edge_it == edge_by_child.end()) {
        continue;
      }
      std::vector<std::string> labels = v34_character_names(tree.maps[edge_it->second]);
      // Append mapped labels while collapsing state continuations across segment and edge boundaries.
      for (const std::string& label : labels) {
        // Add a step only when the state differs from the current path tail.
        if (states.empty() || states.back() != label) {
          states.push_back(label);
        }
      }
    }
    std::ostringstream key;
    key << "|";
    // Encode the collapsed state history as the canonical `|STATE|...` path key.
    for (int i = 0; i < static_cast<int>(states.size()); ++i) {
      // Delimit every state after the first path step.
      if (i > 0) {
        key << "|";
      }
      key << states[i];
    }
    counts[key.str()] += 1;
  }

  std::vector<std::string> names;
  // Collect path labels before applying the legacy-compatible semantic sort order.
  for (const auto& item : counts) {
    names.push_back(item.first);
  }
  std::stable_sort(names.begin(), names.end(), [&](const std::string& left, const std::string& right) {
    int left_rank = v34_state_rank(v34_last_step(left), state_levels);
    int right_rank = v34_state_rank(v34_last_step(right), state_levels);
    // Primary ordering groups paths by their terminal state's declared rank.
    if (left_rank != right_rank) {
      return left_rank < right_rank;
    }
    int left_count = v34_count_steps(left);
    int right_count = v34_count_steps(right);
    // Within a terminal state, shorter transition sequences precede longer ones.
    if (left_count != right_count) {
      return left_count < right_count;
    }
    return left < right;
  });

  std::vector<int> values;
  values.reserve(names.size());
  // Emit counts in the sorted path order established above.
  for (const std::string& name : names) {
    values.push_back(counts[name]);
  }
  return v34_named_table(names, values);
}

/**
 * Group named TT/vector fields by suffix or substring for summary tables.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static Rcpp::List v34_group_named_vector_by_suffix(
    Rcpp::IntegerVector vec,
    Rcpp::CharacterVector suffixes) {
  Rcpp::CharacterVector vec_names = vec.attr("names");
  Rcpp::List out(suffixes.size());
  Rcpp::CharacterVector out_names(suffixes.size());

  // Build one grouped named vector for every requested terminal-state suffix.
  for (int s = 0; s < suffixes.size(); ++s) {
    std::string suffix = v34_as_string(suffixes[s]);
    std::vector<std::string> names;
    std::vector<int> values;
    // Select source entries whose names end in this exact suffix.
    for (int i = 0; i < vec.size(); ++i) {
      std::string name = v34_as_string(vec_names[i]);
      // Exact terminal suffix equality assigns this path count to the requested state group.
      if (name.size() >= suffix.size() &&
          name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0) {
        names.push_back(name);
        values.push_back(vec[i]);
      }
    }
    out[s] = v34_named_table(names, values);
    out_names[s] = suffix;
  }
  out.attr("names") = out_names;
  return out;
}

/**
 * Build tip-count summary-statistics tables from public trees.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static Rcpp::List v34_tip_stats(
    const V34TreeView& tree,
    const std::string& suffix,
    const std::vector<std::string>& state_levels) {
  Rcpp::IntegerVector ntips_bystate = v34_ntips_bystate(tree);
  Rcpp::IntegerVector ntips_bytrans = v34_ntips_bytrans(tree, state_levels);
  Rcpp::CharacterVector bystate_names = ntips_bystate.attr("names");
  Rcpp::List bytrans_bystate = v34_group_named_vector_by_suffix(ntips_bytrans, bystate_names);

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("ntips") = static_cast<int>(tree.tip_label.size()),
    Rcpp::Named("ntips.bystate") = ntips_bystate,
    Rcpp::Named("ntips.bytrans") = ntips_bytrans,
    Rcpp::Named("ntips.bytrans.bystate") = bytrans_bystate
  );
  Rcpp::CharacterVector names = out.attr("names");
  // Add the tree suffix to every tip-stat field without changing its value or nested shape.
  for (int i = 0; i < names.size(); ++i) {
    names[i] = v34_as_string(names[i]) + "." + suffix;
  }
  out.attr("names") = names;
  return out;
}

/**
 * Compute matrix row summaries while preserving R numeric-vector shape.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static Rcpp::NumericVector v34_matrix_row_sums(Rcpp::NumericMatrix mat) {
  Rcpp::CharacterVector row_names = Rcpp::as<Rcpp::List>(mat.attr("dimnames"))[0];
  Rcpp::NumericVector out(mat.nrow());
  // Sum each source-state row across destination states.
  for (int r = 0; r < mat.nrow(); ++r) {
    double total = 0.0;
    // Accumulate every destination count for this source state.
    for (int c = 0; c < mat.ncol(); ++c) {
      total += mat(r, c);
    }
    out[r] = total;
  }
  out.attr("names") = row_names;
  return out;
}

/**
 * Build a stable matrix-column lookup used by C++ TT and summary-statistics code.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static Rcpp::NumericVector v34_matrix_col_sums(Rcpp::NumericMatrix mat) {
  Rcpp::CharacterVector col_names = Rcpp::as<Rcpp::List>(mat.attr("dimnames"))[1];
  Rcpp::NumericVector out(mat.ncol());
  // Sum each destination-state column across source states.
  for (int c = 0; c < mat.ncol(); ++c) {
    double total = 0.0;
    // Accumulate every source count arriving at this destination state.
    for (int r = 0; r < mat.nrow(); ++r) {
      total += mat(r, c);
    }
    out[c] = total;
  }
  out.attr("names") = col_names;
  return out;
}

/**
 * Extract or transform event labels for V34 summary-statistics tables.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static Rcpp::NumericVector v34_event_vector(Rcpp::NumericMatrix mat) {
  Rcpp::List dimnames = mat.attr("dimnames");
  Rcpp::CharacterVector row_names = dimnames[0];
  Rcpp::CharacterVector col_names = dimnames[1];
  int n = mat.nrow() * mat.ncol();
  Rcpp::NumericVector out(n);
  Rcpp::CharacterVector names(n);
  int pos = 0;
  // Flatten in column-major R order so each event name remains paired with its matrix cell.
  for (int c = 0; c < mat.ncol(); ++c) {
    // Emit every source state under the current destination state.
    for (int r = 0; r < mat.nrow(); ++r) {
      out[pos] = mat(r, c);
      names[pos] = v34_as_string(row_names[r]) + "|" + v34_as_string(col_names[c]);
      ++pos;
    }
  }
  out.attr("names") = names;
  return out;
}

/**
 * Group named TT/vector fields by suffix or substring for summary tables.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static Rcpp::List v34_group_numeric_by_contains(
    Rcpp::NumericVector vec,
    Rcpp::CharacterVector patterns,
    bool by_to_state) {
  Rcpp::CharacterVector vec_names = vec.attr("names");
  Rcpp::List out(patterns.size());
  Rcpp::CharacterVector out_names(patterns.size());
  // Build one grouped event vector for every requested source or destination state.
  for (int p = 0; p < patterns.size(); ++p) {
    std::string state = v34_as_string(patterns[p]);
    std::string pattern = by_to_state ? "|" + state : state + "|";
    std::vector<std::string> names;
    std::vector<double> values;
    // Select event entries containing the legacy directional state pattern.
    for (int i = 0; i < vec.size(); ++i) {
      std::string name = v34_as_string(vec_names[i]);
      // A matching event contributes to this state's grouped SS view.
      if (name.find(pattern) != std::string::npos) {
        names.push_back(name);
        values.push_back(vec[i]);
      }
    }
    out[p] = v34_named_numeric(names, values);
    out_names[p] = pattern;
  }
  out.attr("names") = out_names;
  return out;
}

/**
 * Build transition-count summary-statistics tables from public trees.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static Rcpp::List v34_trans_stats(
    Rcpp::NumericMatrix transition_matrix,
    const std::string& suffix) {
  Rcpp::NumericVector by_to_state = v34_matrix_col_sums(transition_matrix);
  Rcpp::NumericVector by_from_state = v34_matrix_row_sums(transition_matrix);
  Rcpp::NumericVector by_event = v34_event_vector(transition_matrix);
  Rcpp::CharacterVector state_names = by_to_state.attr("names");

  double total = 0.0;
  // Sum the canonical event vector to recover total transitions for this tree projection.
  for (double value : by_event) {
    total += value;
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("ntrans.total") = total,
    Rcpp::Named("ntrans.bytostate") = by_to_state,
    Rcpp::Named("ntrans.byfromstate") = by_from_state,
    Rcpp::Named("ntrans.byevent") = by_event,
    Rcpp::Named("ntrans.byevent.bytostate") = v34_group_numeric_by_contains(by_event, state_names, true),
    Rcpp::Named("ntrans.byevent.byfromstate") = v34_group_numeric_by_contains(by_event, state_names, false)
  );
  Rcpp::CharacterVector names = out.attr("names");
  // Add the projection suffix to every transition-stat field without changing values.
  for (int i = 0; i < names.size(); ++i) {
    names[i] = v34_as_string(names[i]) + "." + suffix;
  }
  out.attr("names") = names;
  return out;
}

/**
 * Concatenate named R lists while preserving names used by public summaries.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static Rcpp::List v34_concat_named_lists(
    const Rcpp::List& first,
    const Rcpp::List& second,
    const Rcpp::List& third,
    const std::string& first_prefix,
    const std::string& second_prefix,
    const std::string& third_prefix) {
  Rcpp::CharacterVector first_names = first.attr("names");
  Rcpp::CharacterVector second_names = second.attr("names");
  Rcpp::CharacterVector third_names = third.attr("names");
  int n = first.size() + second.size() + third.size();
  Rcpp::List out(n);
  Rcpp::CharacterVector names(n);
  int pos = 0;

  // Append phylogeny fields first under their explicit projection prefix.
  for (int i = 0; i < first.size(); ++i) {
    out[pos] = first[i];
    names[pos] = first_prefix + "." + v34_as_string(first_names[i]);
    ++pos;
  }
  // Append scenario fields second under their explicit projection prefix.
  for (int i = 0; i < second.size(); ++i) {
    out[pos] = second[i];
    names[pos] = second_prefix + "." + v34_as_string(second_names[i]);
    ++pos;
  }
  // Append transition-tree fields last under their explicit projection prefix.
  for (int i = 0; i < third.size(); ++i) {
    out[pos] = third[i];
    names[pos] = third_prefix + "." + v34_as_string(third_names[i]);
    ++pos;
  }

  out.attr("names") = names;
  return out;
}

/**
 * Build tree-length summary-statistics vectors and state partitions.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static Rcpp::NumericVector v34_treelength_stats(
    const V34TreeView& tree,
    const std::string& suffix) {
  int ntips = tree.tip_label.size();
  double total = 0.0;
  double tip_total = 0.0;
  // Partition physical branch length into terminal-tip and internal exposure.
  for (int i = 0; i < tree.edge_length.size(); ++i) {
    total += tree.edge_length[i];
    // Edges ending at public tips contribute to terminal branch length.
    if (tree.edge(i, 1) <= ntips) {
      tip_total += tree.edge_length[i];
    }
  }
  double internal_total = total - tip_total;
  std::vector<std::string> names = {
    "tl.tl." + suffix,
    "tl.internal_bl." + suffix,
    "tl.tip_bl." + suffix,
    "tl.tip_tl_ratio." + suffix,
    "tl.internal_tl_ratio." + suffix
  };
  std::vector<double> values = {
    total,
    internal_total,
    tip_total,
    total == 0.0 ? NA_REAL : tip_total / total,
    total == 0.0 ? NA_REAL : internal_total / total
  };
  return v34_named_numeric(names, values);
}

/**
 * Build tree-length summary-statistics vectors and state partitions.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static Rcpp::List v34_treelength_bystate_stats(
    const V34TreeView& tree,
    const std::string& suffix,
    const std::vector<std::string>& states_set) {
  int ntips = tree.tip_label.size();
  std::vector<double> total(states_set.size(), 0.0);
  std::vector<double> internal_total(states_set.size(), 0.0);
  std::vector<double> tip_total(states_set.size(), 0.0);
  std::map<std::string, int> state_index;
  // Cache each requested state's output position for mapped-edge accumulation.
  for (int i = 0; i < static_cast<int>(states_set.size()); ++i) {
    state_index[states_set[i]] = i;
  }

  std::vector<std::string> tree_states = v34_dimnames_colnames(tree.mapped_edge);
  // Project every tree mapped-edge state column into the shared SS state domain.
  for (int c = 0; c < static_cast<int>(tree_states.size()); ++c) {
    auto idx_it = state_index.find(tree_states[c]);
    // Ignore states outside the requested compatibility domain.
    if (idx_it == state_index.end()) {
      continue;
    }
    int target = idx_it->second;
    // Partition this state's exposure by terminal versus internal edge role.
    for (int r = 0; r < tree.mapped_edge.nrow(); ++r) {
      double value = tree.mapped_edge(r, c);
      total[target] += value;
      // Terminal edges contribute to tip exposure for this state.
      if (tree.edge(r, 1) <= ntips) {
        tip_total[target] += value;
      } else {
        // Nonterminal edges contribute to internal exposure for this state.
        internal_total[target] += value;
      }
    }
  }

  std::vector<double> tip_ratio(states_set.size(), NA_REAL);
  std::vector<double> internal_ratio(states_set.size(), NA_REAL);
  // Derive terminal/internal proportions only for states with positive total exposure.
  for (int i = 0; i < static_cast<int>(states_set.size()); ++i) {
    // Zero-exposure states retain `NA` ratios because no denominator exists.
    if (total[i] != 0.0) {
      tip_ratio[i] = tip_total[i] / total[i];
      internal_ratio[i] = internal_total[i] / total[i];
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("tl_bystate." + suffix) = v34_named_numeric(states_set, total),
    Rcpp::Named("internal_bl_bystate." + suffix) = v34_named_numeric(states_set, internal_total),
    Rcpp::Named("tip_bl_bystate." + suffix) = v34_named_numeric(states_set, tip_total),
    Rcpp::Named("tip_tl_ratio_bystate." + suffix) = v34_named_numeric(states_set, tip_ratio),
    Rcpp::Named("internal_tl_ratio_bystate." + suffix) = v34_named_numeric(states_set, internal_ratio)
  );
}

/**
 * Collect reusable tree and TT summary parts before table assembly.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static V34TreeSummaryParts v34_summary_parts(
    const V34TreeView& tree,
    const std::string& suffix,
    const std::vector<std::string>& states_set,
    bool include_edge_boundaries) {
  V34TreeSummaryParts out;
  out.tipstats = v34_tip_stats(tree, suffix, states_set);
  out.transtats = v34_trans_stats(
    v34_transition_matrix(tree, states_set, include_edge_boundaries),
    suffix
  );
  out.treelength_tree = v34_treelength_stats(tree, suffix);
  out.treelength_bystate = v34_treelength_bystate_stats(tree, suffix, states_set);
  return out;
}

/**
 * Serialize tree-level compatibility accumulators at the same boundary that
 * finalizes a C++ PST tree.
 *
 * The payload is deliberately compact: SS needs grouped tip counts,
 * transition counts, branch lengths, state labels, and terminal-edge roles.
 * It does not retain maps, topology matrices, or the public tree itself.
 */
static Rcpp::List v34_pack_tree_summary_parts(
    const V34TreeView& tree,
    const V34TreeSummaryParts& parts,
    const std::vector<std::string>& states) {
  Rcpp::LogicalVector terminal_edge(tree.edge.nrow());
  int ntips = tree.tip_label.size();
  // Classify each finalized tree edge once while its topology is already hot
  // in materialization. The compact role vector later lets path-size SS
  // separate terminal stops from through support without rereading the tree.
  for (int edge_id = 0; edge_id < tree.edge.nrow(); ++edge_id) {
    terminal_edge[edge_id] = tree.edge(edge_id, 1) <= ntips;
  }
  return Rcpp::List::create(
    Rcpp::Named("tipstats") = parts.tipstats,
    Rcpp::Named("transtats") = parts.transtats,
    Rcpp::Named("treelength_tree") = parts.treelength_tree,
    Rcpp::Named("treelength_bystate") = parts.treelength_bystate,
    Rcpp::Named("state_levels") = states,
    Rcpp::Named("terminal_edge") = terminal_edge
  );
}

/**
 * Restore the compact metric groups captured during PST materialization.
 *
 * Input owns only summary tables, never public maps or topology. Missing groups
 * are a construction error because silently rebuilding from trees would break
 * V34's C++ ownership contract.
 */
static V34TreeSummaryParts v34_unpack_tree_summary_parts(
    const Rcpp::List& packed) {
  // Every compatibility group is required before SS can reproduce the V30
  // field shapes from traversal/materialization-owned data.
  if (!packed.containsElementNamed("tipstats") ||
      !packed.containsElementNamed("transtats") ||
      !packed.containsElementNamed("treelength_tree") ||
      !packed.containsElementNamed("treelength_bystate")) {
    Rcpp::stop("V34 SS tree summary parts are incomplete");
  }
  V34TreeSummaryParts out;
  out.tipstats = Rcpp::as<Rcpp::List>(packed["tipstats"]);
  out.transtats = Rcpp::as<Rcpp::List>(packed["transtats"]);
  out.treelength_tree = Rcpp::as<Rcpp::NumericVector>(
    packed["treelength_tree"]
  );
  out.treelength_bystate = Rcpp::as<Rcpp::List>(
    packed["treelength_bystate"]
  );
  return out;
}

/**
 * Read one required scalar from a named compatibility summary group.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static double v34_list_numeric_scalar(Rcpp::List x, const std::string& name) {
  // A missing compatibility scalar is represented as `NA`, never inferred from another field.
  if (!x.containsElementNamed(name.c_str())) {
    return NA_REAL;
  }
  Rcpp::NumericVector value = x[name];
  return value.size() ? value[0] : NA_REAL;
}

/**
 * Map compatibility summary groups to the saved original-metric vector.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static Rcpp::NumericVector v34_original_metrics(const Rcpp::List& summary) {
  Rcpp::List tipstats = summary["tipstats"];
  Rcpp::List transtats = summary["transtats"];
  Rcpp::List treelengthstats = summary["treelengthstats"];
  Rcpp::List treelength_tree = treelengthstats["tree"];
  Rcpp::List p_tree = treelength_tree["P"];
  Rcpp::List s_tree = treelength_tree["S"];
  Rcpp::List t_tree = treelength_tree["T"];
  Rcpp::List p_trans = transtats["P"];

  double s_tips = v34_list_numeric_scalar(tipstats, "S.ntips.S");
  double t_tips = v34_list_numeric_scalar(tipstats, "T.ntips.T");
  double s_diversity = v34_list_numeric_scalar(s_tree, "tl.tl.S");
  double p_diversity = v34_list_numeric_scalar(p_tree, "tl.tl.P");
  double ntrans = v34_list_numeric_scalar(p_trans, "ntrans.total.P");
  double t_diversity = v34_list_numeric_scalar(t_tree, "tl.tl.T");
  double total_s_internal_nodes = s_tips - 1.0;

  return Rcpp::NumericVector::create(
    Rcpp::Named("s_tips") = s_tips,
    Rcpp::Named("t_tips") = t_tips,
    Rcpp::Named("s_diversity") = s_diversity,
    Rcpp::Named("p_diversity") = p_diversity,
    Rcpp::Named("ntrans") = ntrans,
    Rcpp::Named("t_diversity") = t_diversity,
    Rcpp::Named("total_S_internal_nodes") = total_s_internal_nodes,
    // These tree-wide empirical metrics intentionally preserve the original
    // R/V15 division semantics. When `ntrans` is zero, positive numerators
    // become `Inf`; that is an observed legacy value used for equality checks.
    Rcpp::Named("t_diversityXntrans") = t_diversity / ntrans,
    Rcpp::Named("s_internalXntrans") = total_s_internal_nodes / ntrans,
    Rcpp::Named("s_internalXtdiversity") = t_diversity == 0.0 ? NA_REAL : total_s_internal_nodes / t_diversity
  );
}

/**
 * Distinguish biological value columns from totals and time-axis metadata.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static bool v34_is_value_column(const std::string& name) {
  return name != "total" && name != "time" && name != "summary_time_vec";
}

/**
 * Build a stable matrix-column lookup used by C++ TT and summary-statistics code.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static std::vector<std::string> v34_matrix_value_colnames(SEXP matrix_like) {
  // An absent optional TT matrix contributes no value columns to dependent SS tables.
  if (Rf_isNull(matrix_like)) {
    return std::vector<std::string>();
  }
  Rcpp::NumericMatrix mat(matrix_like);
  Rcpp::List dimnames = mat.attr("dimnames");
  // Unlabelled columns cannot be aligned semantically and are therefore excluded.
  if (dimnames.size() < 2 || Rf_isNull(dimnames[1])) {
    return std::vector<std::string>();
  }
  Rcpp::CharacterVector colnames = dimnames[1];
  std::vector<std::string> out;
  // Keep biological state/path/event columns while dropping totals and time metadata.
  for (int i = 0; i < colnames.size(); ++i) {
    std::string name = v34_as_string(colnames[i]);
    // Only canonical value columns participate in by-state/by-path aggregation.
    if (v34_is_value_column(name)) {
      out.push_back(name);
    }
  }
  return out;
}

/**
 * Read or build a matrix-shaped value while preserving V34 row/time alignment.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static double v34_matrix_final_value(SEXP matrix_like, const std::string& column, double missing = 0.0) {
  // A missing optional TT surface yields the caller's explicit fallback.
  if (Rf_isNull(matrix_like)) {
    return missing;
  }
  Rcpp::NumericMatrix mat(matrix_like);
  // An empty matrix has no final committed value.
  if (mat.nrow() == 0 || mat.ncol() == 0) {
    return missing;
  }
  Rcpp::List dimnames = mat.attr("dimnames");
  // Semantic lookup requires labelled columns; otherwise return the fallback.
  if (dimnames.size() < 2 || Rf_isNull(dimnames[1])) {
    return missing;
  }
  Rcpp::CharacterVector colnames = dimnames[1];
  // Find the requested semantic column and read its last cumulative TT row.
  for (int col = 0; col < colnames.size(); ++col) {
    // Exact label equality identifies the requested state/path/event series.
    if (v34_as_string(colnames[col]) == column) {
      return mat(mat.nrow() - 1, col);
    }
  }
  return missing;
}

/**
 * Read one value from a named numeric vector with an explicit missing fallback.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static double v34_named_numeric_value(SEXP vector_like, const std::string& name, double missing = 0.0) {
  // A missing optional named vector yields the caller's explicit fallback.
  if (Rf_isNull(vector_like)) {
    return missing;
  }
  Rcpp::NumericVector values(vector_like);
  Rcpp::CharacterVector names = values.attr("names");
  // Resolve by semantic name rather than relying on version-specific vector order.
  for (int i = 0; i < values.size(); ++i) {
    // Return only an in-range value whose name matches exactly.
    if (i < names.size() && v34_as_string(names[i]) == name) {
      return values[i];
    }
  }
  return missing;
}

/**
 * Format the tree-wide table from the differently named original metrics.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static Rcpp::DataFrame v34_tree_wide_table(const Rcpp::NumericVector& original_metrics) {
  Rcpp::DataFrame out = Rcpp::DataFrame::create(
    Rcpp::Named("scope") = Rcpp::CharacterVector::create("tree_wide"),
    Rcpp::Named("s_tips") = v34_named_numeric_value(original_metrics, "s_tips", NA_REAL),
    Rcpp::Named("t_tips") = v34_named_numeric_value(original_metrics, "t_tips", NA_REAL),
    Rcpp::Named("s_length") = v34_named_numeric_value(original_metrics, "s_diversity", NA_REAL),
    Rcpp::Named("p_length") = v34_named_numeric_value(original_metrics, "p_diversity", NA_REAL),
    Rcpp::Named("ntrans") = v34_named_numeric_value(original_metrics, "ntrans", NA_REAL),
    Rcpp::Named("t_length") = v34_named_numeric_value(original_metrics, "t_diversity", NA_REAL),
    Rcpp::Named("s_internal") = v34_named_numeric_value(original_metrics, "total_S_internal_nodes", NA_REAL),
    Rcpp::Named("t_length_ntrans") = v34_named_numeric_value(original_metrics, "t_diversityXntrans", NA_REAL),
    Rcpp::Named("s_internal_ntrans") = v34_named_numeric_value(original_metrics, "s_internalXntrans", NA_REAL),
    Rcpp::Named("s_internal_t_length") = v34_named_numeric_value(original_metrics, "s_internalXtdiversity", NA_REAL),
    Rcpp::Named("stringsAsFactors") = false
  );
  return out;
}

/**
 * Read original public metric vectors and add derived transition counts.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static Rcpp::NumericVector v34_original_metrics_with_ntrans(
    Rcpp::NumericVector original_metrics,
    double ntrans) {
  // Preserve captured metrics when no finite traversal transition count is available.
  if (!R_finite(ntrans)) {
    return original_metrics;
  }

  original_metrics["ntrans"] = ntrans;
  double t_diversity = v34_named_numeric_value(original_metrics, "t_diversity", NA_REAL);
  double total_s_internal_nodes = v34_named_numeric_value(original_metrics, "total_S_internal_nodes", NA_REAL);
  original_metrics["t_diversityXntrans"] = t_diversity / ntrans;
  original_metrics["s_internalXntrans"] = total_s_internal_nodes / ntrans;
  original_metrics["s_internalXtdiversity"] =
    t_diversity == 0.0 ? NA_REAL : total_s_internal_nodes / t_diversity;
  return original_metrics;
}

/**
 * Read a finite vector position while treating absent/nonfinite support as missing.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static double v34_map_value_at(Rcpp::NumericVector values, int index, double missing = 0.0) {
  // Out-of-range segment positions contribute the explicit missing-support fallback.
  if (index < 0 || index >= values.size()) {
    return missing;
  }
  double value = values[index];
  return R_finite(value) ? value : missing;
}

/**
 * Assemble part of the V34 summary-statistics payload from public tree and TT objects.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static V34PathSizeSummary v34_path_size_summary(
    SEXP maybe_path_size_maps,
    const Rcpp::LogicalVector& terminal_edge) {
  V34PathSizeSummary out;
  // Missing optional path-size output yields empty dependent support summaries.
  if (Rf_isNull(maybe_path_size_maps)) {
    return out;
  }

  Rcpp::List path_size_maps(maybe_path_size_maps);
  // Both lineage and scenario path maps are required to preserve the paired by-path SS contract.
  if (!path_size_maps.containsElementNamed("lineage_size_path_maps") ||
      !path_size_maps.containsElementNamed("scenario_size_path_maps")) {
    return out;
  }

  Rcpp::List lineage_maps = path_size_maps["lineage_size_path_maps"];
  Rcpp::List scenario_maps = path_size_maps["scenario_size_path_maps"];
  int nedge = std::min(
    static_cast<int>(terminal_edge.size()),
    std::min(static_cast<int>(lineage_maps.size()), static_cast<int>(scenario_maps.size()))
  );

  // Classify each path-size edge by terminal role and aggregate its segment support by path label.
  for (int edge_i = 0; edge_i < nedge; ++edge_i) {
    Rcpp::NumericVector lineage_map = lineage_maps[edge_i];
    Rcpp::NumericVector scenario_map = scenario_maps[edge_i];
    std::vector<std::string> lineage_labels = v34_character_names(lineage_map);
    std::vector<std::string> scenario_labels = v34_character_names(scenario_map);

    // Terminal transition-tree edges end in STOPPED support; their final map
    // entry is a stop count while every preceding entry is through support.
    if (terminal_edge[edge_i] && !lineage_labels.empty()) {
      int lineage_last = static_cast<int>(lineage_labels.size()) - 1;
      int scenario_last = static_cast<int>(scenario_labels.size()) - 1;
      // Sum lineage support before the final STOPPED segment as through support.
      for (int segment_i = 0; segment_i < lineage_last; ++segment_i) {
        const std::string& label = lineage_labels[segment_i];
        // Every nonterminal segment is a distinct through-support contribution;
        // aggregate repeated path labels instead of retaining only their maximum.
        out.lineage_through[label] +=
          v34_map_value_at(lineage_map, segment_i);
      }
      // Sum scenario support before the final STOPPED segment as through support.
      for (int segment_i = 0; segment_i < scenario_last; ++segment_i) {
        const std::string& label = scenario_labels[segment_i];
        // Scenario-through support follows the same role-aware segment sum.
        out.scenario_through[label] +=
          v34_map_value_at(scenario_map, segment_i);
      }
      out.lineage_stopped[lineage_labels[lineage_last]] += v34_map_value_at(lineage_map, lineage_last);
      // A nonempty scenario map contributes its final segment to stopped support.
      if (scenario_last >= 0) {
        out.scenario_stopped[scenario_labels[scenario_last]] += v34_map_value_at(scenario_map, scenario_last);
      }
    } else {
      // Every segment on a nonterminal edge represents through support rather than a stop.
      for (int segment_i = 0; segment_i < static_cast<int>(lineage_labels.size()); ++segment_i) {
        const std::string& label = lineage_labels[segment_i];
        // All segments on a nonterminal transition-tree edge are through
        // carriers and therefore contribute additively to the path total.
        out.lineage_through[label] +=
          v34_map_value_at(lineage_map, segment_i);
      }
      // Aggregate all nonterminal scenario support by canonical path label.
      for (int segment_i = 0; segment_i < static_cast<int>(scenario_labels.size()); ++segment_i) {
        const std::string& label = scenario_labels[segment_i];
        // Sum repeated scenario path labels across nonterminal edges.
        out.scenario_through[label] +=
          v34_map_value_at(scenario_map, segment_i);
      }
    }
  }

  return out;
}

/**
 * Translate integer ids to public labels through an R lookup table.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static double v34_map_lookup(const std::map<std::string, double>& values, const std::string& key) {
  auto it = values.find(key);
  return it == values.end() ? 0.0 : it->second;
}

/**
 * Fetch an optional named list component without manufacturing a placeholder.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static SEXP v34_list_element_or_null(Rcpp::List x, const std::string& name) {
  // Preserve `NULL` when the optional component is absent; callers decide the fallback semantics.
  if (!x.containsElementNamed(name.c_str())) {
    return R_NilValue;
  }
  return x[name];
}

/**
 * Fetch an optional TT matrix field while preserving NULL for absent families.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static SEXP v34_tt_matrix_or_null(
    Rcpp::Nullable<Rcpp::List> maybe_tt,
    const std::string& group,
    const std::string& component) {
  // SS can format tree-only fields without TT; absent TT yields no dependent matrix.
  if (maybe_tt.isNull()) {
    return R_NilValue;
  }
  Rcpp::List tt(maybe_tt);
  // A missing TT family provides no dependent SS values.
  if (!tt.containsElementNamed(group.c_str())) {
    return R_NilValue;
  }
  Rcpp::List tt_group = Rcpp::as<Rcpp::List>(tt[group]);
  // A family lacking the requested matrix component contributes no values.
  if (!tt_group.containsElementNamed(component.c_str())) {
    return R_NilValue;
  }
  return tt_group[component];
}

/**
 * Read one state-specific branch-length total captured during tree
 * materialization.
 *
 * The named vector is the exact output of `v34_treelength_bystate_stats`; a
 * missing state has zero exposure, while a missing vector indicates a broken
 * C++ materialization contract and is rejected.
 */
static double v34_summary_part_state_total(
    const V34TreeSummaryParts& parts,
    const std::string& suffix,
    const std::string& state) {
  std::string field = "tl_bystate." + suffix;
  // The materializer must provide the complete state-length vector for each of
  // P, S, and T; rebuilding it from a public map is intentionally disallowed.
  if (!parts.treelength_bystate.containsElementNamed(field.c_str())) {
    Rcpp::stop("V34 SS tree summary parts omit state branch lengths");
  }
  return v34_named_numeric_value(
    parts.treelength_bystate[field],
    state,
    0.0
  );
}

/**
 * Parse, count, or summarize path-history labels for transition/path statistics.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static std::string v34_path_terminal_state(const std::string& path) {
  std::vector<std::string> pieces = v34_split_path(path);
  // Scan backward to skip root-prefix empties and return the final biological state step.
  for (int i = static_cast<int>(pieces.size()) - 1; i >= 0; --i) {
    // The first nonempty component from the end is the path's terminal state.
    if (!pieces[i].empty()) {
      return pieces[i];
    }
  }
  return std::string();
}

/**
 * Extract or transform event labels for V34 summary-statistics tables.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static std::string v34_event_from_state(const std::string& event) {
  std::size_t pos = event.find("->");
  return pos == std::string::npos ? std::string() : event.substr(0, pos);
}

/**
 * Extract or transform event labels for V34 summary-statistics tables.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static std::string v34_event_to_state(const std::string& event) {
  std::size_t pos = event.find("->");
  return pos == std::string::npos ? std::string() : event.substr(pos + 2);
}

/**
 * Compute a ratio while representing an undefined zero denominator as `NA`.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static double v34_safe_ratio(double numerator, double denominator) {
  return denominator == 0.0 ? NA_REAL : numerator / denominator;
}

/**
 * Build primary lineage summary tables grouped by state.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static Rcpp::List v34_primary_by_state_tables(
    const V34TreeSummaryParts& p,
    const V34TreeSummaryParts& s,
    const V34TreeSummaryParts& t,
    const std::vector<std::string>& states,
    Rcpp::Nullable<Rcpp::List> tt) {
  SEXP tlen_state = v34_tt_matrix_or_null(tt, "tlen", "state");
  SEXP puniq_arr = v34_tt_matrix_or_null(tt, "puniq", "arr");
  SEXP puniq_leave = v34_tt_matrix_or_null(tt, "puniq", "leave");
  SEXP ldif_arr = v34_tt_matrix_or_null(tt, "ldif", "arr");
  SEXP ldif_leave = v34_tt_matrix_or_null(tt, "ldif", "leave");
  SEXP ctt_arr = v34_tt_matrix_or_null(tt, "ctt", "arr");
  SEXP ctt_leave = v34_tt_matrix_or_null(tt, "ctt", "leave");

  Rcpp::CharacterVector state(states.size());
  Rcpp::NumericVector phylo_time(states.size());
  Rcpp::NumericVector scenario_time(states.size());
  Rcpp::NumericVector transition_time(states.size());
  Rcpp::NumericVector integrated_lineage_time(states.size());
  Rcpp::NumericVector scenario_tip_count(states.size());
  Rcpp::NumericVector transition_tip_count(states.size());
  Rcpp::NumericVector path_unique_arrivals(states.size());
  Rcpp::NumericVector path_unique_leaves(states.size());
  Rcpp::NumericVector lineage_differentiating_arrivals(states.size());
  Rcpp::NumericVector lineage_differentiating_leaves(states.size());
  Rcpp::NumericVector all_transition_arrivals(states.size());
  Rcpp::NumericVector all_transition_leaves(states.size());
  Rcpp::NumericVector punique_arr_per_all_trans_arr(states.size());
  Rcpp::NumericVector punique_leave_per_all_trans_leave(states.size());
  Rcpp::NumericVector lindiff_arr_per_all_trans_arr(states.size());
  Rcpp::NumericVector lindiff_leave_per_all_trans_leave(states.size());
  Rcpp::NumericVector all_trans_arr_per_phylo_length(states.size());
  Rcpp::NumericVector all_trans_leave_per_phylo_length(states.size());
  Rcpp::NumericVector all_trans_arr_per_scenario_length(states.size());
  Rcpp::NumericVector all_trans_leave_per_scenario_length(states.size());
  Rcpp::NumericVector all_trans_arr_per_transition_length(states.size());
  Rcpp::NumericVector all_trans_leave_per_transition_length(states.size());
  Rcpp::NumericVector punique_arr_per_phylo_length(states.size());
  Rcpp::NumericVector punique_leave_per_phylo_length(states.size());
  Rcpp::NumericVector punique_arr_per_scenario_length(states.size());
  Rcpp::NumericVector punique_leave_per_scenario_length(states.size());
  Rcpp::NumericVector punique_arr_per_transition_length(states.size());
  Rcpp::NumericVector punique_leave_per_transition_length(states.size());
  Rcpp::NumericVector lindiff_arr_per_phylo_length(states.size());
  Rcpp::NumericVector lindiff_leave_per_phylo_length(states.size());
  Rcpp::NumericVector lindiff_arr_per_scenario_length(states.size());
  Rcpp::NumericVector lindiff_leave_per_scenario_length(states.size());
  Rcpp::NumericVector lindiff_arr_per_transition_length(states.size());
  Rcpp::NumericVector lindiff_leave_per_transition_length(states.size());

  SEXP scenario_tips = v34_list_element_or_null(s.tipstats, "ntips.bystate.S");
  SEXP transition_tips = v34_list_element_or_null(t.tipstats, "ntips.bystate.T");

  // Build one raw and normalized SS row per declared state from canonical
  // tree-length parts and final online TT counters.
  for (int i = 0; i < static_cast<int>(states.size()); ++i) {
    const std::string& label = states[i];
    state[i] = label;
    phylo_time[i] = v34_summary_part_state_total(p, "P", label);
    scenario_time[i] = v34_summary_part_state_total(s, "S", label);
    transition_time[i] = v34_summary_part_state_total(t, "T", label);
    integrated_lineage_time[i] = v34_matrix_final_value(tlen_state, label);
    scenario_tip_count[i] = v34_named_numeric_value(scenario_tips, label);
    transition_tip_count[i] = v34_named_numeric_value(transition_tips, label);
    path_unique_arrivals[i] = v34_matrix_final_value(puniq_arr, label);
    path_unique_leaves[i] = v34_matrix_final_value(puniq_leave, label);
    lineage_differentiating_arrivals[i] = v34_matrix_final_value(ldif_arr, label);
    lineage_differentiating_leaves[i] = v34_matrix_final_value(ldif_leave, label);
    all_transition_arrivals[i] = v34_matrix_final_value(ctt_arr, label);
    all_transition_leaves[i] = v34_matrix_final_value(ctt_leave, label);
    punique_arr_per_all_trans_arr[i] = v34_safe_ratio(path_unique_arrivals[i], all_transition_arrivals[i]);
    punique_leave_per_all_trans_leave[i] = v34_safe_ratio(path_unique_leaves[i], all_transition_leaves[i]);
    lindiff_arr_per_all_trans_arr[i] = v34_safe_ratio(lineage_differentiating_arrivals[i], all_transition_arrivals[i]);
    lindiff_leave_per_all_trans_leave[i] = v34_safe_ratio(lineage_differentiating_leaves[i], all_transition_leaves[i]);
    all_trans_arr_per_phylo_length[i] = v34_safe_ratio(all_transition_arrivals[i], phylo_time[i]);
    all_trans_leave_per_phylo_length[i] = v34_safe_ratio(all_transition_leaves[i], phylo_time[i]);
    all_trans_arr_per_scenario_length[i] = v34_safe_ratio(all_transition_arrivals[i], scenario_time[i]);
    all_trans_leave_per_scenario_length[i] = v34_safe_ratio(all_transition_leaves[i], scenario_time[i]);
    all_trans_arr_per_transition_length[i] = v34_safe_ratio(all_transition_arrivals[i], transition_time[i]);
    all_trans_leave_per_transition_length[i] = v34_safe_ratio(all_transition_leaves[i], transition_time[i]);
    punique_arr_per_phylo_length[i] = v34_safe_ratio(path_unique_arrivals[i], phylo_time[i]);
    punique_leave_per_phylo_length[i] = v34_safe_ratio(path_unique_leaves[i], phylo_time[i]);
    punique_arr_per_scenario_length[i] = v34_safe_ratio(path_unique_arrivals[i], scenario_time[i]);
    punique_leave_per_scenario_length[i] = v34_safe_ratio(path_unique_leaves[i], scenario_time[i]);
    punique_arr_per_transition_length[i] = v34_safe_ratio(path_unique_arrivals[i], transition_time[i]);
    punique_leave_per_transition_length[i] = v34_safe_ratio(path_unique_leaves[i], transition_time[i]);
    lindiff_arr_per_phylo_length[i] = v34_safe_ratio(lineage_differentiating_arrivals[i], phylo_time[i]);
    lindiff_leave_per_phylo_length[i] = v34_safe_ratio(lineage_differentiating_leaves[i], phylo_time[i]);
    lindiff_arr_per_scenario_length[i] = v34_safe_ratio(lineage_differentiating_arrivals[i], scenario_time[i]);
    lindiff_leave_per_scenario_length[i] = v34_safe_ratio(lineage_differentiating_leaves[i], scenario_time[i]);
    lindiff_arr_per_transition_length[i] = v34_safe_ratio(lineage_differentiating_arrivals[i], transition_time[i]);
    lindiff_leave_per_transition_length[i] = v34_safe_ratio(lineage_differentiating_leaves[i], transition_time[i]);
  }

  Rcpp::DataFrame by_state = Rcpp::DataFrame::create(
    Rcpp::Named("state") = state,
    Rcpp::Named("p_len") = phylo_time,
    Rcpp::Named("s_len") = scenario_time,
    Rcpp::Named("t_len") = transition_time,
    Rcpp::Named("tlen") = integrated_lineage_time,
    Rcpp::Named("s_tips") = scenario_tip_count,
    Rcpp::Named("t_tips") = transition_tip_count,
    Rcpp::Named("all_trans_arr") = all_transition_arrivals,
    Rcpp::Named("all_trans_leave") = all_transition_leaves,
    Rcpp::Named("punique_arr") = path_unique_arrivals,
    Rcpp::Named("punique_leave") = path_unique_leaves,
    Rcpp::Named("lindiff_arr") = lineage_differentiating_arrivals,
    Rcpp::Named("lindiff_leave") = lineage_differentiating_leaves,
    Rcpp::Named("stringsAsFactors") = false
  );

  Rcpp::DataFrame by_state_norm = Rcpp::DataFrame::create(
    Rcpp::Named("state") = state,
    Rcpp::Named("punique_arr|all_trans_arr") = punique_arr_per_all_trans_arr,
    Rcpp::Named("punique_leave|all_trans_leave") = punique_leave_per_all_trans_leave,
    Rcpp::Named("lindiff_arr|all_trans_arr") = lindiff_arr_per_all_trans_arr,
    Rcpp::Named("lindiff_leave|all_trans_leave") = lindiff_leave_per_all_trans_leave,
    Rcpp::Named("all_trans_arr|p_len") = all_trans_arr_per_phylo_length,
    Rcpp::Named("all_trans_leave|p_len") = all_trans_leave_per_phylo_length,
    Rcpp::Named("all_trans_arr|s_len") = all_trans_arr_per_scenario_length,
    Rcpp::Named("all_trans_leave|s_len") = all_trans_leave_per_scenario_length,
    Rcpp::Named("all_trans_arr|t_len") = all_trans_arr_per_transition_length,
    Rcpp::Named("all_trans_leave|t_len") = all_trans_leave_per_transition_length,
    Rcpp::Named("punique_arr|p_len") = punique_arr_per_phylo_length,
    Rcpp::Named("punique_leave|p_len") = punique_leave_per_phylo_length,
    Rcpp::Named("punique_arr|s_len") = punique_arr_per_scenario_length,
    Rcpp::Named("punique_leave|s_len") = punique_leave_per_scenario_length,
    Rcpp::Named("punique_arr|t_len") = punique_arr_per_transition_length,
    Rcpp::Named("punique_leave|t_len") = punique_leave_per_transition_length,
    Rcpp::Named("lindiff_arr|p_len") = lindiff_arr_per_phylo_length,
    Rcpp::Named("lindiff_leave|p_len") = lindiff_leave_per_phylo_length,
    Rcpp::Named("lindiff_arr|s_len") = lindiff_arr_per_scenario_length,
    Rcpp::Named("lindiff_leave|s_len") = lindiff_leave_per_scenario_length,
    Rcpp::Named("lindiff_arr|t_len") = lindiff_arr_per_transition_length,
    Rcpp::Named("lindiff_leave|t_len") = lindiff_leave_per_transition_length,
    Rcpp::Named("stringsAsFactors") = false
  );
  by_state_norm.attr("names") = Rcpp::CharacterVector::create(
    "state",
    "punique_arr|all_trans_arr",
    "punique_leave|all_trans_leave",
    "lindiff_arr|all_trans_arr",
    "lindiff_leave|all_trans_leave",
    "all_trans_arr|p_len",
    "all_trans_leave|p_len",
    "all_trans_arr|s_len",
    "all_trans_leave|s_len",
    "all_trans_arr|t_len",
    "all_trans_leave|t_len",
    "punique_arr|p_len",
    "punique_leave|p_len",
    "punique_arr|s_len",
    "punique_leave|s_len",
    "punique_arr|t_len",
    "punique_leave|t_len",
    "lindiff_arr|p_len",
    "lindiff_leave|p_len",
    "lindiff_arr|s_len",
    "lindiff_leave|s_len",
    "lindiff_arr|t_len",
    "lindiff_leave|t_len"
  );

  return Rcpp::List::create(
    Rcpp::Named("by_state") = by_state,
    Rcpp::Named("by_state_norm") = by_state_norm
  );
}

/**
 * Build primary lineage summary tables grouped by path.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static Rcpp::DataFrame v34_primary_by_path(
    Rcpp::Nullable<Rcpp::List> tt,
    SEXP path_size_maps,
    const Rcpp::LogicalVector& terminal_edge) {
  // By-path summary contract:
  // `lin_through`/`scn_through` are through/internal path counts, while
  // `lin_tip`/`scn_tip` are terminal path counts. The fallback from tip to
  // through is allowed only when a path has no separate through carrier because
  // every reached lineage or scenario terminates there.
  SEXP ltt_trans = v34_tt_matrix_or_null(tt, "ltt", "trans");
  SEXP stt_trans = v34_tt_matrix_or_null(tt, "stt", "trans");
  SEXP cltt_trans = v34_tt_matrix_or_null(tt, "cltt", "trans");
  SEXP cstt_trans = v34_tt_matrix_or_null(tt, "cstt", "trans");
  SEXP seq_trans = v34_tt_matrix_or_null(tt, "tt", "trans");
  SEXP tlen_trans = v34_tt_matrix_or_null(tt, "tlen", "trans");
  SEXP puniq_trans = v34_tt_matrix_or_null(tt, "puniq", "trans");
  SEXP ldif_trans = v34_tt_matrix_or_null(tt, "ldif", "trans");
  V34PathSizeSummary path_sizes = v34_path_size_summary(
    path_size_maps,
    terminal_edge
  );

  std::set<std::string> label_set;
  // Seed the canonical path domain from online lineage-through TT columns.
  for (const std::string& label : v34_matrix_value_colnames(ltt_trans)) label_set.insert(label);
  // Add paths represented by online scenario-through TT columns.
  for (const std::string& label : v34_matrix_value_colnames(stt_trans)) label_set.insert(label);
  // Cumulative TT matrices are the canonical source for through counts.
  for (const std::string& label : v34_matrix_value_colnames(cltt_trans)) label_set.insert(label);
  for (const std::string& label : v34_matrix_value_colnames(cstt_trans)) label_set.insert(label);
  // Add paths represented by active transition-sequence TT columns.
  for (const std::string& label : v34_matrix_value_colnames(seq_trans)) label_set.insert(label);
  // Add paths carrying integrated lineage exposure.
  for (const std::string& label : v34_matrix_value_colnames(tlen_trans)) label_set.insert(label);
  // Add paths with a scored first public appearance.
  for (const std::string& label : v34_matrix_value_colnames(puniq_trans)) label_set.insert(label);
  // Add paths with lineage-differentiating transition counts.
  for (const std::string& label : v34_matrix_value_colnames(ldif_trans)) label_set.insert(label);
  // Add paths present only in traversal-owned lineage size support.
  for (const auto& item : path_sizes.lineage_through) label_set.insert(item.first);
  // Add paths present only in traversal-owned scenario size support.
  for (const auto& item : path_sizes.scenario_through) label_set.insert(item.first);
  // Add terminal-only lineage paths that may have no active TT interval.
  for (const auto& item : path_sizes.lineage_stopped) label_set.insert(item.first);
  // Add terminal-only scenario paths that may have no active TT interval.
  for (const auto& item : path_sizes.scenario_stopped) label_set.insert(item.first);
  std::vector<std::string> labels(label_set.begin(), label_set.end());

  Rcpp::CharacterVector path(labels.size());
  Rcpp::CharacterVector state(labels.size());
  Rcpp::NumericVector lineage_through(labels.size());
  Rcpp::NumericVector scenario_through(labels.size());
  Rcpp::NumericVector lineage_tip(labels.size());
  Rcpp::NumericVector scenario_tip(labels.size());
  Rcpp::NumericVector integrated_path_time(labels.size());
  Rcpp::NumericVector path_unique_arrivals(labels.size());
  Rcpp::NumericVector lineage_differentiating_arrivals(labels.size());

  // Format one by-path row per semantic path, reading final TT values and role-aware size support.
  for (int i = 0; i < static_cast<int>(labels.size()); ++i) {
    const std::string& label = labels[i];
    path[i] = label;
    state[i] = v34_path_terminal_state(label);
    lineage_tip[i] = v34_map_lookup(path_sizes.lineage_stopped, label);
    scenario_tip[i] = v34_map_lookup(path_sizes.scenario_stopped, label);
    lineage_through[i] = v34_matrix_final_value(cltt_trans, label);
    scenario_through[i] = v34_matrix_final_value(cstt_trans, label);
    integrated_path_time[i] = v34_matrix_final_value(tlen_trans, label);
    path_unique_arrivals[i] = v34_matrix_final_value(puniq_trans, label);
    lineage_differentiating_arrivals[i] = v34_matrix_final_value(ldif_trans, label);
  }

  Rcpp::DataFrame out = Rcpp::DataFrame::create(
    Rcpp::Named("path") = path,
    Rcpp::Named("state") = state,
    Rcpp::Named("lin_through") = lineage_through,
    Rcpp::Named("scn_through") = scenario_through,
    Rcpp::Named("lin_tip") = lineage_tip,
    Rcpp::Named("scn_tip") = scenario_tip,
    Rcpp::Named("tlen") = integrated_path_time,
    Rcpp::Named("punique") = path_unique_arrivals,
    Rcpp::Named("lindiff") = lineage_differentiating_arrivals,
    Rcpp::Named("stringsAsFactors") = false
  );
  return out;
}

/**
 * Extract or transform event labels for V34 summary-statistics tables.
 *
 * Inputs are read only; outputs preserve canonical labels, dimensions, and R-facing ids.
 */
static Rcpp::DataFrame v34_primary_by_event(Rcpp::Nullable<Rcpp::List> tt) {
  SEXP ctt_event = v34_tt_matrix_or_null(tt, "ctt", "event");
  SEXP ldif_event = v34_tt_matrix_or_null(tt, "ldif", "event");
  SEXP puniq_event = v34_tt_matrix_or_null(tt, "puniq", "event");
  std::set<std::string> label_set;
  // Seed the event domain from all scored biological transition types.
  for (const std::string& label : v34_matrix_value_colnames(ctt_event)) label_set.insert(label);
  // Add event types appearing only in lineage-differentiation counts.
  for (const std::string& label : v34_matrix_value_colnames(ldif_event)) label_set.insert(label);
  // Add event types appearing only in first-path-appearance counts.
  for (const std::string& label : v34_matrix_value_colnames(puniq_event)) label_set.insert(label);
  std::vector<std::string> labels(label_set.begin(), label_set.end());

  Rcpp::CharacterVector event(labels.size());
  Rcpp::CharacterVector from(labels.size());
  Rcpp::CharacterVector to(labels.size());
  Rcpp::NumericVector change_count(labels.size());
  Rcpp::NumericVector lineage_differentiating_count(labels.size());
  Rcpp::NumericVector path_unique_count(labels.size());

  // Format one by-event row from the final cumulative TT counters for each semantic event label.
  for (int i = 0; i < static_cast<int>(labels.size()); ++i) {
    const std::string& label = labels[i];
    event[i] = label;
    from[i] = v34_event_from_state(label);
    to[i] = v34_event_to_state(label);
    change_count[i] = v34_matrix_final_value(ctt_event, label);
    lineage_differentiating_count[i] = v34_matrix_final_value(ldif_event, label);
    path_unique_count[i] = v34_matrix_final_value(puniq_event, label);
  }

  Rcpp::DataFrame out = Rcpp::DataFrame::create(
    Rcpp::Named("event") = event,
    Rcpp::Named("from") = from,
    Rcpp::Named("to") = to,
    Rcpp::Named("all_trans") = change_count,
    Rcpp::Named("lindiff") = lineage_differentiating_count,
    Rcpp::Named("punique") = path_unique_count,
    Rcpp::Named("stringsAsFactors") = false
  );
  return out;
}

/**
 * Assemble all primary lineage summary tables into the public SS payload.
 *
 * Inputs are read only and failures stop before a partial SS table crosses the R boundary.
 */
static Rcpp::List v34_primary_summary_tables(
    const Rcpp::NumericVector& tree_wide,
    const V34TreeSummaryParts& p,
    const V34TreeSummaryParts& s,
    const V34TreeSummaryParts& t,
    const std::vector<std::string>& states,
    Rcpp::Nullable<Rcpp::List> tt,
    SEXP path_size_maps,
    const Rcpp::LogicalVector& terminal_edge) {
  Rcpp::List by_state_tables = v34_primary_by_state_tables(
    p,
    s,
    t,
    states,
    tt
  );
  return Rcpp::List::create(
    Rcpp::Named("tree_wide") = v34_tree_wide_table(tree_wide),
    Rcpp::Named("by_state") = by_state_tables["by_state"],
    Rcpp::Named("by_state_norm") = by_state_tables["by_state_norm"],
    Rcpp::Named("by_path") = v34_primary_by_path(
      tt,
      path_size_maps,
      terminal_edge
    ),
    Rcpp::Named("by_event") = v34_primary_by_event(tt),
    Rcpp::Named("metadata") = Rcpp::List::create(
      Rcpp::Named("by_state_table") = "raw counts and lengths by terminal state",
      Rcpp::Named("by_state_norm_table") = "directional arrival and leaving ratios by terminal state",
      Rcpp::Named("ratio_denominator") = "column names encode denominator after |",
      Rcpp::Named("path_specific_definition") = "canonical_tt_path_lookup",
      // Saved all-SS references require this metadata scalar. It documents
      // combined event semantics but does not create an `event_family` object.
      Rcpp::Named("event_family_status") = "combined_only_clado_anagenetic_reserved"
    )
  );
}

/**
 * Assemble the stable V34 SS object from compact tree-materialization metrics
 * and traversal-owned TT/path support.
 *
 * This is the single SS construction path for both the public online boundary
 * and the retained compatibility wrapper. It never receives a public tree and
 * therefore cannot rediscover biological values from completed maps.
 */
static Rcpp::List v34_summary_from_parts(
    const V34TreeSummaryParts& p,
    const V34TreeSummaryParts& s,
    const V34TreeSummaryParts& t,
    const std::vector<std::string>& states_set,
    const Rcpp::LogicalVector& trans_terminal_edge,
    SEXP path_size_maps,
    Rcpp::Nullable<Rcpp::List> tt) {
  Rcpp::List tipstats = v34_concat_named_lists(
    p.tipstats,
    s.tipstats,
    t.tipstats,
    "P",
    "S",
    "T"
  );
  Rcpp::List transtats = Rcpp::List::create(
    Rcpp::Named("P") = p.transtats,
    Rcpp::Named("S") = s.transtats,
    Rcpp::Named("T") = t.transtats
  );
  Rcpp::List treelength_tree = Rcpp::List::create(
    Rcpp::Named("P") = p.treelength_tree,
    Rcpp::Named("S") = s.treelength_tree,
    Rcpp::Named("T") = t.treelength_tree
  );
  Rcpp::List treelength_bystate = Rcpp::List::create(
    Rcpp::Named("P") = p.treelength_bystate,
    Rcpp::Named("S") = s.treelength_bystate,
    Rcpp::Named("T") = t.treelength_bystate
  );
  Rcpp::List summary = Rcpp::List::create(
    Rcpp::Named("tipstats") = tipstats,
    Rcpp::Named("transtats") = transtats,
    Rcpp::Named("treelengthstats") = Rcpp::List::create(
      Rcpp::Named("tree") = treelength_tree,
      Rcpp::Named("bystate") = treelength_bystate
    )
  );

  // CTT is authoritative for the public transition count because it includes
  // every biological edge event, including cladogenetic state departures that
  // are intentionally absent from the phylogeny map-only compatibility table.
  Rcpp::NumericVector original_metrics = v34_original_metrics(summary);
  original_metrics = v34_original_metrics_with_ntrans(
    original_metrics,
    v34_matrix_final_value(
      v34_tt_matrix_or_null(tt, "ctt", "event"),
      "total",
      NA_REAL
    )
  );
  Rcpp::List primary_tables = v34_primary_summary_tables(
    original_metrics,
    p,
    s,
    t,
    states_set,
    tt,
    path_size_maps,
    trans_terminal_edge
  );

  // Preserve the V30-shaped field order while recording that all biological
  // inputs were captured before this pure table-formatting step.
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("original_metrics") = original_metrics,
    Rcpp::Named("tree_wide") = primary_tables["tree_wide"],
    Rcpp::Named("by_state") = primary_tables["by_state"],
    Rcpp::Named("by_state_norm") = primary_tables["by_state_norm"],
    Rcpp::Named("by_path") = primary_tables["by_path"],
    Rcpp::Named("by_event") = primary_tables["by_event"],
    Rcpp::Named("summary") = summary,
    Rcpp::Named("table_metadata") = primary_tables["metadata"],
    Rcpp::Named("metadata") = Rcpp::List::create(
      Rcpp::Named("n_states") = static_cast<int>(states_set.size()),
      Rcpp::Named("summary_engine") = "cpp",
      Rcpp::Named("summary_contract") = "anolis_empirical_pst_summary_stats_flatten_false",
      Rcpp::Named("primary_summary_surface") = "ss",
      Rcpp::Named("construction_source") =
        "tt_tree_parts_and_traversal_path_size_maps"
    )
  );
  out.attr("class") = Rcpp::CharacterVector::create(
    "pst_v34_summary_stats",
    "list"
  );
  return out;
}

} // namespace

/**
 * Capture exact compatibility metrics while C++ finalizes one public PST tree.
 *
 * This function is called from the tree-output stage, not from SS finalization.
 * The returned payload owns no public tree or dense scenario matrix.
 */
Rcpp::List pst_v34_capture_tree_summary_parts_cpp(
    Rcpp::List tree,
    const std::string& suffix,
    const Rcpp::CharacterVector& state_levels,
    bool include_edge_boundaries) {
  V34TreeView view = v34_tree_view(tree);
  std::vector<std::string> states = Rcpp::as<std::vector<std::string> >(
    state_levels
  );
  V34TreeSummaryParts parts = v34_summary_parts(
    view,
    suffix,
    states,
    include_edge_boundaries
  );
  return v34_pack_tree_summary_parts(view, parts, states);
}

/**
 * Assemble V34 SS exclusively from precomputed C++ tree parts and online TT.
 *
 * Inputs are read only. Missing P/S/T parts or terminal-edge roles stop before
 * an incomplete SS object can cross the R boundary.
 */
Rcpp::List pst_v34_build_summary_stats_from_parts_cpp(
    Rcpp::List tree_parts,
    SEXP path_size_maps,
    Rcpp::Nullable<Rcpp::List> tt) {
  // All three PST projections are needed because the compatibility summary
  // retains separate phylogeny, scenario, and transition metrics.
  if (!tree_parts.containsElementNamed("P") ||
      !tree_parts.containsElementNamed("S") ||
      !tree_parts.containsElementNamed("T")) {
    Rcpp::stop("V34 SS requires P, S, and T tree summary parts");
  }
  Rcpp::List p_packed = tree_parts["P"];
  Rcpp::List s_packed = tree_parts["S"];
  Rcpp::List t_packed = tree_parts["T"];
  // Transition terminal-edge roles are required to distinguish stopped path
  // support from through support in `ss$by_path`.
  if (!t_packed.containsElementNamed("terminal_edge") ||
      !p_packed.containsElementNamed("state_levels")) {
    Rcpp::stop("V34 SS tree summary parts omit identity metadata");
  }
  V34TreeSummaryParts p = v34_unpack_tree_summary_parts(p_packed);
  V34TreeSummaryParts s = v34_unpack_tree_summary_parts(s_packed);
  V34TreeSummaryParts t = v34_unpack_tree_summary_parts(t_packed);
  std::vector<std::string> states = Rcpp::as<std::vector<std::string> >(
    p_packed["state_levels"]
  );
  Rcpp::LogicalVector terminal_edge = t_packed["terminal_edge"];
  return v34_summary_from_parts(
    p,
    s,
    t,
    states,
    terminal_edge,
    path_size_maps,
    tt
  );
}
