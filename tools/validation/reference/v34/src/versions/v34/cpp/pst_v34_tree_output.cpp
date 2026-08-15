#include "pst_v34_tree_output.hpp"

#include "pst_v34_batch.hpp"
#include "pst_v34_ss_output.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct V34TransitionMap {
  std::vector<std::string> labels;
  std::vector<double> durations;
};

struct V34ScenarioEdge {
  int parent;
  int child;
  V34TransitionMap map;
  V34TransitionMap path_map;
  std::set<int> raw_edge_ids;
};

struct V34ScenarioTree {
  std::vector<V34ScenarioEdge> edges;
  std::vector<std::string> tip_labels;
  int nnode;
  bool has_path_maps = false;
};

/**
 * Join represented tip labels in stable public tip order after scenario-tip merging.
 *
 * Inputs are read only and failures stop before a partially formatted public tree is returned.
 */
static std::string pst_v34_join_tip_labels(
    const std::vector<int>& tip_indices,
    const Rcpp::CharacterVector& tip_labels) {
  std::ostringstream out;
  int emitted = 0;
  // Visit represented biological tips in public order and append each valid
  // label exactly once to the merged scenario-tip label.
  for (int index : tip_indices) {
    // Ignore missing or out-of-domain tip ids because they do not identify a
    // biological descendant that may appear in the public merged label.
    if (index == NA_INTEGER || index < 1 || index > tip_labels.size()) {
      continue;
    }
    // Very large merged labels are line-broken without changing their comma-
    // separated descendant identity.
    if (emitted > 0 && emitted % 9999 == 0) {
      out << "\n";
    } else if (emitted > 0) {
      // Ordinary later descendants are separated from the preceding label.
      out << ",";
    }
    out << Rcpp::as<std::string>(tip_labels[index - 1]);
    ++emitted;
  }
  return out.str();
}

/**
 * Create an internal transition-map segment vector from labels and durations.
 *
 * Inputs are read only and failures stop before a partially formatted public tree is returned.
 */
static V34TransitionMap pst_v34_make_transition_map(
    const std::vector<int>& state_ids,
    const std::vector<int>& lineage_counts,
    const std::vector<double>& durations,
    const std::map<int, std::string>& state_labels) {
  V34TransitionMap out;
  // Convert each traversal segment into its public `STATE_LINEAGES` identity;
  // the three input vectors are segment-parallel by construction.
  for (int i = 0; i < static_cast<int>(state_ids.size()); ++i) {
    auto label_it = state_labels.find(state_ids[i]);
    std::string state = label_it == state_labels.end() ? std::string("NA") : label_it->second;
    std::string label = state + "_" + std::to_string(lineage_counts[i]);
    double duration = durations[i];
    // Adjacent segments with the same state and lineage support are one public
    // scenario run, so only their duration changes.
    if (!out.labels.empty() && out.labels.back() == label) {
      out.durations.back() += duration;
    } else {
      // A changed state or lineage count starts a new scenario-map segment.
      out.labels.push_back(label);
      out.durations.push_back(duration);
    }
  }
  return out;
}

/**
 * Convert an internal transition map into a named R numeric vector.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static Rcpp::NumericVector pst_v34_transition_map_to_vector(const V34TransitionMap& edge_map) {
  Rcpp::NumericVector values(edge_map.durations.size());
  Rcpp::CharacterVector names(edge_map.labels.size());
  // Copy each aligned label/duration pair into the named-vector representation
  // expected by an R `simmap` edge.
  for (int i = 0; i < static_cast<int>(edge_map.durations.size()); ++i) {
    values[i] = edge_map.durations[i];
    names[i] = edge_map.labels[i];
  }
  values.attr("names") = names;
  return values;
}

/**
 * Convert an internal transition map into a named R numeric vector.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static Rcpp::NumericVector pst_v34_scenario_edge_map_to_vector(const V34TransitionMap& edge_map) {
  return pst_v34_transition_map_to_vector(edge_map);
}

/**
 * Create an internal scenario-map segment vector from states, sizes, and durations.
 *
 * Inputs are read only and failures stop before a partially formatted public tree is returned.
 */
static std::pair<V34TransitionMap, V34TransitionMap> pst_v34_make_scenario_maps(
    const Rcpp::IntegerVector& state_ids,
    const Rcpp::IntegerVector& path_ids,
    const Rcpp::IntegerVector& sizes,
    const Rcpp::NumericVector& durations,
    const std::map<int, std::string>& state_labels,
    const std::map<int, std::string>& path_labels) {
  V34TransitionMap state_out;
  V34TransitionMap path_out;
  int n = durations.size();
  if (state_ids.size() != n ||
      path_ids.size() != n ||
      sizes.size() != n ||
      n < 1) {
    Rcpp::stop("V34 scenario record vectors are not nonempty and parallel");
  }
  // Convert the committed state, size, and duration vectors in their original
  // segment order so topology and scenario support remain aligned.
  for (int i = 0; i < n; ++i) {
    double duration = durations[i];
    if (!R_finite(duration) || duration < 0.0) {
      Rcpp::stop("V34 scenario segment duration is nonfinite or negative");
    }
    int state_id = state_ids[i];
    int path_id = path_ids[i];
    int size = sizes[i];
    auto label_it = state_labels.find(state_id);
    std::string state = label_it == state_labels.end() ? std::string("NA") : label_it->second;
    std::string label = state + "_" + std::to_string(size);
    std::map<int, std::string>::const_iterator path_it = path_labels.find(path_id);
    if (path_it == path_labels.end()) {
      Rcpp::stop("V34 scenario segment contains an unregistered path id");
    }
    const std::string& path_label = path_it->second;
    state_out.labels.push_back(label);
    state_out.durations.push_back(duration);
    path_out.labels.push_back(path_label);
    path_out.durations.push_back(duration);
  }
  return std::make_pair(state_out, path_out);
}

/**
 * Parse, count, or summarize path-history labels for transition/path statistics.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static std::string pst_v34_escape_path_state(const std::string& state) {
  std::string escaped;
  escaped.reserve(state.size());
  // Escape every byte that could be confused with the public path delimiter;
  // all other state-label bytes are preserved verbatim.
  for (char ch : state) {
    // Percent must be escaped first so later delimiter escapes remain reversible.
    if (ch == '%') {
      escaped += "%25";
    } else if (ch == '|') {
      // Escape the path delimiter so one biological state cannot become two path steps.
      escaped += "%7C";
    } else {
      // Ordinary state-label bytes are copied unchanged.
      escaped.push_back(ch);
    }
  }
  return escaped;
}

/**
 * Format a root-to-state path label from state-label components.
 *
 * Inputs are read only and failures stop before a partially formatted public tree is returned.
 */
static std::string pst_v34_path_label(
    const std::vector<std::string>& components,
    const std::string& sep) {
  std::string out = sep;
  // Append root-to-current state components in traversal order; the resulting
  // label is the canonical key shared by PST, TT, and SS path surfaces.
  for (int i = 0; i < static_cast<int>(components.size()); ++i) {
    // Every step after the first needs an additional separator between states.
    if (i > 0) {
      out += sep;
    }
    out += pst_v34_escape_path_state(components[i]);
  }
  return out;
}

/**
 * Return phylogeny edge ids in stable root-first order for tree serialization.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static std::vector<int> pst_v34_preorder_edges(const Rcpp::IntegerMatrix& edge, int root) {
  int nedge = edge.nrow();
  int max_node = root;
  // Find the node-id domain before allocating the parent-to-child-edge index.
  for (int i = 0; i < nedge; ++i) {
    max_node = std::max(max_node, edge(i, 0));
    max_node = std::max(max_node, edge(i, 1));
  }

  std::vector<std::vector<int> > children_by_parent(max_node + 1);
  // Index every phylogeny edge under its parent while retaining 1-based public edge ids.
  for (int edge_id = 1; edge_id <= nedge; ++edge_id) {
    int parent = edge(edge_id - 1, 0);
    // Only parents inside the allocated node domain can own an outgoing edge;
    // malformed parents are detected later when not all edges are reachable.
    if (parent >= 0 && parent <= max_node) {
      children_by_parent[parent].push_back(edge_id);
    }
  }
  // Sort each sibling set by child id so serialization is independent of input edge order.
  for (std::vector<int>& edges : children_by_parent) {
    std::sort(edges.begin(), edges.end(), [&](int left, int right) {
      return edge(left - 1, 1) < edge(right - 1, 1);
    });
  }

  std::vector<int> queue(1, root);
  std::vector<int> ordered;
  // Expand the root-first node queue, appending each outgoing edge before its descendants.
  for (int pos = 0; pos < static_cast<int>(queue.size()); ++pos) {
    int parent = queue[pos];
    // Skip a malformed queued node; the reachability check below converts any
    // resulting missing edge into an explicit topology error.
    if (parent < 0 || parent >= static_cast<int>(children_by_parent.size())) {
      continue;
    }
    const std::vector<int>& edges = children_by_parent[parent];
    // Emit this parent's sorted child edges and queue their child nodes.
    for (int edge_id : edges) {
      ordered.push_back(edge_id);
      queue.push_back(edge(edge_id - 1, 1));
    }
  }
  // A rooted phylogeny must make every edge reachable exactly once; otherwise
  // deterministic PST serialization would silently omit topology.
  if (static_cast<int>(ordered.size()) != nedge) {
    Rcpp::stop("Could not traverse every simmap edge from the root");
  }
  return ordered;
}

/**
 * Split an encoded path or label string into stable components.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static std::vector<std::string> pst_v34_split_csv(const std::string& label) {
  std::vector<std::string> out;
  std::string current;
  // Scan the encoded descendant label once, preserving component order.
  for (char ch : label) {
    // A comma closes the current represented-tip component.
    if (ch == ',') {
      // Empty components carry no tip identity and are intentionally discarded.
      if (!current.empty()) {
        out.push_back(current);
      }
      current.clear();
    } else {
      // Non-delimiter bytes belong to the current tip label.
      current.push_back(ch);
    }
  }
  // Flush the final component because a valid encoded label need not end in a comma.
  if (!current.empty()) {
    out.push_back(current);
  }
  return out;
}

/**
 * Join label components into the encoded public label form.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static std::string pst_v34_join_csv(const std::vector<std::string>& labels) {
  std::ostringstream out;
  // Re-encode represented biological tips in their established public order.
  for (int i = 0; i < static_cast<int>(labels.size()); ++i) {
    // Separators occur only between labels, never before the first descendant.
    if (i > 0) {
      out << ",";
    }
    out << labels[i];
  }
  return out.str();
}

/**
 * Parse or format scenario labels of the form state plus lineage count.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static std::string pst_v34_scenario_state(const std::string& label) {
  std::size_t pos = label.find_last_of('_');
  // Labels without a size suffix are already state-only and remain unchanged.
  if (pos == std::string::npos) {
    return label;
  }
  return label.substr(0, pos);
}

/**
 * Parse or format scenario labels of the form state plus lineage count.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static int pst_v34_scenario_size(const std::string& label) {
  std::size_t pos = label.find_last_of('_');
  // A missing or empty suffix cannot provide lineage support; return `NA` so
  // callers cannot mistake the malformed label for a real zero-size scenario.
  if (pos == std::string::npos || pos + 1 >= label.size()) {
    return NA_INTEGER;
  }
  return std::atoi(label.substr(pos + 1).c_str());
}

/**
 * Parse or format scenario labels of the form state plus lineage count.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static std::string pst_v34_scenario_label(const std::string& state, int size) {
  return state + "_" + std::to_string(size);
}

/**
 * Drop zero-duration entries and merge adjacent identical labels in a transition map.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static V34TransitionMap pst_v34_compact_transition_map(const V34TransitionMap& edge_map) {
  V34TransitionMap out;
  // Visit mapped segments in biological time order, retaining only positive finite exposure.
  for (int i = 0; i < static_cast<int>(edge_map.durations.size()); ++i) {
    double duration = edge_map.durations[i];
    // Nonfinite and zero-duration entries create no public branch exposure.
    if (!R_finite(duration) || duration == 0.0) {
      continue;
    }
    const std::string& label = edge_map.labels[i];
    // Identical adjacent labels are one continuous biological run, not a transition.
    if (!out.labels.empty() && out.labels.back() == label) {
      out.durations.back() += duration;
    } else {
      // A changed label begins a new mapped segment at the same chronological position.
      out.labels.push_back(label);
      out.durations.push_back(duration);
    }
  }
  return out;
}

/** Compact scenario and path maps jointly so segment boundaries cannot drift. */
static void pst_v34_compact_scenario_edge_maps(
    V34ScenarioEdge& edge,
    double time_tolerance) {
  if (edge.map.labels.size() != edge.path_map.labels.size() ||
      edge.map.durations.size() != edge.path_map.durations.size() ||
      edge.map.labels.size() != edge.map.durations.size()) {
    Rcpp::stop("V34 scenario state/path maps lost segment alignment");
  }
  V34TransitionMap state_out;
  V34TransitionMap path_out;
  // Retain the shared chronological boundary unless both identities match.
  for (int segment_id = 0;
       segment_id < static_cast<int>(edge.map.durations.size());
       ++segment_id) {
    double state_duration = edge.map.durations[segment_id];
    double path_duration = edge.path_map.durations[segment_id];
    if (!R_finite(state_duration) || !R_finite(path_duration) ||
        std::fabs(state_duration - path_duration) > time_tolerance) {
      Rcpp::stop("V34 scenario state/path durations diverged");
    }
    if (state_duration == 0.0) {
      continue;
    }
    bool same_identity =
      !state_out.labels.empty() &&
      state_out.labels.back() == edge.map.labels[segment_id] &&
      path_out.labels.back() == edge.path_map.labels[segment_id];
    if (same_identity) {
      state_out.durations.back() += state_duration;
      path_out.durations.back() += path_duration;
    } else {
      state_out.labels.push_back(edge.map.labels[segment_id]);
      state_out.durations.push_back(state_duration);
      path_out.labels.push_back(edge.path_map.labels[segment_id]);
      path_out.durations.push_back(path_duration);
    }
  }
  edge.map = state_out;
  edge.path_map = path_out;
}

/** Concatenate maps without compacting boundaries owned by a parallel map. */
static V34TransitionMap pst_v34_concat_maps_raw(
    const std::vector<V34TransitionMap>& maps) {
  V34TransitionMap out;
  for (const V34TransitionMap& map : maps) {
    out.labels.insert(out.labels.end(), map.labels.begin(), map.labels.end());
    out.durations.insert(
      out.durations.end(), map.durations.begin(), map.durations.end()
    );
  }
  return out;
}

/**
 * Sum one compact mapped edge without changing its segment boundaries.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static double pst_v34_map_sum(const V34TransitionMap& edge_map) {
  double total = 0.0;
  // Sum every retained segment to recover the physical edge length.
  for (double value : edge_map.durations) {
    total += value;
  }
  return total;
}

/**
 * Concatenate mapped-edge vectors while preserving segment order.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static V34TransitionMap pst_v34_concat_maps(const std::vector<V34TransitionMap>& maps) {
  V34TransitionMap out;
  // Append source maps in edge-path order before compacting their shared boundaries.
  for (const V34TransitionMap& edge_map : maps) {
    // Copy each label/duration pair together so map identity cannot drift from exposure.
    for (int i = 0; i < static_cast<int>(edge_map.durations.size()); ++i) {
      out.labels.push_back(edge_map.labels[i]);
      out.durations.push_back(edge_map.durations[i]);
    }
  }
  return pst_v34_compact_transition_map(out);
}

/**
 * Find the largest node id before topology compaction allocates replacements.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static int pst_v34_max_scenario_node(const V34ScenarioTree& tree) {
  int max_node = 0;
  // Inspect both endpoints of every scenario edge to bound future node allocation.
  for (const V34ScenarioEdge& edge : tree.edges) {
    max_node = std::max(max_node, edge.parent);
    max_node = std::max(max_node, edge.child);
  }
  return max_node;
}

/**
 * Restore ape-compatible compact node numbering after topology edits.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static void pst_v34_renumber_internal_nodes(V34ScenarioTree& tree) {
  int ntips = static_cast<int>(tree.tip_labels.size());
  std::set<int> all_nodes;
  std::vector<int> parents_order;
  std::set<int> parents_seen;
  std::set<int> children;
  // Collect the node domain and first parent order from every retained edge;
  // this order makes root selection deterministic after compaction.
  for (const V34ScenarioEdge& edge : tree.edges) {
    all_nodes.insert(edge.parent);
    all_nodes.insert(edge.child);
    children.insert(edge.child);
    // Record each internal parent once while preserving its first edge occurrence.
    if (!parents_seen.count(edge.parent)) {
      parents_order.push_back(edge.parent);
      parents_seen.insert(edge.parent);
    }
  }
  std::vector<int> internal_old;
  // Separate internal ids from the fixed public tip-id range `1..ntips`.
  for (int node : all_nodes) {
    // Any endpoint outside the tip range is an internal topology node.
    if (node < 1 || node > ntips) {
      internal_old.push_back(node);
    }
  }
  // A tree containing no internal endpoint needs no node remapping.
  if (internal_old.empty()) {
    tree.nnode = 0;
    return;
  }

  std::vector<int> root_candidates;
  // Root candidates are parents that never occur as a child.
  for (int parent : parents_order) {
    // Only a parent without an incoming edge can be the scenario-tree root.
    if (!children.count(parent)) {
      root_candidates.push_back(parent);
    }
  }
  std::sort(root_candidates.begin(), root_candidates.end());
  int root_old = root_candidates.empty() ? internal_old[0] : root_candidates[0];
  std::vector<int> other_internal;
  // Retain all nonroot internal nodes for contiguous numbering after the root.
  for (int node : internal_old) {
    // The root has the mandatory `ntips + 1` id and is excluded from this suffix.
    if (node != root_old) {
      other_internal.push_back(node);
    }
  }
  std::sort(other_internal.begin(), other_internal.end());

  std::map<int, int> node_map;
  // Public tip ids are stable and therefore map to themselves.
  for (int tip = 1; tip <= ntips; ++tip) {
    node_map[tip] = tip;
  }
  node_map[root_old] = ntips + 1;
  // Assign the remaining internal nodes contiguously after the root.
  for (int i = 0; i < static_cast<int>(other_internal.size()); ++i) {
    node_map[other_internal[i]] = ntips + 2 + i;
  }
  // Rewrite both endpoints together so edge maps remain attached to their original branches.
  for (V34ScenarioEdge& edge : tree.edges) {
    edge.parent = node_map[edge.parent];
    edge.child = node_map[edge.child];
  }
  tree.nnode = static_cast<int>(other_internal.size()) + 1;
}

/**
 * Restore ape-compatible compact node numbering after topology edits.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static void pst_v34_renumber_tree_nodes(V34ScenarioTree& tree) {
  std::vector<int> parents_order;
  std::set<int> parents_seen;
  std::vector<int> children_order;
  std::set<int> children_seen;
  std::set<int> all_nodes_set;
  // Discover first-occurrence parent/child order and the complete node domain
  // after topology edits that may have converted former internals into leaves.
  for (const V34ScenarioEdge& edge : tree.edges) {
    all_nodes_set.insert(edge.parent);
    all_nodes_set.insert(edge.child);
    // Preserve each parent's first occurrence for deterministic root selection.
    if (!parents_seen.count(edge.parent)) {
      parents_order.push_back(edge.parent);
      parents_seen.insert(edge.parent);
    }
    // Preserve each child's first occurrence for deterministic leaf ordering.
    if (!children_seen.count(edge.child)) {
      children_order.push_back(edge.child);
      children_seen.insert(edge.child);
    }
  }
  std::set<int> parent_set(parents_order.begin(), parents_order.end());
  std::vector<int> leaf_old;
  // Children that never parent another edge are the current tree leaves.
  for (int child : children_order) {
    // Exclude internal children because they still own outgoing topology.
    if (!parent_set.count(child)) {
      leaf_old.push_back(child);
    }
  }
  std::sort(leaf_old.begin(), leaf_old.end());
  std::set<int> leaf_set(leaf_old.begin(), leaf_old.end());
  std::vector<int> all_nodes(all_nodes_set.begin(), all_nodes_set.end());
  std::vector<int> internal_old;
  // Everything not classified as a leaf remains an internal node.
  for (int node : all_nodes) {
    // Internal nodes occupy the id range following all leaves.
    if (!leaf_set.count(node)) {
      internal_old.push_back(node);
    }
  }

  std::map<int, int> node_map;
  // Renumber leaves contiguously in deterministic old-id order.
  for (int i = 0; i < static_cast<int>(leaf_old.size()); ++i) {
    node_map[leaf_old[i]] = i + 1;
  }
  // Internal-node numbering is needed only when retained topology has a parent node.
  if (!internal_old.empty()) {
    std::set<int> children_set(children_order.begin(), children_order.end());
    std::vector<int> root_candidates;
    // Identify roots as parents without incoming edges.
    for (int parent : parents_order) {
      // A parent absent from the child set has no ancestor in this retained tree.
      if (!children_set.count(parent)) {
        root_candidates.push_back(parent);
      }
    }
    std::sort(root_candidates.begin(), root_candidates.end());
    int root_old = root_candidates.empty() ? internal_old[0] : root_candidates[0];
    node_map[root_old] = static_cast<int>(leaf_old.size()) + 1;
    std::vector<int> other_internal;
    // Collect nonroot internals for deterministic ids after `nleaves + 1`.
    for (int node : internal_old) {
      // Exclude the root because its required public id has already been assigned.
      if (node != root_old) {
        other_internal.push_back(node);
      }
    }
    std::sort(other_internal.begin(), other_internal.end());
    // Assign all remaining internal ids contiguously.
    for (int i = 0; i < static_cast<int>(other_internal.size()); ++i) {
      node_map[other_internal[i]] = static_cast<int>(leaf_old.size()) + 2 + i;
    }
    tree.nnode = static_cast<int>(other_internal.size()) + 1;
  } else {
    // A degenerate leaf-only tree reports no internal nodes.
    tree.nnode = 0;
  }

  std::vector<std::string> old_tip_labels = tree.tip_labels;
  std::vector<std::string> new_tip_labels;
  // Reorder tip labels to match the newly assigned leaf ids.
  for (int old_leaf : leaf_old) {
    // Existing biological leaves retain their represented-tip labels.
    if (old_leaf >= 1 && old_leaf <= static_cast<int>(old_tip_labels.size())) {
      new_tip_labels.push_back(old_tip_labels[old_leaf - 1]);
    } else {
      // A topology-created leaf has no original biological label; retain an
      // empty placeholder rather than reusing another tip's identity.
      new_tip_labels.push_back("");
    }
  }
  tree.tip_labels = new_tip_labels;
  // Apply the same node map to every edge endpoint without reordering edge maps.
  for (V34ScenarioEdge& edge : tree.edges) {
    edge.parent = node_map[edge.parent];
    edge.child = node_map[edge.child];
  }
}

/**
 * Convert an R simmap list into the internal C++ scenario-tree representation.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static V34ScenarioTree pst_v34_scenario_tree_from_r(
    Rcpp::List tree,
    double time_tolerance) {
  Rcpp::IntegerMatrix edge = tree["edge"];
  Rcpp::List maps = tree["maps"];
  bool has_path_maps = tree.containsElementNamed("path_maps");
  Rcpp::List path_maps = has_path_maps ?
    Rcpp::as<Rcpp::List>(tree["path_maps"]) : Rcpp::List();
  if (has_path_maps && path_maps.size() != maps.size()) {
    Rcpp::stop("V34 scenario state/path maps are not edge-parallel");
  }
  Rcpp::CharacterVector tip_labels = tree["tip.label"];
  Rcpp::IntegerVector raw_scenario_edge_ids =
    tree.containsElementNamed("raw_scenario_edge_ids") ?
      Rcpp::as<Rcpp::IntegerVector>(tree["raw_scenario_edge_ids"]) :
      Rcpp::IntegerVector();
  if (raw_scenario_edge_ids.size() > 0 &&
      raw_scenario_edge_ids.size() != edge.nrow()) {
    Rcpp::stop("V34 raw scenario-edge ids are not edge-parallel");
  }
  V34ScenarioTree out;
  out.has_path_maps = has_path_maps;
  out.nnode = Rcpp::as<int>(tree["Nnode"]);
  out.tip_labels.reserve(tip_labels.size());
  // Copy public tip labels into owned C++ strings before topology mutation.
  for (int i = 0; i < tip_labels.size(); ++i) {
    out.tip_labels.push_back(Rcpp::as<std::string>(tip_labels[i]));
  }
  out.edges.reserve(edge.nrow());
  // Import each edge with its segment-parallel map and compact only adjacent identical runs.
  for (int i = 0; i < edge.nrow(); ++i) {
    Rcpp::NumericVector edge_map = maps[i];
    Rcpp::CharacterVector names = edge_map.names();
    Rcpp::NumericVector edge_path_map = has_path_maps ?
      Rcpp::as<Rcpp::NumericVector>(path_maps[i]) : edge_map;
    Rcpp::CharacterVector path_names = edge_path_map.names();
    if (has_path_maps && edge_path_map.size() != edge_map.size()) {
      Rcpp::stop("V34 scenario state/path maps are not segment-parallel");
    }
    V34TransitionMap map;
    V34TransitionMap path_map;
    // Copy each named duration as one internal map segment.
    for (int j = 0; j < edge_map.size(); ++j) {
      map.labels.push_back(Rcpp::as<std::string>(names[j]));
      map.durations.push_back(edge_map[j]);
      path_map.labels.push_back(Rcpp::as<std::string>(path_names[j]));
      path_map.durations.push_back(edge_path_map[j]);
    }
    out.edges.push_back(V34ScenarioEdge{
      edge(i, 0),
      edge(i, 1),
      map,
      path_map,
      std::set<int>{
        raw_scenario_edge_ids.size() > 0 ?
          raw_scenario_edge_ids[i] : i + 1
      }
    });
    pst_v34_compact_scenario_edge_maps(
      out.edges.back(),
      time_tolerance
    );
  }
  return out;
}

/**
 * Remove a zero-duration terminal child when its parent has another child.
 *
 * Such a leaf is created only to represent a lineage ending in the same batch
 * that another lineage continues or changes path. Its traversal provenance is
 * transferred to the entering parent edge so terminal metadata can still name
 * the canonical scenario occupied at the endpoint.
 */
static void pst_v34_legacy_prune_endpoint_leaf_cpp(
    V34ScenarioTree& tree,
    double time_tolerance) {
  bool changed = true;
  while (changed) {
    changed = false;
    int ntips = static_cast<int>(tree.tip_labels.size());
    for (int edge_id = 0;
         edge_id < static_cast<int>(tree.edges.size());
         ++edge_id) {
      const V34ScenarioEdge& terminal_edge = tree.edges[edge_id];
      if (terminal_edge.child < 1 || terminal_edge.child > ntips ||
          pst_v34_map_sum(terminal_edge.map) > time_tolerance) {
        continue;
      }
      int sibling_count = 0;
      int parent_edge_id = -1;
      for (int candidate_id = 0;
           candidate_id < static_cast<int>(tree.edges.size());
           ++candidate_id) {
        if (tree.edges[candidate_id].parent == terminal_edge.parent) {
          ++sibling_count;
        }
        if (tree.edges[candidate_id].child == terminal_edge.parent) {
          parent_edge_id = candidate_id;
        }
      }
      if (sibling_count < 2 || parent_edge_id < 0) {
        continue;
      }
      tree.edges[parent_edge_id].raw_edge_ids.insert(
        terminal_edge.raw_edge_ids.begin(),
        terminal_edge.raw_edge_ids.end()
      );
      tree.edges.erase(tree.edges.begin() + edge_id);
      pst_v34_renumber_tree_nodes(tree);
      changed = true;
      break;
    }
  }
}

/**
 * Format lineage-size maps directly from native `STATE_SIZE` scenario labels.
 *
 * Each output vector remains segment-parallel to its scenario edge map. State
 * names are separated from integer lineage sizes only at this public boundary;
 * traversal and canonicalization continue to use the combined label as one
 * scenario identity. Malformed labels retain size zero rather than borrowing
 * support from another segment.
 */
static Rcpp::List pst_v34_scenario_lineage_size_maps(
    const V34ScenarioTree& tree) {
  Rcpp::List lineage_size_maps(tree.edges.size());
  // Convert every canonical scenario edge without changing edge or segment order.
  for (int edge_id = 0; edge_id < static_cast<int>(tree.edges.size()); ++edge_id) {
    const V34TransitionMap& edge_map = tree.edges[edge_id].map;
    Rcpp::NumericVector sizes(edge_map.labels.size());
    Rcpp::CharacterVector states(edge_map.labels.size());
    // Split each committed scenario label into the state key and displayed lineage size.
    for (int segment_id = 0;
         segment_id < static_cast<int>(edge_map.labels.size());
         ++segment_id) {
      const std::string& label = edge_map.labels[segment_id];
      states[segment_id] = pst_v34_scenario_state(label);
      sizes[segment_id] = pst_v34_scenario_size(label);
    }
    sizes.attr("names") = states;
    lineage_size_maps[edge_id] = sizes;
  }
  return lineage_size_maps;
}

/**
 * Convert the internal C++ scenario-tree representation back to an R simmap list.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static Rcpp::List pst_v34_scenario_tree_to_r(const V34ScenarioTree& tree) {
  int nedge = static_cast<int>(tree.edges.size());
  Rcpp::IntegerMatrix edge(nedge, 2);
  Rcpp::NumericVector edge_length(nedge);
  Rcpp::List maps(nedge);
  Rcpp::List path_maps(nedge);
  Rcpp::List raw_scenario_edge_ids(nedge);
  Rcpp::CharacterVector rownames(nedge);
  // Serialize every internal edge, preserving edge/map order and recomputing
  // physical edge length from its compact scenario map.
  for (int i = 0; i < nedge; ++i) {
    edge(i, 0) = tree.edges[i].parent;
    edge(i, 1) = tree.edges[i].child;
    edge_length[i] = pst_v34_map_sum(tree.edges[i].map);
    maps[i] = pst_v34_transition_map_to_vector(tree.edges[i].map);
    path_maps[i] = pst_v34_transition_map_to_vector(tree.edges[i].path_map);
    raw_scenario_edge_ids[i] = Rcpp::wrap(tree.edges[i].raw_edge_ids);
    rownames[i] = std::to_string(tree.edges[i].parent) + "," + std::to_string(tree.edges[i].child);
  }
  edge.attr("dimnames") = Rcpp::List::create(
    R_NilValue,
    Rcpp::CharacterVector::create("parent", "child")
  );
  Rcpp::CharacterVector tip_labels(tree.tip_labels.size());
  // Emit tip labels in the same id order used by serialized edge children.
  for (int i = 0; i < static_cast<int>(tree.tip_labels.size()); ++i) {
    SET_STRING_ELT(tip_labels, i, Rf_mkCharCE(tree.tip_labels[i].c_str(), CE_UTF8));
  }
  Rcpp::NumericMatrix mapped_edge(nedge, 0);
  mapped_edge.attr("dimnames") = Rcpp::List::create(rownames, Rcpp::CharacterVector(0));
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("edge") = edge,
    Rcpp::Named("edge.length") = edge_length,
    Rcpp::Named("tip.label") = tip_labels,
    Rcpp::Named("Nnode") = tree.nnode,
    Rcpp::Named("maps") = maps,
    Rcpp::Named("mapped.edge") = mapped_edge
  );
  if (tree.has_path_maps) {
    Rcpp::CharacterVector edge_type(nedge);
    std::map<int, int> edge_by_child;
    // Index final canonical edges by child node for parent-edge classification.
    for (int edge_id = 0; edge_id < nedge; ++edge_id) {
      edge_by_child[tree.edges[edge_id].child] = edge_id;
    }
    for (int edge_id = 0; edge_id < nedge; ++edge_id) {
      std::map<int, int>::const_iterator parent_edge =
        edge_by_child.find(tree.edges[edge_id].parent);
      if (parent_edge == edge_by_child.end()) {
        edge_type[edge_id] = "Root";
        continue;
      }
      const V34TransitionMap& parent_paths =
        tree.edges[parent_edge->second].path_map;
      const V34TransitionMap& child_paths = tree.edges[edge_id].path_map;
      if (parent_paths.labels.empty() || child_paths.labels.empty()) {
        Rcpp::stop("V34 final scenario edge lacks path identity for classification");
      }
      edge_type[edge_id] =
        parent_paths.labels.back() == child_paths.labels.front() ?
          "Stay" : "Leave";
    }
    out["edge_type"] = edge_type;
    out["path_maps"] = path_maps;
  }
  // Scenario lineage sizes are native formatting output, not a legacy R-side
  // reconstruction from the returned tree.
  out["size_maps"] = Rcpp::List::create(
    Rcpp::Named("lineage_size_maps") = pst_v34_scenario_lineage_size_maps(tree)
  );
  out["raw_scenario_edge_ids"] = raw_scenario_edge_ids;
  out.attr("pst_v34_raw_scenario_edge_ids") = raw_scenario_edge_ids;
  out.attr("class") = Rcpp::CharacterVector::create("simmap", "phylo");
  return out;
}

/**
 * Remove zero-length scenario-tree edges created by cleanup while preserving mapped durations.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static void pst_v34_contract_zero_length_edges_cpp(
    V34ScenarioTree& tree,
    double time_tolerance) {
  int ntips = static_cast<int>(tree.tip_labels.size());
  bool searching = true;
  // Contract one eligible zero-length internal bridge per pass; restarting
  // avoids stale vector indices after an edge erase.
  while (searching) {
    searching = false;
    // Search retained edges for a zero-exposure internal child that can be bypassed.
    for (int edge_index = 0; edge_index < static_cast<int>(tree.edges.size()); ++edge_index) {
      // Positive-length edges and terminal tips are biologically observable and cannot be contracted.
      if (pst_v34_map_sum(tree.edges[edge_index].map) > time_tolerance ||
          tree.edges[edge_index].child <= ntips) {
        continue;
      }
      int parent_node = tree.edges[edge_index].parent;
      int child_node = tree.edges[edge_index].child;
      std::vector<int> child_edges;
      // Find every outgoing branch owned by the bridge's child node.
      for (int i = 0; i < static_cast<int>(tree.edges.size()); ++i) {
        // Only direct child edges are rewired when the bridge node disappears.
        if (tree.edges[i].parent == child_node) {
          child_edges.push_back(i);
        }
      }
      // A zero-length internal node without descendants is malformed; stop
      // contraction rather than deleting unexplained topology.
      if (child_edges.empty()) {
        return;
      }
      V34TransitionMap bridge_map = tree.edges[edge_index].map;
      // A unary bridge contributes its map prefix to its sole continuation.
      if (child_edges.size() == 1) {
        int child_edge_id = child_edges[0];
        tree.edges[child_edge_id].map = pst_v34_concat_maps_raw(
          std::vector<V34TransitionMap>{bridge_map, tree.edges[child_edge_id].map}
        );
        tree.edges[child_edge_id].path_map = pst_v34_concat_maps_raw(
          std::vector<V34TransitionMap>{
            tree.edges[edge_index].path_map,
            tree.edges[child_edge_id].path_map
          }
        );
        tree.edges[child_edge_id].parent = parent_node;
        tree.edges[child_edge_id].raw_edge_ids.insert(
          tree.edges[edge_index].raw_edge_ids.begin(),
          tree.edges[edge_index].raw_edge_ids.end()
        );
      } else {
        // A branching zero-length bridge contributes no duration; reparent all
        // daughters while preserving their individual maps.
        for (int child_edge_id : child_edges) {
          tree.edges[child_edge_id].parent = parent_node;
          tree.edges[child_edge_id].raw_edge_ids.insert(
            tree.edges[edge_index].raw_edge_ids.begin(),
            tree.edges[edge_index].raw_edge_ids.end()
          );
        }
      }
      tree.edges.erase(tree.edges.begin() + edge_index);
      searching = true;
      break;
    }
  }
  // Recompact every surviving map after bridge concatenation, then restore ape node ids.
  for (V34ScenarioEdge& edge : tree.edges) {
    pst_v34_compact_scenario_edge_maps(edge, time_tolerance);
  }
  pst_v34_renumber_internal_nodes(tree);
}

/**
 * Remove unary internal nodes created by compaction and repair topology metadata.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static void pst_v34_contract_unary_internal_edges_cpp(
    V34ScenarioTree& tree,
    double time_tolerance) {
  int ntips = static_cast<int>(tree.tip_labels.size());
  int max_node = pst_v34_max_scenario_node(tree);
  std::vector<int> parents;
  std::vector<int> children;
  std::set<int> parent_set;
  std::set<int> child_set;
  // Build parent/child sets and degree inputs from the current scenario topology.
  for (const V34ScenarioEdge& edge : tree.edges) {
    parents.push_back(edge.parent);
    children.push_back(edge.child);
    parent_set.insert(edge.parent);
    child_set.insert(edge.child);
  }
  std::vector<int> roots;
  // Identify roots as parent nodes with no incoming edge.
  for (int parent : parent_set) {
    // Only a parent absent from the child set can anchor the retained tree.
    if (!child_set.count(parent)) {
      roots.push_back(parent);
    }
  }
  // Without a root the topology is not safely contractible; leave it unchanged for validation upstream.
  if (roots.empty()) {
    return;
  }
  int root_node = roots[0];
  std::vector<int> outdegree(max_node + 1, 0);
  std::vector<int> child_edge_for_parent(max_node + 1, -1);
  // Count outgoing edges and remember each unary node's sole child edge.
  for (int i = 0; i < static_cast<int>(tree.edges.size()); ++i) {
    outdegree[tree.edges[i].parent] += 1;
    child_edge_for_parent[tree.edges[i].parent] = i;
  }
  std::vector<bool> is_unary(max_node + 1, false);
  // Mark nonroot internal nodes with exactly one continuation for contraction.
  for (int node = ntips + 1; node <= max_node; ++node) {
    is_unary[node] = outdegree[node] == 1 && node != root_node;
  }
  std::vector<V34ScenarioEdge> new_edges;
  // Start output chains only at non-unary parents; unary-parent edges are consumed by their ancestor chain.
  for (int edge_id = 0; edge_id < static_cast<int>(tree.edges.size()); ++edge_id) {
    int parent = tree.edges[edge_id].parent;
    // Skip an edge whose parent will be absorbed into an earlier concatenated chain.
    if (parent >= 0 && parent < static_cast<int>(is_unary.size()) && is_unary[parent]) {
      continue;
    }
    int child = tree.edges[edge_id].child;
    std::vector<V34TransitionMap> edge_maps;
    std::vector<V34TransitionMap> edge_path_maps;
    std::set<int> raw_edge_ids = tree.edges[edge_id].raw_edge_ids;
    edge_maps.push_back(tree.edges[edge_id].map);
    edge_path_maps.push_back(tree.edges[edge_id].path_map);
    // Follow consecutive unary internal nodes, appending their maps in biological time order.
    while (child > ntips &&
           child < static_cast<int>(is_unary.size()) &&
           is_unary[child]) {
      int child_edge_id = child_edge_for_parent[child];
      // A marked unary node must own a child edge; stop safely if topology metadata disagrees.
      if (child_edge_id < 0) {
        break;
      }
      const V34ScenarioEdge& child_edge = tree.edges[child_edge_id];
      bool duplicate_endpoint_contribution = false;
      if (child_edge.child <= ntips &&
          pst_v34_map_sum(child_edge.map) <= time_tolerance &&
          !edge_maps.empty() &&
          !edge_maps.back().labels.empty() &&
          !child_edge.map.labels.empty() &&
          !edge_path_maps.empty() &&
          !edge_path_maps.back().labels.empty() &&
          !child_edge.path_map.labels.empty()) {
        // Prefix hoisting can leave one structural zero-duration terminal tip
        // after its equivalent terminal-only siblings are removed. Its label
        // carries one tip's private support, while the accumulated parent
        // already carries the complete survivor/endpoint support. Appending it
        // would invent a zero-duration support drop on the public edge. Omit
        // only when state and path are unchanged; a real endpoint transition
        // remains an explicit mapped step.
        duplicate_endpoint_contribution =
          pst_v34_scenario_state(edge_maps.back().labels.back()) ==
            pst_v34_scenario_state(child_edge.map.labels.front()) &&
          edge_path_maps.back().labels.back() ==
            child_edge.path_map.labels.front();
      }
      if (!duplicate_endpoint_contribution) {
        edge_maps.push_back(child_edge.map);
        edge_path_maps.push_back(child_edge.path_map);
      }
      raw_edge_ids.insert(
        child_edge.raw_edge_ids.begin(),
        child_edge.raw_edge_ids.end()
      );
      child = child_edge.child;
    }
    new_edges.push_back(V34ScenarioEdge{
      parent,
      child,
      pst_v34_concat_maps_raw(edge_maps),
      pst_v34_concat_maps_raw(edge_path_maps),
      raw_edge_ids
    });
    pst_v34_compact_scenario_edge_maps(
      new_edges.back(),
      time_tolerance
    );
  }
  tree.edges = new_edges;
  // Compact concatenated chains and restore valid public node numbering.
  for (V34ScenarioEdge& edge : tree.edges) {
    pst_v34_compact_scenario_edge_maps(edge, time_tolerance);
  }
  pst_v34_renumber_internal_nodes(tree);
}

/**
 * Move identical leading sibling segments onto the parent edge to compact scenario trees.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static void pst_v34_legacy_shift_common_prefix_cpp(
    V34ScenarioTree& tree,
    double time_tolerance) {
  bool changed = true;
  // Repeat until no sibling group shares a positive leading interval; each
  // pass may expose a new common prefix after the previous hoist.
  while (changed) {
    changed = false;
    std::vector<int> parent_nodes;
    std::set<int> parent_seen;
    // Collect each current parent once so sibling groups are evaluated independently.
    for (const V34ScenarioEdge& edge : tree.edges) {
      // Preserve first occurrence to make group iteration deterministic.
      if (!parent_seen.count(edge.parent)) {
        parent_nodes.push_back(edge.parent);
        parent_seen.insert(edge.parent);
      }
    }
    // Evaluate outgoing scenario edges for each parent without mixing unrelated lineages.
    for (int parent_node : parent_nodes) {
      std::map<std::string, std::vector<int> > groups;
      // Group daughters by their leading biological state; only equal states can share a prefix.
      for (int edge_id = 0; edge_id < static_cast<int>(tree.edges.size()); ++edge_id) {
        // Ignore unrelated edges and empty maps because neither contributes a daughter prefix here.
        if (tree.edges[edge_id].parent != parent_node || tree.edges[edge_id].map.labels.empty()) {
          continue;
        }
        if (tree.edges[edge_id].path_map.labels.empty()) {
          Rcpp::stop("V34 scenario sibling prefix lacks a canonical path");
        }
        groups[
          pst_v34_scenario_state(tree.edges[edge_id].map.labels[0]) + "\n" +
          tree.edges[edge_id].path_map.labels[0]
        ].push_back(edge_id);
      }
      // Test each same-state daughter group for a positive common leading duration.
      for (const auto& group_pair : groups) {
        const std::vector<int>& group_edges = group_pair.second;
        // A single daughter has no sibling-shared segment to hoist.
        if (group_edges.size() <= 1) {
          continue;
        }
        double common_duration = tree.edges[group_edges[0]].map.durations[0];
        int shared_size = 0;
        // The hoisted duration is the shortest leading run; its lineage size is
        // the sum of daughter scenario supports over that shared interval.
        for (int edge_id : group_edges) {
          common_duration = std::min(common_duration, tree.edges[edge_id].map.durations[0]);
          shared_size += pst_v34_scenario_size(tree.edges[edge_id].map.labels[0]);
        }
        // Nonpositive common exposure is not a biological branch segment.
        if (common_duration <= time_tolerance) {
          continue;
        }
        std::string shared_state = pst_v34_scenario_state(tree.edges[group_edges[0]].map.labels[0]);
        int new_node = pst_v34_max_scenario_node(tree) + 1;
        V34TransitionMap shared_map;
        shared_map.labels.push_back(pst_v34_scenario_label(shared_state, shared_size));
        shared_map.durations.push_back(common_duration);
        V34TransitionMap shared_path_map;
        shared_path_map.labels.push_back(
          tree.edges[group_edges[0]].path_map.labels[0]
        );
        shared_path_map.durations.push_back(common_duration);
        tree.edges.push_back(V34ScenarioEdge{
          parent_node,
          new_node,
          shared_map,
          shared_path_map,
          std::set<int>()
        });
        for (int edge_id : group_edges) {
          tree.edges.back().raw_edge_ids.insert(
            tree.edges[edge_id].raw_edge_ids.begin(),
            tree.edges[edge_id].raw_edge_ids.end()
          );
        }
        // Reparent every participating daughter under the new shared scenario
        // edge and subtract exactly the hoisted exposure from its private map.
        for (int edge_id : group_edges) {
          tree.edges[edge_id].parent = new_node;
          tree.edges[edge_id].map.durations[0] -= common_duration;
          tree.edges[edge_id].path_map.durations[0] -= common_duration;
          pst_v34_compact_scenario_edge_maps(
            tree.edges[edge_id],
            time_tolerance
          );
        }
        changed = true;
        break;
      }
      // Restart from the updated topology after one successful hoist; stored edge ids are now stale.
      if (changed) {
        break;
      }
    }
  }
  // Remove any zero remainder exposed by the final hoist and restore public node ids.
  for (V34ScenarioEdge& edge : tree.edges) {
    pst_v34_compact_scenario_edge_maps(edge, time_tolerance);
  }
  pst_v34_renumber_internal_nodes(tree);
}

/**
 * Test exact mapped-segment identity before sibling scenario histories merge.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static bool pst_v34_maps_identical(const V34TransitionMap& left, const V34TransitionMap& right) {
  // Different segment counts cannot represent the same complete scenario trajectory.
  if (left.labels.size() != right.labels.size() ||
      left.durations.size() != right.durations.size()) {
    return false;
  }
  // Compare every aligned scenario identity and exposure exactly; merging must
  // not erase a state, size, or duration distinction.
  for (int i = 0; i < static_cast<int>(left.labels.size()); ++i) {
    // Any label or duration difference makes the histories biologically distinct.
    if (left.labels[i] != right.labels[i] || left.durations[i] != right.durations[i]) {
      return false;
    }
  }
  return true;
}

/**
 * Multiply encoded scenario sizes when identical represented tips are merged.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static V34TransitionMap pst_v34_multiply_support_map(const V34TransitionMap& edge_map, int multiplicity) {
  // A single represented lineage or an empty trajectory needs no support scaling.
  if (multiplicity <= 1 || edge_map.labels.empty()) {
    return edge_map;
  }
  V34TransitionMap out = edge_map;
  // Multiply each segment's lineage count because all merged tips share the full mapped history.
  for (std::string& label : out.labels) {
    std::string state = pst_v34_scenario_state(label);
    int size = pst_v34_scenario_size(label);
    label = pst_v34_scenario_label(state, size * multiplicity);
  }
  return out;
}

/**
 * Merge sibling tips with identical maps while preserving represented-tip multiplicity.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static bool pst_v34_legacy_collapse_peer_leaves_cpp(V34ScenarioTree& tree) {
  bool any_changed = false;
  // Merge one identical sibling-tip group per pass; node and edge ids are
  // rebuilt after each merge, so candidates must then be rediscovered.
  while (true) {
    int ntips = static_cast<int>(tree.tip_labels.size());
    std::vector<int> tip_edges;
    // Collect edges that terminate directly at current public tips.
    for (int i = 0; i < static_cast<int>(tree.edges.size()); ++i) {
      // Only terminal children participate in sibling-tip merging.
      if (tree.edges[i].child <= ntips) {
        tip_edges.push_back(i);
      }
    }
    // Fewer than two terminal edges cannot form a merge group.
    if (tip_edges.size() <= 1) {
      return any_changed;
    }
    std::vector<std::string> current_tip_order;
    std::set<std::string> seen_tips;
    // Recover the global represented-tip order before any composite labels are rebuilt.
    for (const std::string& label : tree.tip_labels) {
      // Expand each possibly merged label into its biological tip identities.
      for (const std::string& tip : pst_v34_split_csv(label)) {
        // Preserve the first occurrence of each represented tip and suppress duplicates.
        if (!seen_tips.count(tip)) {
          current_tip_order.push_back(tip);
          seen_tips.insert(tip);
        }
      }
    }
    std::map<int, std::vector<int> > parent_groups;
    // Partition terminal edges by parent because only true siblings may collapse together.
    for (int edge_id : tip_edges) {
      parent_groups[tree.edges[edge_id].parent].push_back(edge_id);
    }
    bool changed = false;
    // Evaluate each sibling set independently for complete mapped-history identity.
    for (const auto& group_pair : parent_groups) {
      std::vector<int> remaining = group_pair.second;
      // Partition remaining siblings by exact map identity until every candidate is classified.
      while (!remaining.empty()) {
        int reference_edge = remaining[0];
        std::vector<int> matching_edges;
        std::vector<int> not_matching;
        // Compare each sibling with the current reference trajectory.
        for (int candidate_edge : remaining) {
          // Exact state/size/duration equality places the edge in this merge class.
          if (pst_v34_maps_identical(
                tree.edges[reference_edge].map,
                tree.edges[candidate_edge].map
              ) &&
              pst_v34_maps_identical(
                tree.edges[reference_edge].path_map,
                tree.edges[candidate_edge].path_map
              )) {
            matching_edges.push_back(candidate_edge);
          } else {
            // Distinct histories remain available for a later comparison class.
            not_matching.push_back(candidate_edge);
          }
        }
        // A class with multiple siblings represents one public scenario trajectory.
        if (matching_edges.size() > 1) {
          int keep_edge = matching_edges[0];
          int keep_tip = tree.edges[keep_edge].child;
          std::vector<std::string> merged_tips;
          std::set<std::string> merged_seen;
          // Gather every biological tip represented by the equivalent sibling edges.
          for (int edge_id : matching_edges) {
            int tip_node = tree.edges[edge_id].child;
            // Skip malformed tip ids rather than attaching an unrelated label to the merged edge.
            if (tip_node < 1 || tip_node > static_cast<int>(tree.tip_labels.size())) {
              continue;
            }
            // Expand any already merged terminal label into biological tip identities.
            for (const std::string& tip : pst_v34_split_csv(tree.tip_labels[tip_node - 1])) {
              // Each represented tip appears once in the new composite label.
              if (!merged_seen.count(tip)) {
                merged_tips.push_back(tip);
                merged_seen.insert(tip);
              }
            }
          }
          std::map<std::string, int> order;
          // Build a rank map so the merged label retains original public tip order.
          for (int i = 0; i < static_cast<int>(current_tip_order.size()); ++i) {
            order[current_tip_order[i]] = i;
          }
          std::stable_sort(merged_tips.begin(), merged_tips.end(), [&](const std::string& a, const std::string& b) {
            return order[a] < order[b];
          });
          // Replace only a valid retained tip label with the ordered composite identity.
          if (keep_tip >= 1 && keep_tip <= static_cast<int>(tree.tip_labels.size())) {
            tree.tip_labels[keep_tip - 1] = pst_v34_join_csv(merged_tips);
          }
          tree.edges[keep_edge].map = pst_v34_multiply_support_map(
            tree.edges[keep_edge].map,
            static_cast<int>(matching_edges.size())
          );
          for (int edge_id : matching_edges) {
            if (edge_id == keep_edge) {
              continue;
            }
            tree.edges[keep_edge].raw_edge_ids.insert(
              tree.edges[edge_id].raw_edge_ids.begin(),
              tree.edges[edge_id].raw_edge_ids.end()
            );
          }
          std::set<int> remove_set(matching_edges.begin() + 1, matching_edges.end());
          std::vector<V34ScenarioEdge> kept_edges;
          // Remove redundant sibling edges while retaining every unrelated branch and map.
          for (int i = 0; i < static_cast<int>(tree.edges.size()); ++i) {
            // Keep the representative edge and all edges outside this merge class.
            if (!remove_set.count(i)) {
              kept_edges.push_back(tree.edges[i]);
            }
          }
          tree.edges = kept_edges;
          pst_v34_renumber_tree_nodes(tree);
          changed = true;
          any_changed = true;
          break;
        }
        remaining = not_matching;
      }
      // Stop using parent-group indices after topology has been renumbered.
      if (changed) {
        break;
      }
    }
    // No merge in a full pass means the sibling-tip topology is canonical.
    if (!changed) {
      return any_changed;
    }
  }
}

/**
 * Merge terminal trajectories that encode the same public scenario history.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static bool pst_v34_legacy_collapse_endpoint_histories_cpp(
    V34ScenarioTree& tree,
    const std::map<std::string, int>& trajectory_group_by_tip,
    double time_tolerance) {
  bool any_changed = false;
  // Merge one duplicate root-to-tip scenario trajectory per pass. Topology is
  // renumbered after each merge, so parent/child indexes are rebuilt each time.
  while (true) {
    int ntips = static_cast<int>(tree.tip_labels.size());
    std::map<int, int> parent_by_child;
    std::map<int, int> edge_by_child;
    // Index each edge by child for topology-only shared-prefix reconstruction.
    for (int edge_id = 0; edge_id < static_cast<int>(tree.edges.size()); ++edge_id) {
      parent_by_child[tree.edges[edge_id].child] = tree.edges[edge_id].parent;
      edge_by_child[tree.edges[edge_id].child] = edge_id;
    }

    std::map<int, std::vector<int> > leaves_by_trajectory;
    // Assign every current leaf to the traversal-owned terminal trajectory
    // group shared by all biological tips represented in its public label.
    for (int leaf = 1; leaf <= ntips; ++leaf) {
      auto edge_it = edge_by_child.find(leaf);
      // A leaf without an incoming edge is not a valid terminal trajectory candidate.
      if (edge_it == edge_by_child.end()) {
        continue;
      }
      std::vector<std::string> represented_tips = pst_v34_split_csv(
        tree.tip_labels[leaf - 1]
      );
      int trajectory_group = 0;
      bool group_is_consistent = !represented_tips.empty();
      // Verify every represented biological tip carries the same committed
      // trajectory group; sibling merges may already have combined labels.
      for (const std::string& tip : represented_tips) {
        std::map<std::string, int>::const_iterator group =
          trajectory_group_by_tip.find(tip);
        // Unknown tips or disagreeing groups make this leaf ineligible for a
        // broader whole-trajectory merge.
        if (group == trajectory_group_by_tip.end() ||
            (trajectory_group != 0 && trajectory_group != group->second)) {
          group_is_consistent = false;
          break;
        }
        trajectory_group = group->second;
      }
      // Ineligible leaves receive a unique negative key and therefore cannot
      // merge accidentally; positive keys are traversal-owned equivalence ids.
      int group_key = group_is_consistent && trajectory_group > 0 ?
        trajectory_group : -leaf;
      leaves_by_trajectory[group_key].push_back(leaf);
    }

    std::vector<int> merge_leaves;
    // Select the first deterministic trajectory class containing multiple leaves.
    for (const auto& item : leaves_by_trajectory) {
      // Only duplicate trajectory classes require a topology merge.
      if (item.second.size() > 1) {
        merge_leaves = item.second;
        break;
      }
    }
    // No duplicate class means terminal trajectory canonicalization is complete.
    if (merge_leaves.empty()) {
      return any_changed;
    }

    std::vector<std::string> current_tip_order;
    std::set<std::string> seen_tips;
    // Recover global biological tip order before composing the retained leaf label.
    for (const std::string& label : tree.tip_labels) {
      // Expand each current public label into represented biological tips.
      for (const std::string& tip : pst_v34_split_csv(label)) {
        // Preserve each tip's first public occurrence and suppress duplicates.
        if (!seen_tips.count(tip)) {
          current_tip_order.push_back(tip);
          seen_tips.insert(tip);
        }
      }
    }
    std::map<std::string, int> order;
    // Build stable tip ranks for the merged composite label.
    for (int i = 0; i < static_cast<int>(current_tip_order.size()); ++i) {
      order[current_tip_order[i]] = i;
    }

    int keep_leaf = merge_leaves[0];
    int keep_tip_multiplicity =
      keep_leaf >= 1 && keep_leaf <= static_cast<int>(tree.tip_labels.size()) ?
        static_cast<int>(pst_v34_split_csv(tree.tip_labels[keep_leaf - 1]).size()) :
        0;
    std::vector<std::string> merged_tips;
    std::set<std::string> merged_seen;
    // Gather represented tips from every leaf in the duplicate trajectory class.
    for (int leaf : merge_leaves) {
      // Ignore malformed leaf ids instead of reusing another leaf's biological identity.
      if (leaf < 1 || leaf > static_cast<int>(tree.tip_labels.size())) {
        continue;
      }
      // Expand previously merged labels before constructing the new union.
      for (const std::string& tip : pst_v34_split_csv(tree.tip_labels[leaf - 1])) {
        // Include each biological tip exactly once.
        if (!merged_seen.count(tip)) {
          merged_tips.push_back(tip);
          merged_seen.insert(tip);
        }
      }
    }
    std::stable_sort(merged_tips.begin(), merged_tips.end(), [&](const std::string& a, const std::string& b) {
      return order[a] < order[b];
    });
    // Store the ordered composite label only on the valid representative leaf.
    if (keep_leaf >= 1 && keep_leaf <= static_cast<int>(tree.tip_labels.size())) {
      tree.tip_labels[keep_leaf - 1] = pst_v34_join_csv(merged_tips);
    }

    std::vector<std::vector<int> > merge_paths;
    // Reconstruct each duplicate leaf's root-to-edge-id path for shared-prefix analysis.
    for (int leaf : merge_leaves) {
      std::vector<int> path_reverse;
      int node = leaf;
      // Walk from leaf to root through the child index, then reverse to root-first order.
      while (edge_by_child.count(node)) {
        int edge_id = edge_by_child[node];
        path_reverse.push_back(edge_id);
        node = parent_by_child[node];
      }
      std::reverse(path_reverse.begin(), path_reverse.end());
      merge_paths.push_back(path_reverse);
    }
    int shared_prefix = 0;
    // Shared-prefix comparison is meaningful only when at least one path was reconstructed.
    if (!merge_paths.empty()) {
      bool prefix_matches = true;
      // Advance while every trajectory owns the same edge at this root-first position.
      while (prefix_matches && shared_prefix < static_cast<int>(merge_paths[0].size())) {
        int edge_id = merge_paths[0][shared_prefix];
        // Compare the reference edge with the corresponding edge in every other path.
        for (int path_id = 1; path_id < static_cast<int>(merge_paths.size()); ++path_id) {
          // A shorter path or different edge ends the topology-shared prefix.
          if (shared_prefix >= static_cast<int>(merge_paths[path_id].size()) ||
              merge_paths[path_id][shared_prefix] != edge_id) {
            prefix_matches = false;
            break;
          }
        }
        // Count this position only when all duplicate trajectories share it.
        if (prefix_matches) {
          ++shared_prefix;
        }
      }
    }
    // Scale only private suffix edges; shared-prefix support already represents all descendants.
    if (!merge_paths.empty()) {
      const std::vector<int>& keep_path = merge_paths[0];
      bool aligned_private_support = true;
      int private_edge_count = static_cast<int>(keep_path.size()) - shared_prefix;
      for (const std::vector<int>& merge_path : merge_paths) {
        if (static_cast<int>(merge_path.size()) - shared_prefix !=
            private_edge_count) {
          aligned_private_support = false;
          break;
        }
        for (int offset = 0;
             aligned_private_support && offset < private_edge_count;
             ++offset) {
          int keep_edge_id = keep_path[shared_prefix + offset];
          int merge_edge_id = merge_path[shared_prefix + offset];
          if (keep_edge_id < 0 ||
              keep_edge_id >= static_cast<int>(tree.edges.size()) ||
              merge_edge_id < 0 ||
              merge_edge_id >= static_cast<int>(tree.edges.size())) {
            aligned_private_support = false;
            break;
          }
          const V34ScenarioEdge& keep_edge = tree.edges[keep_edge_id];
          const V34ScenarioEdge& merge_edge = tree.edges[merge_edge_id];
          if (keep_edge.map.labels.size() != merge_edge.map.labels.size() ||
              keep_edge.map.durations != merge_edge.map.durations ||
              !pst_v34_maps_identical(
                keep_edge.path_map,
                merge_edge.path_map
              )) {
            aligned_private_support = false;
            break;
          }
          for (int segment_id = 0;
               segment_id < static_cast<int>(keep_edge.map.labels.size());
               ++segment_id) {
            if (pst_v34_scenario_state(
                  keep_edge.map.labels[segment_id]
                ) !=
                pst_v34_scenario_state(
                  merge_edge.map.labels[segment_id]
                )) {
              aligned_private_support = false;
              break;
            }
          }
        }
      }

      // Equal-time duplicate suffixes carry parallel committed size runs. Sum
      // those runs segment by segment instead of applying one multiplicity to
      // a possibly unequal-terminal history. Different timing or path layouts
      // deliberately retain the representative's committed survivor support.
      if (aligned_private_support && merge_paths.size() > 1) {
        for (int offset = 0; offset < private_edge_count; ++offset) {
          int keep_edge_id = keep_path[shared_prefix + offset];
          V34ScenarioEdge& keep_edge = tree.edges[keep_edge_id];
          for (int segment_id = 0;
               segment_id < static_cast<int>(keep_edge.map.labels.size());
               ++segment_id) {
            int combined_size = 0;
            for (const std::vector<int>& merge_path : merge_paths) {
              int merge_edge_id = merge_path[shared_prefix + offset];
              combined_size += pst_v34_scenario_size(
                tree.edges[merge_edge_id].map.labels[segment_id]
              );
            }
            keep_edge.map.labels[segment_id] = pst_v34_scenario_label(
              pst_v34_scenario_state(keep_edge.map.labels[segment_id]),
              combined_size
            );
          }
        }
      }

      // Equal-endpoint duplicates can have different private edge layouts
      // after topology cleanup even though they cover the same complete
      // trajectory and duration. In that case the representative's committed
      // size run is per represented member group, so scale it by the exact
      // biological-tip ratio. Never apply this fallback to unequal terminal
      // durations: those histories require their existing piecewise support.
      if (!aligned_private_support && merge_paths.size() > 1 &&
          keep_tip_multiplicity > 0) {
        std::vector<double> private_durations;
        for (const std::vector<int>& merge_path : merge_paths) {
          double duration = 0.0;
          for (int path_index = shared_prefix;
               path_index < static_cast<int>(merge_path.size());
               ++path_index) {
            int edge_id = merge_path[path_index];
            if (edge_id < 0 || edge_id >= static_cast<int>(tree.edges.size())) {
              Rcpp::stop("V34 duplicate terminal trajectory has an invalid private edge");
            }
            duration += pst_v34_map_sum(tree.edges[edge_id].map);
          }
          private_durations.push_back(duration);
        }
        bool equal_terminal_duration = true;
        for (double duration : private_durations) {
          if (std::abs(duration - private_durations[0]) > time_tolerance) {
            equal_terminal_duration = false;
            break;
          }
        }
        if (equal_terminal_duration) {
          int merged_tip_multiplicity = static_cast<int>(merged_tips.size());
          for (int path_index = shared_prefix;
               path_index < static_cast<int>(keep_path.size());
               ++path_index) {
            V34ScenarioEdge& keep_edge = tree.edges[keep_path[path_index]];
            for (std::string& label : keep_edge.map.labels) {
              int original_size = pst_v34_scenario_size(label);
              int scaled_numerator = original_size * merged_tip_multiplicity;
              if (scaled_numerator % keep_tip_multiplicity != 0) {
                Rcpp::stop("V34 equal-time terminal support ratio is non-integral");
              }
              label = pst_v34_scenario_label(
                pst_v34_scenario_state(label),
                scaled_numerator / keep_tip_multiplicity
              );
            }
          }
        }
      }

      std::set<int> merged_private_raw_edge_ids;
      for (const std::vector<int>& merge_path : merge_paths) {
        for (int path_index = shared_prefix;
             path_index < static_cast<int>(merge_path.size());
             ++path_index) {
          int edge_id = merge_path[path_index];
          if (edge_id >= 0 && edge_id < static_cast<int>(tree.edges.size())) {
            merged_private_raw_edge_ids.insert(
              tree.edges[edge_id].raw_edge_ids.begin(),
              tree.edges[edge_id].raw_edge_ids.end()
            );
          }
        }
      }
      // Preserve traversal-committed support on the retained suffix. Group
      // members may terminate at different times, so no uniform multiplicity
      // can be applied without overcounting later private segments.
      for (int path_index = shared_prefix; path_index < static_cast<int>(keep_path.size()); ++path_index) {
        int edge_id = keep_path[path_index];
        // Retain valid edge ids from the pre-merge topology and merge only
        // their traversal provenance.
        if (edge_id >= 0 && edge_id < static_cast<int>(tree.edges.size())) {
          tree.edges[edge_id].raw_edge_ids.insert(
            merged_private_raw_edge_ids.begin(),
            merged_private_raw_edge_ids.end()
          );
        }
      }
    }

    std::set<int> remove_leaves(merge_leaves.begin() + 1, merge_leaves.end());
    std::vector<V34ScenarioEdge> kept_edges;
    // Drop redundant terminal edges and retain the representative plus all unrelated topology.
    for (const V34ScenarioEdge& edge : tree.edges) {
      // Remove only edges terminating at nonrepresentative duplicate leaves.
      if (!remove_leaves.count(edge.child)) {
        kept_edges.push_back(edge);
      }
    }
    tree.edges = kept_edges;
    pst_v34_contract_unary_internal_edges_cpp(tree, time_tolerance);
    pst_v34_renumber_tree_nodes(tree);
    any_changed = true;
  }
}

/**
 * Restore public tip-label order after C++ scenario-tree compaction merges tips.
 *
 * Inputs are read only and failures stop before a partially formatted public tree is returned.
 */
static void pst_v34_reorder_tip_labels_cpp(
    V34ScenarioTree& tree,
    const Rcpp::CharacterVector& original_tip_labels) {
  std::map<std::string, int> order;
  // Rank original biological labels once so composite scenario tips can be restored to input order.
  for (int i = 0; i < original_tip_labels.size(); ++i) {
    order[Rcpp::as<std::string>(original_tip_labels[i])] = i;
  }
  // Reorder every possibly composite scenario-tip label without changing its represented-tip set.
  for (std::string& label : tree.tip_labels) {
    std::vector<std::string> tips = pst_v34_split_csv(label);
    std::stable_sort(tips.begin(), tips.end(), [&](const std::string& a, const std::string& b) {
      auto ai = order.find(a);
      auto bi = order.find(b);
      int av = ai == order.end() ? std::numeric_limits<int>::max() : ai->second;
      int bv = bi == order.end() ? std::numeric_limits<int>::max() : bi->second;
      return av < bv;
    });
    label = pst_v34_join_csv(tips);
  }
}

struct V34TerminalEventRow {
  long long time_key;
  int phylo_tip_id;
  std::string phylo_tip_label;
  double terminal_time;
  int phylo_edge_id;
  int scenario_edge_id;
  int path_id;
  std::string state;
};

/** Serialize one public terminal-event row per input phylogenetic tip. */
static Rcpp::DataFrame pst_v34_terminal_events_from_tree(
    const V34ScenarioTree& tree,
    const Rcpp::List& summary,
    double time_tolerance) {
  Rcpp::DataFrame terminal_records = Rcpp::as<Rcpp::DataFrame>(
    summary["terminal_records"]
  );
  Rcpp::CharacterVector input_tip_labels = summary["tip_labels"];
  if (terminal_records.nrows() != input_tip_labels.size()) {
    Rcpp::stop("V34 terminal records are not aligned one-to-one with input tips");
  }
  Rcpp::IntegerVector tip_id = terminal_records["tip_id"];
  Rcpp::CharacterVector tip_label = terminal_records["tip_label"];
  Rcpp::IntegerVector phylo_edge_id = terminal_records["final_phylo_edge_id"];
  if (!terminal_records.containsElementNamed("stable_scenario_id")) {
    Rcpp::stop("V34 terminal records lack stable scenario identity");
  }
  Rcpp::IntegerVector stable_scenario_id =
    terminal_records["stable_scenario_id"];
  Rcpp::IntegerVector path_id = terminal_records["final_path_id"];
  Rcpp::NumericVector terminal_time = terminal_records["terminal_time"];

  Rcpp::DataFrame path_lookup = Rcpp::as<Rcpp::DataFrame>(
    summary["path_lookup"]
  );
  Rcpp::IntegerVector lookup_path_id = path_lookup["path_id"];
  Rcpp::IntegerVector terminal_state_id = path_lookup["terminal_state_id"];
  std::map<int, int> state_id_by_path;
  for (int row = 0; row < lookup_path_id.size(); ++row) {
    state_id_by_path[lookup_path_id[row]] = terminal_state_id[row];
  }
  Rcpp::DataFrame state_lookup = Rcpp::as<Rcpp::DataFrame>(
    summary["state_lookup"]
  );
  Rcpp::IntegerVector lookup_state_id = state_lookup["state_id"];
  Rcpp::CharacterVector lookup_state = state_lookup["state"];
  std::map<int, std::string> state_by_id;
  for (int row = 0; row < lookup_state_id.size(); ++row) {
    state_by_id[lookup_state_id[row]] =
      Rcpp::as<std::string>(lookup_state[row]);
  }

  int nedge = static_cast<int>(tree.edges.size());
  std::map<int, int> public_edge_id_by_stable_id;
  for (int edge_id = 0; edge_id < nedge; ++edge_id) {
    if (tree.edges[edge_id].raw_edge_ids.size() != 1) {
      Rcpp::stop(
        "V34 direct scenario edge lacks one exact stable identity"
      );
    }
    int stable_id = *tree.edges[edge_id].raw_edge_ids.begin();
    if (stable_id < 1 ||
        !public_edge_id_by_stable_id
          .insert(std::make_pair(stable_id, edge_id))
          .second) {
      Rcpp::stop("V34 direct scenario tree has duplicate stable identity");
    }
  }
  int max_node = pst_v34_max_scenario_node(tree);
  std::vector<int> edge_by_child(max_node + 1, -1);
  std::vector<std::vector<int> > outgoing(max_node + 1);
  for (int edge_id = 0; edge_id < nedge; ++edge_id) {
    int parent = tree.edges[edge_id].parent;
    int child = tree.edges[edge_id].child;
    if (parent < 1 || child < 1 || parent > max_node || child > max_node ||
        edge_by_child[child] != -1) {
      Rcpp::stop("V34 terminal-event topology is malformed");
    }
    edge_by_child[child] = edge_id;
    outgoing[parent].push_back(edge_id);
  }
  int root = -1;
  for (int node = 1; node <= max_node; ++node) {
    if (!outgoing[node].empty() && edge_by_child[node] == -1) {
      if (root != -1) {
        Rcpp::stop("V34 terminal-event topology has multiple roots");
      }
      root = node;
    }
  }
  if (root < 1) {
    Rcpp::stop("V34 terminal-event topology has no root");
  }

  std::vector<double> node_time(
    max_node + 1,
    std::numeric_limits<double>::quiet_NaN()
  );
  std::vector<int> queue(1, root);
  node_time[root] = 0.0;
  for (std::size_t queue_id = 0; queue_id < queue.size(); ++queue_id) {
    int parent = queue[queue_id];
    for (int edge_id : outgoing[parent]) {
      int child = tree.edges[edge_id].child;
      double length = pst_v34_map_sum(tree.edges[edge_id].map);
      if (!std::isfinite(length) || length < 0.0) {
        Rcpp::stop("V34 terminal-event scenario edge has invalid length");
      }
      node_time[child] = node_time[parent] + length;
      queue.push_back(child);
    }
  }

  double time_scale = pst_v34_batch_time_scale(time_tolerance);
  std::vector<long long> edge_start_key(nedge);
  std::vector<long long> edge_end_key(nedge);
  for (int edge_id = 0; edge_id < nedge; ++edge_id) {
    edge_start_key[edge_id] = pst_v34_batch_time_key(
      node_time[tree.edges[edge_id].parent],
      time_scale
    );
    edge_end_key[edge_id] = pst_v34_batch_time_key(
      node_time[tree.edges[edge_id].child],
      time_scale
    );
  }

  std::vector<V34TerminalEventRow> events;
  events.reserve(terminal_records.nrows());
  long long maximum_terminal_key = std::numeric_limits<long long>::min();
  for (int row = 0; row < terminal_records.nrows(); ++row) {
    if (tip_id[row] < 1 || tip_id[row] > input_tip_labels.size() ||
        phylo_edge_id[row] < 1 || stable_scenario_id[row] < 1 ||
        path_id[row] < 1 || !R_finite(terminal_time[row]) ||
        !state_id_by_path.count(path_id[row]) ||
        !state_by_id.count(state_id_by_path[path_id[row]])) {
      Rcpp::stop("V34 terminal record contains invalid endpoint identity");
    }
    long long time_key = pst_v34_batch_time_key(
      terminal_time[row],
      time_scale
    );
    std::map<int, int>::const_iterator public_edge =
      public_edge_id_by_stable_id.find(stable_scenario_id[row]);
    if (public_edge == public_edge_id_by_stable_id.end()) {
      std::ostringstream diagnostic;
      diagnostic
        << "V34 terminal record has no public edge for its stable identity"
        << " [tip=" << Rcpp::as<std::string>(tip_label[row])
        << ", stable_scenario_id=" << stable_scenario_id[row]
        << "]";
      Rcpp::stop(diagnostic.str());
    }
    int canonical_edge_id = public_edge->second;
    if (time_key < edge_start_key[canonical_edge_id] ||
        time_key > edge_end_key[canonical_edge_id]) {
      std::ostringstream diagnostic;
      diagnostic
        << "V34 terminal time lies outside its exact stable scenario edge"
        << " [tip=" << Rcpp::as<std::string>(tip_label[row])
        << ", stable_scenario_id=" << stable_scenario_id[row]
        << ", terminal_time_key=" << time_key
        << ", edge_time_key=" << edge_start_key[canonical_edge_id]
        << ".." << edge_end_key[canonical_edge_id]
        << "]";
      Rcpp::stop(diagnostic.str());
    }
    events.push_back(V34TerminalEventRow{
      time_key,
      tip_id[row],
      Rcpp::as<std::string>(tip_label[row]),
      terminal_time[row],
      phylo_edge_id[row],
      canonical_edge_id + 1,
      path_id[row],
      state_by_id[state_id_by_path[path_id[row]]]
    });
    maximum_terminal_key = std::max(maximum_terminal_key, time_key);
  }
  std::sort(events.begin(), events.end(), [](
      const V34TerminalEventRow& left,
      const V34TerminalEventRow& right) {
    if (left.time_key != right.time_key) {
      return left.time_key < right.time_key;
    }
    if (left.scenario_edge_id != right.scenario_edge_id) {
      return left.scenario_edge_id < right.scenario_edge_id;
    }
    return left.phylo_tip_id < right.phylo_tip_id;
  });

  int nevent = static_cast<int>(events.size());
  Rcpp::IntegerVector event_id(nevent);
  Rcpp::IntegerVector event_group_id(nevent);
  Rcpp::IntegerVector out_tip_id(nevent);
  Rcpp::CharacterVector out_tip_label(nevent);
  Rcpp::NumericVector out_terminal_time(nevent);
  Rcpp::NumericVector out_terminal_time_key(nevent);
  Rcpp::IntegerVector out_phylo_edge_id(nevent);
  Rcpp::IntegerVector out_scenario_edge_id(nevent);
  Rcpp::IntegerVector out_path_id(nevent);
  Rcpp::CharacterVector out_state(nevent);
  Rcpp::LogicalVector ends_before_horizon(nevent);
  int group_id = 0;
  for (int row = 0; row < nevent; ++row) {
    if (row == 0 || events[row].time_key != events[row - 1].time_key ||
        events[row].scenario_edge_id != events[row - 1].scenario_edge_id) {
      ++group_id;
    }
    event_id[row] = row + 1;
    event_group_id[row] = group_id;
    out_tip_id[row] = events[row].phylo_tip_id;
    out_tip_label[row] = events[row].phylo_tip_label;
    out_terminal_time[row] = events[row].terminal_time;
    out_terminal_time_key[row] =
      static_cast<double>(events[row].time_key);
    out_phylo_edge_id[row] = events[row].phylo_edge_id;
    out_scenario_edge_id[row] = events[row].scenario_edge_id;
    out_path_id[row] = events[row].path_id;
    out_state[row] = events[row].state;
    ends_before_horizon[row] = events[row].time_key < maximum_terminal_key;
  }
  Rcpp::DataFrame out = Rcpp::DataFrame::create(
    Rcpp::Named("event_id") = event_id,
    Rcpp::Named("event_group_id") = event_group_id,
    Rcpp::Named("phylo_tip_id") = out_tip_id,
    Rcpp::Named("phylo_tip_label") = out_tip_label,
    Rcpp::Named("terminal_time") = out_terminal_time,
    Rcpp::Named("phylo_edge_id") = out_phylo_edge_id,
    Rcpp::Named("scenario_edge_id") = out_scenario_edge_id,
    Rcpp::Named("path_id") = out_path_id,
    Rcpp::Named("state") = out_state,
    Rcpp::Named("ends_before_horizon") = ends_before_horizon,
    Rcpp::Named("stringsAsFactors") = false
  );
  out.attr("pst_v34_terminal_time_keys") = out_terminal_time_key;
  return out;
}

} // namespace

/**
 * Format the initial C++ scenario tree before canonical cleanup.
 *
 * This is an Rcpp boundary: inputs are R objects and outputs must preserve documented V34 public shapes.
 */
static Rcpp::List pst_v34_legacy_scenario_tree_initial_impl(
    Rcpp::List summary) {
  Rcpp::List records = summary["scenario_edge_records"];
  Rcpp::List edge_parent_nodes = records["edge_parent_nodes"];
  Rcpp::List edge_child_nodes = records["edge_child_nodes"];
  Rcpp::List edge_state_ids = records["edge_state_ids"];
  Rcpp::List edge_path_ids = records["edge_path_ids"];
  Rcpp::List edge_sizes = records["edge_sizes"];
  Rcpp::List edge_durations = records["edge_durations"];
  Rcpp::List edge_tip_ids = records["edge_tip_ids"];
  Rcpp::DataFrame state_lookup = Rcpp::as<Rcpp::DataFrame>(summary["state_lookup"]);
  Rcpp::IntegerVector state_id_column = state_lookup["state_id"];
  Rcpp::CharacterVector state_column = state_lookup["state"];
  Rcpp::CharacterVector tip_labels = summary["tip_labels"];
  Rcpp::DataFrame path_lookup = Rcpp::as<Rcpp::DataFrame>(summary["path_lookup"]);
  Rcpp::IntegerVector path_id_column = path_lookup["path_id"];
  Rcpp::CharacterVector path_label_column = path_lookup["label"];

  std::map<int, std::string> state_labels;
  // Cache state labels by traversal id for scenario-segment serialization.
  for (int i = 0; i < state_id_column.size(); ++i) {
    state_labels[state_id_column[i]] = Rcpp::as<std::string>(state_column[i]);
  }
  std::map<int, std::string> path_labels;
  // Canonical public path labels remain keyed by traversal path id.
  for (int i = 0; i < path_id_column.size(); ++i) {
    path_labels[path_id_column[i]] =
      Rcpp::as<std::string>(path_label_column[i]);
  }

  std::vector<int> edge_ids;
  int nrecord = edge_durations.size();
  // Retain only segment records with valid parent and child endpoints; empty
  // provisional records do not represent public scenario edges.
  for (int edge_id = 1; edge_id <= nrecord; ++edge_id) {
    Rcpp::IntegerVector parents = edge_parent_nodes[edge_id - 1];
    Rcpp::IntegerVector children = edge_child_nodes[edge_id - 1];
    // A record missing either endpoint cannot own a public branch.
    if (parents.size() == 0 || children.size() == 0) {
      continue;
    }
    int parent_tail = parents[parents.size() - 1];
    int child_tail = children[children.size() - 1];
    // Missing endpoint ids are excluded rather than creating an `NA` topology edge.
    if (parent_tail == NA_INTEGER || child_tail == NA_INTEGER) {
      continue;
    }
    edge_ids.push_back(edge_id);
  }

  int nedge = static_cast<int>(edge_ids.size());
  std::vector<std::pair<int, int> > edge_old;
  edge_old.reserve(nedge);
  std::set<int> all_old_set;
  std::vector<int> parent_order;
  std::set<int> parent_seen;
  std::vector<int> child_order;
  std::set<int> child_seen;
  // Collect valid old endpoints and deterministic first-occurrence order for node classification.
  for (int edge_id : edge_ids) {
    Rcpp::IntegerVector parents = edge_parent_nodes[edge_id - 1];
    Rcpp::IntegerVector children = edge_child_nodes[edge_id - 1];
    int parent_tail = parents[parents.size() - 1];
    int child_tail = children[children.size() - 1];
    edge_old.push_back(std::make_pair(parent_tail, child_tail));
    all_old_set.insert(parent_tail);
    all_old_set.insert(child_tail);
    // Record each parent once for deterministic root discovery.
    if (!parent_seen.count(parent_tail)) {
      parent_order.push_back(parent_tail);
      parent_seen.insert(parent_tail);
    }
    // Record each child once for deterministic leaf discovery.
    if (!child_seen.count(child_tail)) {
      child_order.push_back(child_tail);
      child_seen.insert(child_tail);
    }
  }

  std::set<int> parent_set(parent_order.begin(), parent_order.end());
  std::vector<int> leaf_old;
  // Scenario children that never parent another edge are public leaves.
  for (int child : child_order) {
    // Exclude children with outgoing edges because they remain internal nodes.
    if (!parent_set.count(child)) {
      leaf_old.push_back(child);
    }
  }
  std::vector<int> leaf_sorted = leaf_old;
  std::sort(leaf_sorted.begin(), leaf_sorted.end());
  std::set<int> leaf_set(leaf_old.begin(), leaf_old.end());
  std::vector<int> internal_old;
  // Classify every nonleaf endpoint as internal topology.
  for (int node : all_old_set) {
    // Internal ids will be remapped after all leaves.
    if (!leaf_set.count(node)) {
      internal_old.push_back(node);
    }
  }

  std::set<int> child_set(child_order.begin(), child_order.end());
  std::vector<int> root_candidates;
  // Root candidates are internal nodes with no incoming scenario edge.
  for (int node : internal_old) {
    // Only an internal node absent from the child set can anchor the tree.
    if (!child_set.count(node)) {
      root_candidates.push_back(node);
    }
  }
  int root_old = internal_old.empty() ? NA_INTEGER : internal_old[0];
  // A unique topology root supersedes the deterministic fallback internal id.
  if (root_candidates.size() == 1) {
    root_old = root_candidates[0];
  }

  std::map<int, int> node_map;
  // Assign contiguous public tip ids in sorted old-node order.
  for (int i = 0; i < static_cast<int>(leaf_sorted.size()); ++i) {
    node_map[leaf_sorted[i]] = i + 1;
  }
  int root_new = static_cast<int>(leaf_old.size()) + 1;
  // A nondegenerate scenario tree places its root immediately after all tips.
  if (root_old != NA_INTEGER) {
    node_map[root_old] = root_new;
  }
  std::vector<int> other_internal;
  // Collect nonroot internal ids for the remaining contiguous public range.
  for (int node : internal_old) {
    // Exclude the root because its public id is already fixed.
    if (node != root_old) {
      other_internal.push_back(node);
    }
  }
  std::sort(other_internal.begin(), other_internal.end());
  // Assign deterministic ids to all remaining internal nodes.
  for (int i = 0; i < static_cast<int>(other_internal.size()); ++i) {
    node_map[other_internal[i]] = static_cast<int>(leaf_old.size()) + 2 + i;
  }

  Rcpp::IntegerMatrix edge(nedge, 2);
  // Rewrite each retained endpoint pair under the public node map without changing edge order.
  for (int i = 0; i < nedge; ++i) {
    edge(i, 0) = node_map[edge_old[i].first];
    edge(i, 1) = node_map[edge_old[i].second];
  }
  edge.attr("dimnames") = Rcpp::List::create(
    R_NilValue,
    Rcpp::CharacterVector::create("parent", "child")
  );

  Rcpp::List maps(nedge);
  Rcpp::List path_maps(nedge);
  Rcpp::NumericVector edge_lengths(nedge);
  // Materialize each retained edge map from traversal-owned state, size, and duration records.
  for (int i = 0; i < nedge; ++i) {
    int edge_id = edge_ids[i];
    Rcpp::IntegerVector states = edge_state_ids[edge_id - 1];
    Rcpp::IntegerVector paths = edge_path_ids[edge_id - 1];
    Rcpp::IntegerVector sizes = edge_sizes[edge_id - 1];
    Rcpp::NumericVector durations = edge_durations[edge_id - 1];
    std::pair<V34TransitionMap, V34TransitionMap> edge_maps =
      pst_v34_make_scenario_maps(
      states,
      paths,
      sizes,
      durations,
      state_labels,
      path_labels
    );
    const V34TransitionMap& edge_map = edge_maps.first;
    maps[i] = pst_v34_scenario_edge_map_to_vector(edge_map);
    path_maps[i] = pst_v34_transition_map_to_vector(edge_maps.second);
    double total = 0.0;
    // Sum compact segments to produce the physical branch length paired with this map.
    for (double value : edge_map.durations) {
      total += value;
    }
    edge_lengths[i] = total;
  }

  std::map<int, std::string> tip_labels_old;
  // Initialize one represented-tip label slot for every public scenario leaf.
  for (int leaf : leaf_sorted) {
    tip_labels_old[leaf] = "";
  }
  // Attach biological tips using final memberships captured from the live
  // traversal frontier, not a completed scenario-id matrix.
  for (int i = 0; i < nedge; ++i) {
    int edge_id = edge_ids[i];
    std::vector<int> terminal_tip_ids = Rcpp::as<std::vector<int> >(
      edge_tip_ids[edge_id - 1]
    );
    // Only terminal scenario edges receive a represented-tip label.
    if (!terminal_tip_ids.empty()) {
      std::ostringstream label;
      // Join represented biological tips in input-tip order for deterministic labels.
      for (int j = 0; j < static_cast<int>(terminal_tip_ids.size()); ++j) {
        // Separate later represented tips from the first label.
        if (j > 0) {
          label << ",";
        }
        int tip_id = terminal_tip_ids[j];
        label << Rcpp::as<std::string>(tip_labels[tip_id - 1]);
      }
      tip_labels_old[edge_old[i].second] = label.str();
    }
  }

  Rcpp::CharacterVector final_tip_labels(leaf_sorted.size());
  // Emit scenario-tip labels in the same sorted leaf order used by the node map.
  for (int i = 0; i < static_cast<int>(leaf_sorted.size()); ++i) {
    const std::string& label = tip_labels_old[leaf_sorted[i]];
    SET_STRING_ELT(final_tip_labels, i, Rf_mkCharCE(label.c_str(), CE_UTF8));
  }

  int max_new_node = 0;
  // Bound the remapped node domain to derive the final `Nnode` value.
  for (int i = 0; i < nedge; ++i) {
    max_new_node = std::max(max_new_node, edge(i, 0));
    max_new_node = std::max(max_new_node, edge(i, 1));
  }

  Rcpp::NumericMatrix mapped_edge(nedge, 0);
  Rcpp::IntegerVector raw_scenario_edge_ids(edge_ids.begin(), edge_ids.end());
  Rcpp::List tree = Rcpp::List::create(
    Rcpp::Named("edge") = edge,
    Rcpp::Named("edge.length") = edge_lengths,
    Rcpp::Named("tip.label") = final_tip_labels,
    Rcpp::Named("Nnode") = max_new_node - static_cast<int>(leaf_old.size()),
    Rcpp::Named("maps") = maps,
    Rcpp::Named("mapped.edge") = mapped_edge,
    Rcpp::Named("path_maps") = path_maps,
    Rcpp::Named("raw_scenario_edge_ids") = raw_scenario_edge_ids
  );
  tree.attr("class") = Rcpp::CharacterVector::create("simmap", "phylo");
  return tree;
}

/**
 * Serialize the final traversal-owned scenario registry without changing its
 * topology, stable identities, or segment boundaries.
 */
Rcpp::List pst_v34_scenario_tree_initial_impl(Rcpp::List summary) {
  if (!summary.containsElementNamed("scenario_edge_records") ||
      !summary.containsElementNamed("state_lookup") ||
      !summary.containsElementNamed("path_lookup") ||
      !summary.containsElementNamed("tip_labels")) {
    Rcpp::stop("V34 direct scenario serialization lacks required summary fields");
  }
  Rcpp::List records = summary["scenario_edge_records"];
  const char* required_fields[] = {
    "edge_durations", "edge_end_times", "edge_state_ids", "edge_path_ids",
    "edge_sizes", "edge_tip_ids", "edge_parent_scenario_edge_id",
    "edge_formation_time", "edge_type"
  };
  for (const char* field : required_fields) {
    if (!records.containsElementNamed(field)) {
      Rcpp::stop(
        std::string("V34 direct scenario serialization lacks ") + field
      );
    }
  }

  Rcpp::List edge_durations = records["edge_durations"];
  Rcpp::List edge_end_times = records["edge_end_times"];
  Rcpp::List edge_state_ids = records["edge_state_ids"];
  Rcpp::List edge_path_ids = records["edge_path_ids"];
  Rcpp::List edge_sizes = records["edge_sizes"];
  Rcpp::List edge_tip_ids = records["edge_tip_ids"];
  Rcpp::IntegerVector parent_stable_ids =
    records["edge_parent_scenario_edge_id"];
  Rcpp::NumericVector formation_times = records["edge_formation_time"];
  Rcpp::CharacterVector edge_types = records["edge_type"];
  int nedge = edge_durations.size();
  if (nedge < 1 ||
      edge_end_times.size() != nedge ||
      edge_state_ids.size() != nedge ||
      edge_path_ids.size() != nedge ||
      edge_sizes.size() != nedge ||
      edge_tip_ids.size() != nedge ||
      parent_stable_ids.size() != nedge ||
      formation_times.size() != nedge ||
      edge_types.size() != nedge) {
    Rcpp::stop("V34 direct scenario records are not nonempty and edge-parallel");
  }

  std::set<int> referenced_as_parent;
  int root_record_count = 0;
  for (int edge_index = 0; edge_index < nedge; ++edge_index) {
    int stable_id = edge_index + 1;
    int parent_id = parent_stable_ids[edge_index];
    if (parent_id == 0) {
      ++root_record_count;
    } else {
      if (parent_id < 1 || parent_id > nedge || parent_id == stable_id) {
        std::ostringstream diagnostic;
        diagnostic
          << "V34 direct scenario record has an invalid parent stable ID"
          << " [stable_id=" << stable_id
          << ", parent_stable_id=" << parent_id
          << ", record_count=" << nedge
          << "]";
        Rcpp::stop(diagnostic.str());
      }
      referenced_as_parent.insert(parent_id);
    }
    int cursor = stable_id;
    for (int depth = 0; depth <= nedge; ++depth) {
      int parent = parent_stable_ids[cursor - 1];
      if (parent == 0) {
        break;
      }
      if (parent < 1 || parent > nedge || depth == nedge) {
        Rcpp::stop("V34 direct scenario registry contains a parent cycle");
      }
      cursor = parent;
    }
  }
  if (root_record_count < 1) {
    Rcpp::stop("V34 direct scenario registry contains no root record");
  }

  std::vector<int> leaf_stable_ids;
  std::vector<int> internal_stable_ids;
  for (int stable_id = 1; stable_id <= nedge; ++stable_id) {
    if (referenced_as_parent.count(stable_id)) {
      internal_stable_ids.push_back(stable_id);
    } else {
      leaf_stable_ids.push_back(stable_id);
    }
  }
  if (leaf_stable_ids.empty()) {
    Rcpp::stop("V34 direct scenario registry contains no terminal topology edge");
  }

  int ntips = static_cast<int>(leaf_stable_ids.size());
  std::vector<int> child_node_by_stable_id(
    static_cast<std::size_t>(nedge + 1),
    NA_INTEGER
  );
  for (int leaf_index = 0; leaf_index < ntips; ++leaf_index) {
    child_node_by_stable_id[leaf_stable_ids[leaf_index]] = leaf_index + 1;
  }
  int next_internal_node = ntips + 2;
  for (int stable_id : internal_stable_ids) {
    child_node_by_stable_id[stable_id] = next_internal_node++;
  }
  int root_node = ntips + 1;

  Rcpp::DataFrame state_lookup =
    Rcpp::as<Rcpp::DataFrame>(summary["state_lookup"]);
  Rcpp::IntegerVector lookup_state_ids = state_lookup["state_id"];
  Rcpp::CharacterVector lookup_states = state_lookup["state"];
  std::map<int, std::string> state_labels;
  for (int row = 0; row < state_lookup.nrows(); ++row) {
    state_labels[lookup_state_ids[row]] =
      Rcpp::as<std::string>(lookup_states[row]);
  }
  Rcpp::DataFrame path_lookup =
    Rcpp::as<Rcpp::DataFrame>(summary["path_lookup"]);
  Rcpp::IntegerVector lookup_path_ids = path_lookup["path_id"];
  Rcpp::CharacterVector lookup_paths = path_lookup["label"];
  std::map<int, std::string> path_labels;
  for (int row = 0; row < path_lookup.nrows(); ++row) {
    path_labels[lookup_path_ids[row]] =
      Rcpp::as<std::string>(lookup_paths[row]);
  }

  Rcpp::IntegerMatrix edge(nedge, 2);
  Rcpp::NumericVector edge_lengths(nedge);
  Rcpp::List maps(nedge);
  Rcpp::List path_maps(nedge);
  Rcpp::IntegerVector stable_ids(nedge);
  for (int edge_index = 0; edge_index < nedge; ++edge_index) {
    int stable_id = edge_index + 1;
    int parent_id = parent_stable_ids[edge_index];
    edge(edge_index, 0) =
      parent_id == 0 ? root_node : child_node_by_stable_id[parent_id];
    edge(edge_index, 1) = child_node_by_stable_id[stable_id];

    Rcpp::IntegerVector states = edge_state_ids[edge_index];
    Rcpp::IntegerVector paths = edge_path_ids[edge_index];
    Rcpp::IntegerVector sizes = edge_sizes[edge_index];
    Rcpp::NumericVector durations = edge_durations[edge_index];
    std::pair<V34TransitionMap, V34TransitionMap> serialized =
      pst_v34_make_scenario_maps(
        states,
        paths,
        sizes,
        durations,
        state_labels,
        path_labels
      );
    maps[edge_index] =
      pst_v34_scenario_edge_map_to_vector(serialized.first);
    path_maps[edge_index] =
      pst_v34_transition_map_to_vector(serialized.second);
    edge_lengths[edge_index] = pst_v34_map_sum(serialized.first);
    stable_ids[edge_index] = stable_id;
  }
  edge.attr("dimnames") = Rcpp::List::create(
    R_NilValue,
    Rcpp::CharacterVector::create("parent", "child")
  );

  Rcpp::CharacterVector input_tip_labels = summary["tip_labels"];
  Rcpp::CharacterVector scenario_tip_labels(ntips);
  for (int leaf_index = 0; leaf_index < ntips; ++leaf_index) {
    int stable_id = leaf_stable_ids[leaf_index];
    Rcpp::IntegerVector represented_tip_ids = edge_tip_ids[stable_id - 1];
    if (represented_tip_ids.size() < 1) {
      Rcpp::stop(
        "V34 direct scenario leaf lacks traversal-owned tip membership"
      );
    }
    std::ostringstream label;
    for (int member_index = 0;
         member_index < represented_tip_ids.size();
         ++member_index) {
      int tip_id = represented_tip_ids[member_index];
      if (tip_id < 1 || tip_id > input_tip_labels.size()) {
        Rcpp::stop("V34 direct scenario leaf contains an invalid tip ID");
      }
      if (member_index > 0) {
        label << ",";
      }
      label << Rcpp::as<std::string>(input_tip_labels[tip_id - 1]);
    }
    scenario_tip_labels[leaf_index] = label.str();
  }

  Rcpp::NumericMatrix mapped_edge(nedge, 0);
  Rcpp::List tree = Rcpp::List::create(
    Rcpp::Named("edge") = edge,
    Rcpp::Named("edge.length") = edge_lengths,
    Rcpp::Named("tip.label") = scenario_tip_labels,
    Rcpp::Named("Nnode") = nedge + 1 - ntips,
    Rcpp::Named("maps") = maps,
    Rcpp::Named("mapped.edge") = mapped_edge,
    Rcpp::Named("path_maps") = path_maps,
    Rcpp::Named("raw_scenario_edge_ids") = stable_ids
  );
  tree.attr("class") = Rcpp::CharacterVector::create("simmap", "phylo");
  return tree;
}

/**
 * Retained legacy formatter used only by older internal tree materializers.
 *
 * This is an Rcpp boundary: inputs are R objects and outputs must preserve documented V34 public shapes.
 */
static Rcpp::List pst_v34_legacy_canonicalize_scenario_tree_impl(
    Rcpp::List tree,
    Rcpp::List summary,
    double time_tolerance) {
  if (!std::isfinite(time_tolerance) || time_tolerance <= 0.0) {
    Rcpp::stop("V34 scenario-tree time tolerance must be finite and positive");
  }
  V34ScenarioTree scenario_tree = pst_v34_scenario_tree_from_r(
    tree,
    time_tolerance
  );
  pst_v34_legacy_prune_endpoint_leaf_cpp(
    scenario_tree,
    time_tolerance
  );
  pst_v34_contract_zero_length_edges_cpp(scenario_tree, time_tolerance);
  pst_v34_contract_unary_internal_edges_cpp(scenario_tree, time_tolerance);
  pst_v34_legacy_shift_common_prefix_cpp(scenario_tree, time_tolerance);
  pst_v34_legacy_prune_endpoint_leaf_cpp(
    scenario_tree,
    time_tolerance
  );
  pst_v34_contract_zero_length_edges_cpp(scenario_tree, time_tolerance);
  pst_v34_contract_unary_internal_edges_cpp(scenario_tree, time_tolerance);
  std::set<std::string> labels_before_terminal_merge(
    scenario_tree.tip_labels.begin(),
    scenario_tree.tip_labels.end()
  );
  Rcpp::CharacterVector original_tip_labels = summary["tip_labels"];
  Rcpp::IntegerVector trajectory_group_ids =
    summary["terminal_trajectory_group_id_by_tip"];
  // Terminal trajectory groups are traversal-owned and must remain aligned to
  // original biological tip labels before topology formatting consumes them.
  if (trajectory_group_ids.size() != original_tip_labels.size()) {
    Rcpp::stop("V34 terminal trajectory groups are not tip-aligned");
  }
  std::map<std::string, int> trajectory_group_by_tip;
  // Index each biological tip label by its committed state-duration group.
  for (int tip_id = 0; tip_id < original_tip_labels.size(); ++tip_id) {
    trajectory_group_by_tip[
      Rcpp::as<std::string>(original_tip_labels[tip_id])
    ] = trajectory_group_ids[tip_id];
  }
  int pre_merge_edge_count = static_cast<int>(scenario_tree.edges.size());
  pst_v34_legacy_collapse_peer_leaves_cpp(scenario_tree);
  pst_v34_legacy_collapse_endpoint_histories_cpp(
    scenario_tree,
    trajectory_group_by_tip,
    time_tolerance
  );
  // Tip or trajectory merging can expose unary internals; contract them before public serialization.
  if (static_cast<int>(scenario_tree.edges.size()) != pre_merge_edge_count) {
    pst_v34_contract_unary_internal_edges_cpp(scenario_tree, time_tolerance);
  }

  pst_v34_reorder_tip_labels_cpp(scenario_tree, original_tip_labels);
  std::vector<std::string> terminal_trajectory_merged_labels;
  // Record only composite labels newly created by whole-trajectory merging;
  // these labels require summary scenario-id synchronization below.
  for (const std::string& label : scenario_tree.tip_labels) {
    // A new label representing multiple tips identifies a merged terminal trajectory.
    if (labels_before_terminal_merge.count(label) == 0 &&
        pst_v34_split_csv(label).size() > 1) {
      terminal_trajectory_merged_labels.push_back(label);
    }
  }

  Rcpp::DataFrame state_lookup = Rcpp::as<Rcpp::DataFrame>(summary["state_lookup"]);
  Rcpp::IntegerVector state_id_column = state_lookup["state_id"];
  Rcpp::CharacterVector state_column = state_lookup["state"];
  std::map<int, std::string> state_labels;
  // Cache state labels for zero-duration terminal placeholders.
  for (int i = 0; i < state_id_column.size(); ++i) {
    state_labels[state_id_column[i]] = Rcpp::as<std::string>(state_column[i]);
  }
  Rcpp::IntegerVector observed_state_ids = summary["observed_state_ids"];
  std::string fallback_state = observed_state_ids.size() > 0 ?
    state_labels[observed_state_ids[0]] :
    std::string("NA");

  Rcpp::List final_frontier = summary["final_frontier"];
  Rcpp::IntegerVector final_state_ids = final_frontier["state_id"];
  Rcpp::IntegerVector final_path_ids = final_frontier["path_id"];
  Rcpp::DataFrame path_lookup = Rcpp::as<Rcpp::DataFrame>(summary["path_lookup"]);
  Rcpp::IntegerVector lookup_path_ids = path_lookup["path_id"];
  Rcpp::CharacterVector lookup_path_labels = path_lookup["label"];
  std::map<int, std::string> path_labels;
  for (int path_id = 0; path_id < lookup_path_ids.size(); ++path_id) {
    path_labels[lookup_path_ids[path_id]] =
      Rcpp::as<std::string>(lookup_path_labels[path_id]);
  }
  std::string fallback_path = path_labels.empty() ?
    std::string("NA") : path_labels.begin()->second;
  std::map<std::string, int> final_state_by_tip;
  std::map<std::string, int> final_path_by_tip;
  // Key each input tip's live final state by its stable biological label. This
  // frontier was captured at traversal commit and requires no matrix replay.
  for (int tip = 0; tip < original_tip_labels.size(); ++tip) {
    // Final frontier vectors must remain tip-aligned with the input labels.
    if (tip >= final_state_ids.size() || tip >= final_path_ids.size()) {
      Rcpp::stop("V34 final frontier is not aligned to tip labels");
    }
    std::string tip_label = Rcpp::as<std::string>(original_tip_labels[tip]);
    final_state_by_tip[tip_label] = final_state_ids[tip];
    final_path_by_tip[tip_label] = final_path_ids[tip];
  }

  int ntips = static_cast<int>(scenario_tree.tip_labels.size());
  // Give any empty terminal placeholder a zero-duration state/size label so
  // every public edge remains a valid named simmap vector.
  for (V34ScenarioEdge& edge : scenario_tree.edges) {
    // Existing mapped exposure already carries its own scenario identity.
    if (!edge.map.labels.empty()) {
      continue;
    }
    std::string state = fallback_state;
    std::string path = fallback_path;
    // Terminal placeholders can recover their state from the first represented biological tip.
    if (edge.child >= 1 && edge.child <= ntips) {
      std::vector<std::string> tips = pst_v34_split_csv(scenario_tree.tip_labels[edge.child - 1]);
      // An empty represented-tip label leaves the observed-state fallback unchanged.
      if (!tips.empty()) {
        auto state_it = final_state_by_tip.find(tips[0]);
        // Use the tip's committed final state when it exists in the summary matrix.
        if (state_it != final_state_by_tip.end()) {
          auto label_it = state_labels.find(state_it->second);
          // Replace the fallback only with a known state label.
          if (label_it != state_labels.end()) {
            state = label_it->second;
          }
        }
        std::map<std::string, int>::const_iterator path_it =
          final_path_by_tip.find(tips[0]);
        if (path_it != final_path_by_tip.end() &&
            path_labels.count(path_it->second)) {
          path = path_labels[path_it->second];
        }
      }
    }
    edge.map.labels.push_back(pst_v34_scenario_label(state, 1));
    edge.map.durations.push_back(0.0);
    edge.path_map.labels.push_back(path);
    edge.path_map.durations.push_back(0.0);
  }

  Rcpp::List out = pst_v34_scenario_tree_to_r(scenario_tree);
  out["terminal_events"] = pst_v34_terminal_events_from_tree(
    scenario_tree,
    summary,
    time_tolerance
  );
  Rcpp::CharacterVector merged_labels(terminal_trajectory_merged_labels.size());
  // Serialize merged-label diagnostics without changing public tree fields.
  for (int i = 0; i < static_cast<int>(terminal_trajectory_merged_labels.size()); ++i) {
    SET_STRING_ELT(merged_labels, i, Rf_mkCharCE(terminal_trajectory_merged_labels[i].c_str(), CE_UTF8));
  }
  out.attr("pst_v34_terminal_trajectory_merged_labels") = merged_labels;
  return out;
}

/**
 * Finalize direct scenario output by validating stable topology and exposing
 * traversal-owned state, path, support, and edge-type records unchanged.
 */
Rcpp::List pst_v34_canonicalize_scenario_tree_impl(
    Rcpp::List tree,
    Rcpp::List summary,
    double time_tolerance) {
  if (!std::isfinite(time_tolerance) || time_tolerance <= 0.0) {
    Rcpp::stop("V34 scenario-tree time tolerance must be finite and positive");
  }
  if (!tree.containsElementNamed("edge") ||
      !tree.containsElementNamed("edge.length") ||
      !tree.containsElementNamed("maps") ||
      !tree.containsElementNamed("path_maps") ||
      !tree.containsElementNamed("raw_scenario_edge_ids") ||
      !summary.containsElementNamed("scenario_edge_records")) {
    Rcpp::stop("V34 direct scenario finalization lacks required fields");
  }

  Rcpp::List records = summary["scenario_edge_records"];
  const char* required_fields[] = {
    "edge_durations", "edge_end_times", "edge_state_ids", "edge_path_ids",
    "edge_sizes", "edge_parent_scenario_edge_id", "edge_formation_time",
    "edge_type"
  };
  for (const char* field : required_fields) {
    if (!records.containsElementNamed(field)) {
      Rcpp::stop(
        std::string("V34 direct scenario finalization lacks ") + field
      );
    }
  }

  Rcpp::IntegerMatrix edge = tree["edge"];
  Rcpp::NumericVector edge_lengths = tree["edge.length"];
  Rcpp::List encoded_maps = tree["maps"];
  Rcpp::List path_maps = tree["path_maps"];
  Rcpp::IntegerVector stable_ids = tree["raw_scenario_edge_ids"];
  Rcpp::List durations_by_edge = records["edge_durations"];
  Rcpp::List end_times_by_edge = records["edge_end_times"];
  Rcpp::List state_ids_by_edge = records["edge_state_ids"];
  Rcpp::List path_ids_by_edge = records["edge_path_ids"];
  Rcpp::List sizes_by_edge = records["edge_sizes"];
  Rcpp::IntegerVector parent_stable_ids =
    records["edge_parent_scenario_edge_id"];
  Rcpp::NumericVector formation_times = records["edge_formation_time"];
  Rcpp::CharacterVector edge_types = records["edge_type"];
  int nedge = durations_by_edge.size();
  if (edge.ncol() != 2 ||
      edge.nrow() != nedge ||
      edge_lengths.size() != nedge ||
      encoded_maps.size() != nedge ||
      path_maps.size() != nedge ||
      stable_ids.size() != nedge ||
      end_times_by_edge.size() != nedge ||
      state_ids_by_edge.size() != nedge ||
      path_ids_by_edge.size() != nedge ||
      sizes_by_edge.size() != nedge ||
      parent_stable_ids.size() != nedge ||
      formation_times.size() != nedge ||
      edge_types.size() != nedge) {
    Rcpp::stop("V34 direct scenario output is not edge-parallel");
  }

  std::set<int> seen_stable_ids;
  std::map<int, int> public_edge_by_stable_id;
  std::map<int, int> public_edge_by_child_node;
  for (int edge_index = 0; edge_index < nedge; ++edge_index) {
    int stable_id = stable_ids[edge_index];
    if (stable_id < 1 || stable_id > nedge ||
        !seen_stable_ids.insert(stable_id).second) {
      Rcpp::stop("V34 direct scenario output has invalid stable IDs");
    }
    public_edge_by_stable_id[stable_id] = edge_index;
    if (!public_edge_by_child_node
          .insert(std::make_pair(edge(edge_index, 1), edge_index))
          .second) {
      Rcpp::stop("V34 direct scenario topology has duplicate child nodes");
    }
  }

  Rcpp::List public_maps(nedge);
  Rcpp::List lineage_size_maps(nedge);
  Rcpp::List singleton_provenance(nedge);
  for (int edge_index = 0; edge_index < nedge; ++edge_index) {
    int stable_id = stable_ids[edge_index];
    int record_index = stable_id - 1;
    int parent_stable_id = parent_stable_ids[record_index];
    std::map<int, int>::const_iterator topology_parent =
      public_edge_by_child_node.find(edge(edge_index, 0));
    int topology_parent_stable_id = topology_parent ==
      public_edge_by_child_node.end() ?
        0 : stable_ids[topology_parent->second];
    if (parent_stable_id != topology_parent_stable_id) {
      std::ostringstream diagnostic;
      diagnostic
        << "V34 direct scenario parent mismatch [stable_id="
        << stable_id
        << ", record_parent=" << parent_stable_id
        << ", topology_parent=" << topology_parent_stable_id
        << "]";
      Rcpp::stop(diagnostic.str());
    }
    std::string edge_type = Rcpp::as<std::string>(edge_types[record_index]);
    if ((parent_stable_id == 0 && edge_type != "Root") ||
        (parent_stable_id != 0 &&
         edge_type != "Stay" &&
         edge_type != "Leave")) {
      Rcpp::stop("V34 direct scenario record has an invalid edge type");
    }

    Rcpp::NumericVector durations = durations_by_edge[record_index];
    Rcpp::NumericVector end_times = end_times_by_edge[record_index];
    Rcpp::IntegerVector state_ids = state_ids_by_edge[record_index];
    Rcpp::IntegerVector path_ids = path_ids_by_edge[record_index];
    Rcpp::IntegerVector sizes = sizes_by_edge[record_index];
    Rcpp::NumericVector encoded_map = encoded_maps[edge_index];
    Rcpp::NumericVector path_map = path_maps[edge_index];
    Rcpp::CharacterVector encoded_names = encoded_map.names();
    Rcpp::CharacterVector path_names = path_map.names();
    int nsegment = durations.size();
    if (nsegment < 1 ||
        end_times.size() != nsegment ||
        state_ids.size() != nsegment ||
        path_ids.size() != nsegment ||
        sizes.size() != nsegment ||
        encoded_map.size() != nsegment ||
        path_map.size() != nsegment ||
        encoded_names.size() != nsegment ||
        path_names.size() != nsegment ||
        !R_finite(formation_times[record_index])) {
      Rcpp::stop("V34 direct scenario record has misaligned segments");
    }

    Rcpp::NumericVector state_map(nsegment);
    Rcpp::CharacterVector state_names(nsegment);
    Rcpp::NumericVector size_map(nsegment);
    double cursor = formation_times[record_index];
    for (int segment_index = 0;
         segment_index < nsegment;
         ++segment_index) {
      double duration = durations[segment_index];
      if (!R_finite(duration) || duration < 0.0 ||
          !R_finite(end_times[segment_index])) {
        Rcpp::stop("V34 direct scenario record has an invalid segment time");
      }
      cursor += duration;
      if (std::abs(cursor - end_times[segment_index]) > time_tolerance) {
        std::ostringstream diagnostic;
        diagnostic
          << "V34 direct scenario duration/end-time mismatch [stable_id="
          << stable_id
          << ", segment=" << (segment_index + 1)
          << "]";
        Rcpp::stop(diagnostic.str());
      }
      if (std::abs(encoded_map[segment_index] - duration) >
            time_tolerance ||
          std::abs(path_map[segment_index] - duration) >
            time_tolerance) {
        Rcpp::stop(
          "V34 direct scenario public map changed a committed duration"
        );
      }
      state_map[segment_index] = duration;
      state_names[segment_index] = pst_v34_scenario_state(
        Rcpp::as<std::string>(encoded_names[segment_index])
      );
      size_map[segment_index] = sizes[segment_index];
    }
    state_map.attr("names") = state_names;
    size_map.attr("names") = state_names;
    public_maps[edge_index] = state_map;
    lineage_size_maps[edge_index] = size_map;
    singleton_provenance[edge_index] =
      Rcpp::IntegerVector::create(stable_id);
    if (std::abs(edge_lengths[edge_index] -
          std::accumulate(
            durations.begin(),
            durations.end(),
            0.0
          )) > time_tolerance) {
      Rcpp::stop("V34 direct scenario edge length changed committed duration");
    }
  }

  Rcpp::List out = Rcpp::clone(tree);
  out["maps"] = public_maps;
  out["edge_type"] = Rcpp::clone(edge_types);
  out["size_maps"] = Rcpp::List::create(
    Rcpp::Named("lineage_size_maps") = lineage_size_maps
  );
  out["raw_scenario_edge_ids"] = singleton_provenance;
  out.attr("pst_v34_raw_scenario_edge_ids") = singleton_provenance;
  V34ScenarioTree terminal_event_tree =
    pst_v34_scenario_tree_from_r(tree, time_tolerance);
  out["terminal_events"] = pst_v34_terminal_events_from_tree(
    terminal_event_tree,
    summary,
    time_tolerance
  );
  out.attr("pst_v34_terminal_trajectory_merged_labels") =
    Rcpp::CharacterVector(0);
  out.attr("class") = Rcpp::CharacterVector::create("simmap", "phylo");
  return out;
}

/**
 * Retained legacy LDIF formatter for the isolated pre-direct implementation.
 *
 * Traversal scores survivor differentiation before scenario-tree cleanup. Cleanup
 * can contract a provisional raw Leave into a final Stay, so retain exactly the
 * raw event that supports each final public Leave and rebuild the cumulative path
 * counter from those retained events. Final source/destination paths and node
 * time establish identity; raw-edge provenance resolves ties when cleanup has
 * preserved that correspondence. Dense scenario matrices are not replayed.
 */
static void pst_v34_legacy_rewrite_differentiation(
    Rcpp::List& summary,
    const Rcpp::List& scenario_tree) {
  if (!summary.containsElementNamed("online_tt_core") ||
      !summary.containsElementNamed("path_lookup") ||
      !scenario_tree.containsElementNamed("edge") ||
      !scenario_tree.containsElementNamed("edge.length") ||
      !scenario_tree.containsElementNamed("edge_type") ||
      !scenario_tree.containsElementNamed("path_maps") ||
      !scenario_tree.containsElementNamed("raw_scenario_edge_ids")) {
    return;
  }

  Rcpp::List tt_core = Rcpp::clone(Rcpp::as<Rcpp::List>(
    summary["online_tt_core"]
  ));
  if (!tt_core.containsElementNamed("lineage_differentiating_records") ||
      !tt_core.containsElementNamed("ldif_path")) {
    return;
  }

  Rcpp::DataFrame records = Rcpp::as<Rcpp::DataFrame>(
    tt_core["lineage_differentiating_records"]
  );
  const char* required_record_columns[] = {
    "time", "from_state_id", "to_state_id", "phylo_edge_id",
    "source_path_id", "path_id", "scenario_edge_before",
    "scenario_edge_after"
  };
  for (const char* column : required_record_columns) {
    if (!records.containsElementNamed(column)) {
      Rcpp::stop("V34 LDIF record table lacks a required column");
    }
  }

  Rcpp::DataFrame path_lookup = Rcpp::as<Rcpp::DataFrame>(
    summary["path_lookup"]
  );
  if (!path_lookup.containsElementNamed("path_id") ||
      !path_lookup.containsElementNamed("label") ||
      !path_lookup.containsElementNamed("terminal_state_id")) {
    Rcpp::stop("V34 path lookup lacks final-state metadata required by LDIF");
  }
  Rcpp::IntegerVector lookup_path_id = path_lookup["path_id"];
  Rcpp::CharacterVector lookup_label = path_lookup["label"];
  Rcpp::IntegerVector lookup_terminal_state = path_lookup["terminal_state_id"];
  std::map<std::string, int> path_id_by_label;
  std::map<int, int> terminal_state_by_path_id;
  for (int row = 0; row < path_lookup.nrows(); ++row) {
    std::string label = Rcpp::as<std::string>(lookup_label[row]);
    int path_id = lookup_path_id[row];
    int terminal_state_id = lookup_terminal_state[row];
    if (path_id < 1 || terminal_state_id < 1 || path_id_by_label.count(label)) {
      Rcpp::stop("V34 path lookup is invalid or non-unique during LDIF finalization");
    }
    path_id_by_label[label] = path_id;
    terminal_state_by_path_id[path_id] = terminal_state_id;
  }

  Rcpp::IntegerMatrix edge = scenario_tree["edge"];
  Rcpp::NumericVector edge_length = scenario_tree["edge.length"];
  Rcpp::CharacterVector edge_type = scenario_tree["edge_type"];
  Rcpp::List path_maps = scenario_tree["path_maps"];
  Rcpp::List raw_edge_ids = scenario_tree["raw_scenario_edge_ids"];
  if (edge.ncol() != 2 || edge_length.size() != edge.nrow() ||
      edge_type.size() != edge.nrow() ||
      path_maps.size() != edge.nrow() || raw_edge_ids.size() != edge.nrow()) {
    Rcpp::stop("V34 final scenario metadata is not edge-parallel");
  }
  std::map<int, int> edge_by_child;
  for (int edge_id = 0; edge_id < edge.nrow(); ++edge_id) {
    edge_by_child[edge(edge_id, 1)] = edge_id;
  }
  int max_node = 0;
  std::set<int> child_nodes;
  for (int edge_id = 0; edge_id < edge.nrow(); ++edge_id) {
    max_node = std::max(max_node, std::max(edge(edge_id, 0), edge(edge_id, 1)));
    child_nodes.insert(edge(edge_id, 1));
  }
  int root_node = -1;
  for (int edge_id = 0; edge_id < edge.nrow(); ++edge_id) {
    int parent_node = edge(edge_id, 0);
    if (!child_nodes.count(parent_node)) {
      if (root_node >= 0 && root_node != parent_node) {
        Rcpp::stop("V34 final scenario topology has multiple roots");
      }
      root_node = parent_node;
    }
  }
  if (root_node < 0) {
    Rcpp::stop("V34 final scenario topology lacks a root");
  }
  std::vector<double> node_time(max_node + 1, NA_REAL);
  node_time[root_node] = 0.0;
  std::vector<bool> resolved(edge.nrow(), false);
  for (int pass = 0; pass < edge.nrow(); ++pass) {
    bool progressed = false;
    for (int edge_id = 0; edge_id < edge.nrow(); ++edge_id) {
      if (resolved[edge_id] ||
          !R_finite(node_time[edge(edge_id, 0)])) {
        continue;
      }
      node_time[edge(edge_id, 1)] =
        node_time[edge(edge_id, 0)] + edge_length[edge_id];
      resolved[edge_id] = true;
      progressed = true;
    }
    if (!progressed) break;
  }
  if (std::find(resolved.begin(), resolved.end(), false) != resolved.end()) {
    Rcpp::stop("V34 final scenario node times could not be resolved");
  }

  Rcpp::NumericVector record_time = records["time"];
  Rcpp::IntegerVector record_phylo_edge = records["phylo_edge_id"];
  Rcpp::IntegerVector record_source_path = records["source_path_id"];
  Rcpp::IntegerVector record_path = records["path_id"];
  Rcpp::IntegerVector record_scenario_after = records["scenario_edge_after"];
  std::vector<bool> record_used(records.nrows(), false);

  struct V34FinalLdifRecord {
    int raw_record_id;
    int public_edge_id;
    int public_parent_edge_id;
    int source_path_id;
    int path_id;
  };
  std::vector<V34FinalLdifRecord> final_records;

  for (int edge_id = 0; edge_id < edge.nrow(); ++edge_id) {
    if (Rcpp::as<std::string>(edge_type[edge_id]) != "Leave") {
      continue;
    }
    std::map<int, int>::const_iterator parent_it =
      edge_by_child.find(edge(edge_id, 0));
    if (parent_it == edge_by_child.end()) {
      Rcpp::stop("V34 final Leave lacks a parent scenario edge");
    }
    int parent_edge_id = parent_it->second;
    Rcpp::NumericVector parent_path_map = path_maps[parent_edge_id];
    Rcpp::NumericVector child_path_map = path_maps[edge_id];
    Rcpp::CharacterVector parent_names = parent_path_map.names();
    Rcpp::CharacterVector child_names = child_path_map.names();
    if (parent_names.size() < 1 || child_names.size() < 1) {
      Rcpp::stop("V34 final Leave lacks source or destination path identity");
    }
    std::string source_label = Rcpp::as<std::string>(
      parent_names[parent_names.size() - 1]
    );
    std::string destination_label = Rcpp::as<std::string>(child_names[0]);
    if (!path_id_by_label.count(source_label) ||
        !path_id_by_label.count(destination_label)) {
      Rcpp::stop("V34 final Leave references an unknown canonical path");
    }
    int source_path_id = path_id_by_label[source_label];
    int destination_path_id = path_id_by_label[destination_label];
    Rcpp::IntegerVector provenance = raw_edge_ids[edge_id];
    double public_event_time = node_time[edge(edge_id, 0)];

    int selected_record = -1;
    int selected_identity_score = -1;
    double selected_time_distance = R_PosInf;
    for (int record_id = 0; record_id < records.nrows(); ++record_id) {
      if (record_used[record_id]) {
        continue;
      }
      double time_distance = std::abs(record_time[record_id] - public_event_time);
      bool exact_source = record_source_path[record_id] == source_path_id;
      bool exact_destination = record_path[record_id] == destination_path_id;
      bool exact_pair = exact_source && exact_destination;
      // Cleanup may contract the provisional path boundary while retaining the
      // same differentiation node. Admit only an exact path pair or an event at
      // the exact final node time; unrelated records remain ineligible.
      if (!exact_pair && time_distance > 1e-8) {
        continue;
      }
      bool has_provenance = false;
      for (int raw_id : provenance) {
        if (record_scenario_after[record_id] == raw_id) {
          has_provenance = true;
          break;
        }
      }
      int identity_score =
        (exact_source ? 2 : 0) +
        (exact_destination ? 4 : 0) +
        (has_provenance ? 1 : 0);
      if (selected_record < 0 ||
          identity_score > selected_identity_score ||
          (identity_score == selected_identity_score &&
           time_distance < selected_time_distance)) {
        selected_record = record_id;
        selected_identity_score = identity_score;
        selected_time_distance = time_distance;
      }
    }
    if (selected_record < 0) {
      Rcpp::stop("V34 final Leave has no survivor-scored raw LDIF provenance");
    }
    record_used[selected_record] = true;
    final_records.push_back(V34FinalLdifRecord{
      selected_record,
      edge_id,
      parent_edge_id,
      source_path_id,
      destination_path_id
    });
  }

  std::sort(
    final_records.begin(),
    final_records.end(),
    [](const V34FinalLdifRecord& lhs, const V34FinalLdifRecord& rhs) {
      return lhs.raw_record_id < rhs.raw_record_id;
    }
  );

  int final_count = static_cast<int>(final_records.size());
  Rcpp::NumericVector final_time(final_count);
  Rcpp::IntegerVector final_from_state(final_count);
  Rcpp::IntegerVector final_to_state(final_count);
  Rcpp::IntegerVector final_phylo_edge(final_count);
  Rcpp::IntegerVector final_source_path(final_count);
  Rcpp::IntegerVector final_path(final_count);
  Rcpp::IntegerVector final_scenario_before(final_count);
  Rcpp::IntegerVector final_scenario_after(final_count);
  for (int row = 0; row < final_count; ++row) {
    const V34FinalLdifRecord& record = final_records[row];
    final_time[row] = record_time[record.raw_record_id];
    final_from_state[row] = terminal_state_by_path_id[record.source_path_id];
    final_to_state[row] = terminal_state_by_path_id[record.path_id];
    final_phylo_edge[row] = record_phylo_edge[record.raw_record_id];
    final_source_path[row] = record.source_path_id;
    final_path[row] = record.path_id;
    final_scenario_before[row] = record.public_parent_edge_id + 1;
    final_scenario_after[row] = record.public_edge_id + 1;
  }
  Rcpp::DataFrame reconciled_records = Rcpp::DataFrame::create(
    Rcpp::Named("time") = final_time,
    Rcpp::Named("from_state_id") = final_from_state,
    Rcpp::Named("to_state_id") = final_to_state,
    Rcpp::Named("phylo_edge_id") = final_phylo_edge,
    Rcpp::Named("source_path_id") = final_source_path,
    Rcpp::Named("path_id") = final_path,
    Rcpp::Named("scenario_edge_before") = final_scenario_before,
    Rcpp::Named("scenario_edge_after") = final_scenario_after,
    Rcpp::Named("stringsAsFactors") = false
  );
  reconciled_records.attr("event_levels") = records.attr("event_levels");

  Rcpp::NumericMatrix ldif_path = Rcpp::clone(Rcpp::as<Rcpp::NumericMatrix>(
    tt_core["ldif_path"]
  ));
  Rcpp::CharacterVector matrix_names = Rcpp::colnames(ldif_path);
  std::map<std::string, int> matrix_column_by_name;
  int total_column = -1;
  int time_column = -1;
  for (int column = 0; column < matrix_names.size(); ++column) {
    std::string name = Rcpp::as<std::string>(matrix_names[column]);
    matrix_column_by_name[name] = column;
    if (name == "total") total_column = column;
    if (name == "summary_time_vec") time_column = column;
  }
  if (total_column < 0 || time_column < 0) {
    Rcpp::stop("V34 LDIF path matrix lacks total or summary_time_vec");
  }
  for (int row = 0; row < ldif_path.nrow(); ++row) {
    for (int column = 0; column < ldif_path.ncol(); ++column) {
      if (column != time_column) {
        ldif_path(row, column) = 0.0;
      }
    }
  }
  for (const V34FinalLdifRecord& record : final_records) {
    std::string destination_label;
    for (const std::pair<const std::string, int>& item : path_id_by_label) {
      if (item.second == record.path_id) {
        destination_label = item.first;
        break;
      }
    }
    if (!matrix_column_by_name.count(destination_label)) {
      Rcpp::stop("V34 LDIF path matrix lacks a finalized destination path");
    }
    int destination_column = matrix_column_by_name[destination_label];
    double event_time = record_time[record.raw_record_id];
    for (int row = 0; row < ldif_path.nrow(); ++row) {
      if (ldif_path(row, time_column) + 1e-12 >= event_time) {
        ldif_path(row, destination_column) += 1.0;
        ldif_path(row, total_column) += 1.0;
      }
    }
  }

  tt_core["lineage_differentiating_records"] = reconciled_records;
  tt_core["ldif_path"] = ldif_path;
  summary["online_tt_core"] = tt_core;
}

/**
 * Rewrite summary scenario ids so matrices address the canonical public scenario-tree tips.
 *
 * This is an Rcpp boundary: inputs are R objects and outputs must preserve documented V34 public shapes.
 */
static Rcpp::List pst_v34_legacy_canonicalize_summary_scenario_ids_impl(
    Rcpp::List summary,
    Rcpp::List scenario_tree) {
  Rcpp::List out = Rcpp::clone(summary);
  pst_v34_legacy_rewrite_differentiation(out, scenario_tree);
  // Missing summary or label inputs mean no canonical scenario-id rewrite is possible.
  if (!out.containsElementNamed("scenario_summary") ||
      !out.containsElementNamed("tip_labels") ||
      !scenario_tree.containsElementNamed("tip.label")) {
    return out;
  }

  Rcpp::List scenario_summary = Rcpp::clone(Rcpp::as<Rcpp::List>(out["scenario_summary"]));
  // Scenario-id formatting and the traversal-owned final semantic prefix are
  // sufficient. State/path matrices are deliberately not consulted here.
  if (!scenario_summary.containsElementNamed("scenario_edge_ids") ||
      !out.containsElementNamed("online_tt_core")) {
    return out;
  }

  Rcpp::IntegerMatrix scenario_edge_ids = Rcpp::clone(Rcpp::as<Rcpp::IntegerMatrix>(
    scenario_summary["scenario_edge_ids"]
  ));
  Rcpp::List online_tt_core = out["online_tt_core"];
  // Missing compact semantic identities make equivalence unprovable; preserve
  // raw scenario ids rather than falling back to dense history replay.
  if (!online_tt_core.containsElementNamed("final_semantic_prefix_by_tip")) {
    return out;
  }
  Rcpp::IntegerVector final_semantic_prefix_by_tip =
    online_tt_core["final_semantic_prefix_by_tip"];
  Rcpp::CharacterVector original_tip_labels = out["tip_labels"];
  Rcpp::CharacterVector scenario_tip_labels = scenario_tree["tip.label"];
  Rcpp::CharacterVector merged_label_attr = scenario_tree.attr("pst_v34_terminal_trajectory_merged_labels");
  std::set<std::string> terminal_trajectory_merged_labels;
  // Load the exact composite labels created by terminal-trajectory merging.
  for (int i = 0; i < merged_label_attr.size(); ++i) {
    terminal_trajectory_merged_labels.insert(Rcpp::as<std::string>(merged_label_attr[i]));
  }
  // Without whole-trajectory merges, traversal scenario ids already match the public tree.
  if (terminal_trajectory_merged_labels.empty()) {
    return out;
  }

  std::map<std::string, int> row_by_tip;
  // Map biological tip labels back to rows in traversal-owned matrices.
  for (int row = 0; row < original_tip_labels.size(); ++row) {
    row_by_tip[Rcpp::as<std::string>(original_tip_labels[row])] = row;
  }

  // Inspect each public scenario tip and rewrite only explicitly merged terminal trajectories.
  for (int label_id = 0; label_id < scenario_tip_labels.size(); ++label_id) {
    std::string scenario_tip_label = Rcpp::as<std::string>(scenario_tip_labels[label_id]);
    // Ordinary tips and sibling-only merges retain their existing scenario ids.
    if (!terminal_trajectory_merged_labels.count(scenario_tip_label)) {
      continue;
    }
    std::vector<std::string> members = pst_v34_split_csv(
      scenario_tip_label
    );
    // Two-member labels do not require the broad terminal-history id rewrite used here.
    if (members.size() <= 2) {
      continue;
    }

    int primary_row = NA_INTEGER;
    // Select the first represented biological tip with a valid summary row as the canonical row.
    for (const std::string& member : members) {
      auto row_it = row_by_tip.find(member);
      // A known member supplies the candidate state/path history.
      if (row_it != row_by_tip.end()) {
        primary_row = row_it->second;
        break;
      }
    }
    // No known member means the merged label cannot be synchronized safely.
    if (primary_row == NA_INTEGER) {
      continue;
    }

    bool histories_match =
      primary_row >= 0 && primary_row < final_semantic_prefix_by_tip.size() &&
      final_semantic_prefix_by_tip[primary_row] > 0;
    int primary_prefix = histories_match ?
      final_semantic_prefix_by_tip[primary_row] : NA_INTEGER;
    // Verify every represented member shares the same traversal-owned recursive
    // semantic prefix. Equal final prefixes prove equal complete path histories.
    for (const std::string& member : members) {
      auto row_it = row_by_tip.find(member);
      // Unknown represented labels provide no tip identity and are ignored.
      if (row_it == row_by_tip.end()) {
        continue;
      }
      int member_row = row_it->second;
      // Missing or different compact prefixes prove equivalence was not established.
      if (member_row < 0 ||
          member_row >= final_semantic_prefix_by_tip.size() ||
          final_semantic_prefix_by_tip[member_row] != primary_prefix) {
        histories_match = false;
      }
      // Stop as soon as one represented member differs from the canonical history.
      if (!histories_match) {
        break;
      }
    }
    // Preserve original ids when semantic equivalence was not proven.
    if (!histories_match) {
      continue;
    }

    // Assign the canonical scenario-id trajectory to every equivalent represented tip.
    for (const std::string& member : members) {
      auto row_it = row_by_tip.find(member);
      // Skip labels absent from the original tip-row index.
      if (row_it == row_by_tip.end()) {
        continue;
      }
      int member_row = row_it->second;
      // Copy scenario ownership at every committed time while leaving state/path matrices unchanged.
      for (int col = 0; col < scenario_edge_ids.ncol(); ++col) {
        scenario_edge_ids(member_row, col) = scenario_edge_ids(primary_row, col);
      }
    }
  }

  scenario_summary["scenario_edge_ids"] = scenario_edge_ids;
  out["scenario_summary"] = scenario_summary;

  // Keep terminal-record scenario ids synchronized when terminal metadata is present.
  if (out.containsElementNamed("terminal_records")) {
    Rcpp::DataFrame terminal_records = Rcpp::clone(Rcpp::as<Rcpp::DataFrame>(
      out["terminal_records"]
    ));
    // Rewriting requires tip ids, final scenario ids, and at least one committed matrix column.
    if (terminal_records.containsElementNamed("tip_id") &&
        terminal_records.containsElementNamed("final_scenario_edge_id") &&
        scenario_edge_ids.ncol() > 0) {
      Rcpp::IntegerVector tip_id = terminal_records["tip_id"];
      Rcpp::IntegerVector final_scenario_edge_id = terminal_records["final_scenario_edge_id"];
      int final_col = scenario_edge_ids.ncol() - 1;
      // Replace each valid tip's terminal scenario id with its canonical final matrix value.
      for (int i = 0; i < terminal_records.nrows(); ++i) {
        int tip = tip_id[i];
        // Only in-domain biological tip ids may index the final scenario column.
        if (tip != NA_INTEGER && tip >= 1 && tip <= scenario_edge_ids.nrow()) {
          final_scenario_edge_id[i] = scenario_edge_ids(tip - 1, final_col);
        }
      }
      terminal_records["final_scenario_edge_id"] = final_scenario_edge_id;
      out["terminal_records"] = terminal_records;
    }
  }

  return out;
}

/**
 * Preserve the traversal-owned summary exactly. Stable scenario identities and
 * LDIF are committed together before tree serialization, so output formatting
 * has no summary rewrite authority.
 */
Rcpp::List pst_v34_canonicalize_summary_scenario_ids_impl(
    Rcpp::List summary,
    Rcpp::List scenario_tree) {
  if (!scenario_tree.containsElementNamed("raw_scenario_edge_ids")) {
    Rcpp::stop("V34 direct scenario tree lacks stable-ID provenance");
  }
  return Rcpp::clone(summary);
}

/**
 * Rebuild `scenario_mats$scenario_edge_ids` in the finalized scenario-tree edge
 * namespace without changing any other trajectory field.
 *
 * `get_trajectory_obj_v34()` calls this function only through the optional R
 * finalizer after the complete trajectory already exists. The finalized
 * `scenario_tree` supplies topology, edge lengths, and grouped tip labels;
 * `phylo_tip_labels` supplies the matrix-row order; `time_vec` supplies the
 * committed column times; and `existing_ids` supplies only the required matrix
 * dimensions and dimnames. The returned matrix contains rows of
 * `scenario_tree$edge`, uses the post-event child at a node-time boundary, and
 * carries a grouped scenario tip's terminal edge across every still-populated
 * input-matrix cell. Inputs are read only.
 */
Rcpp::IntegerMatrix pst_v34_final_scenario_edge_ids_impl(
    Rcpp::List scenario_tree,
    Rcpp::CharacterVector phylo_tip_labels,
    Rcpp::NumericVector time_vec,
    Rcpp::IntegerMatrix existing_ids,
    double time_tolerance) {
  // Validate and convert the shared tolerance through the same integer time-key
  // implementation used by V34 event batching. This prevents a second floating-
  // point interval convention from entering scenario-matrix reconstruction.
  double time_scale = pst_v34_batch_time_scale(time_tolerance);

  // The finalized scenario tree must expose one length for every two-column
  // edge row and one label for every scenario tip node.
  if (!scenario_tree.containsElementNamed("edge") ||
      !scenario_tree.containsElementNamed("edge.length") ||
      !scenario_tree.containsElementNamed("tip.label")) {
    Rcpp::stop("Final scenario tree lacks edge, edge.length, or tip.label");
  }
  Rcpp::IntegerMatrix edge = scenario_tree["edge"];
  Rcpp::NumericVector edge_length = scenario_tree["edge.length"];
  Rcpp::CharacterVector scenario_tip_labels = scenario_tree["tip.label"];
  int nedge = edge.nrow();
  if (edge.ncol() != 2 || nedge < 1 || edge_length.size() != nedge) {
    Rcpp::stop("Final scenario topology and edge lengths are not edge-parallel");
  }
  Rcpp::List raw_scenario_edge_ids = Rcpp::as<Rcpp::List>(
    scenario_tree.attr("pst_v34_raw_scenario_edge_ids")
  );
  if (raw_scenario_edge_ids.size() != nedge) {
    Rcpp::stop("Final scenario topology lacks edge-parallel raw provenance");
  }

  // Invert edge-parallel provenance once for terminal-only matrix rows. Without
  // this lookup, every populated fallback cell scans every final scenario edge
  // and searches its raw-id vector. Large empirical trees can have many
  // terminal-only tips and thousands of time columns, multiplying the same
  // provenance search into hundreds of millions of comparisons.
  std::map<int, std::vector<int> > final_edges_by_raw_id;
  // One pass represents one finalized scenario edge. Each raw traversal id in
  // that edge's provenance gains this edge as a possible public owner; later
  // cell resolution still applies the original interval and latest-edge rules.
  for (int edge_index = 0; edge_index < nedge; ++edge_index) {
    Rcpp::IntegerVector raw_ids = raw_scenario_edge_ids[edge_index];
    for (int raw_id : raw_ids) {
      if (raw_id != NA_INTEGER) {
        final_edges_by_raw_id[raw_id].push_back(edge_index);
      }
    }
  }

  // The matrix row and column domains must match the supplied phylogeny tips
  // and committed time vector before a replacement matrix is allocated.
  if (existing_ids.nrow() != phylo_tip_labels.size() ||
      existing_ids.ncol() != time_vec.size() ||
      existing_ids.nrow() < 1 || existing_ids.ncol() < 1) {
    Rcpp::stop("Scenario-edge id dimensions do not match phylogeny tips and times");
  }

  int max_node = 0;
  // Scan every final edge once to establish the node domain and reject invalid
  // nonpositive node ids before they are used as vector indices.
  for (int edge_index = 0; edge_index < nedge; ++edge_index) {
    int parent = edge(edge_index, 0);
    int child = edge(edge_index, 1);
    if (parent < 1 || child < 1) {
      Rcpp::stop("Final scenario topology contains a nonpositive node id");
    }
    max_node = std::max(max_node, std::max(parent, child));
  }

  std::vector<char> is_child(max_node + 1, 0);
  std::vector<int> edge_by_child(max_node + 1, -1);
  std::vector<std::vector<int> > outgoing(max_node + 1);
  // Build entering-edge and outgoing-edge indexes. One pass represents one
  // final edge and records the relationships needed for root discovery, node-
  // time propagation, and later tip-to-root walks.
  for (int edge_index = 0; edge_index < nedge; ++edge_index) {
    int parent = edge(edge_index, 0);
    int child = edge(edge_index, 1);
    if (edge_by_child[child] != -1) {
      Rcpp::stop("Final scenario topology gives one node multiple parents");
    }
    is_child[child] = 1;
    edge_by_child[child] = edge_index;
    outgoing[parent].push_back(edge_index);
  }

  int root = -1;
  // Inspect nodes that actually occur as parents. Exactly one such node may
  // lack an entering edge; that node is the scenario-tree root.
  for (int node = 1; node <= max_node; ++node) {
    if (!outgoing[node].empty() && !is_child[node]) {
      if (root != -1) {
        Rcpp::stop("Final scenario topology has multiple roots");
      }
      root = node;
    }
  }
  if (root == -1) {
    Rcpp::stop("Final scenario topology has no root");
  }

  std::vector<double> node_time(
    max_node + 1,
    std::numeric_limits<double>::quiet_NaN()
  );
  std::vector<int> node_queue;
  node_queue.reserve(max_node);
  node_time[root] = 0.0;
  node_queue.push_back(root);
  int visited_edges = 0;

  // Walk parent nodes from the root toward the tips. Each pass resolves every
  // outgoing edge's child time and queues that child for its own descendants.
  // The queue ends when every reachable scenario node has been processed.
  for (std::size_t queue_index = 0;
       queue_index < node_queue.size();
       ++queue_index) {
    int parent = node_queue[queue_index];
    // One inner pass resolves one edge leaving the current parent node.
    for (int edge_index : outgoing[parent]) {
      double length = edge_length[edge_index];
      if (!std::isfinite(length) || length < 0.0) {
        Rcpp::stop("Final scenario tree contains an invalid edge length");
      }
      int child = edge(edge_index, 1);
      double child_time = node_time[parent] + length;
      if (std::isfinite(node_time[child])) {
        Rcpp::stop("Final scenario topology reaches one child more than once");
      }
      node_time[child] = child_time;
      node_queue.push_back(child);
      ++visited_edges;
    }
  }
  if (visited_edges != nedge) {
    Rcpp::stop("Final scenario topology contains unreachable edges");
  }

  std::vector<long long> edge_start_key(nedge);
  std::vector<long long> edge_end_key(nedge);
  // Convert every finalized edge interval to the same integer time keys used by
  // V34 batching. One pass stores the inclusive start and end key for that edge.
  for (int edge_index = 0; edge_index < nedge; ++edge_index) {
    edge_start_key[edge_index] = pst_v34_batch_time_key(
      node_time[edge(edge_index, 0)],
      time_scale
    );
    edge_end_key[edge_index] = pst_v34_batch_time_key(
      node_time[edge(edge_index, 1)],
      time_scale
    );
    if (edge_end_key[edge_index] < edge_start_key[edge_index]) {
      Rcpp::stop("Final scenario edge has a negative integer-key interval");
    }
  }

  std::vector<long long> column_time_key(time_vec.size());
  // Quantize matrix columns once and verify their order. Each pass represents
  // one committed scenario-matrix column shared by every phylogeny-tip row.
  for (int column = 0; column < time_vec.size(); ++column) {
    column_time_key[column] = pst_v34_batch_time_key(
      time_vec[column],
      time_scale
    );
    if (column > 0 && column_time_key[column] < column_time_key[column - 1]) {
      Rcpp::stop("Scenario matrix time vector is not nondecreasing");
    }
  }

  std::map<std::string, int> scenario_tip_node_by_phylo_tip;
  // Decode each finalized scenario-tip label. One pass registers every
  // represented biological tip with the scenario tip node whose incoming edge
  // terminates its final root-to-tip chain.
  for (int tip_index = 0; tip_index < scenario_tip_labels.size(); ++tip_index) {
    std::string scenario_tip_label = Rcpp::as<std::string>(
      scenario_tip_labels[tip_index]
    );
    std::vector<std::string> members = pst_v34_split_csv(scenario_tip_label);
    if (members.empty()) {
      Rcpp::stop("Final scenario tree contains an empty tip label");
    }
    // One inner pass registers one biological tip represented by this scenario
    // tip and rejects duplicate membership across final scenario tips.
    for (const std::string& member : members) {
      if (scenario_tip_node_by_phylo_tip.count(member)) {
        Rcpp::stop("A phylogeny tip appears in multiple final scenario tips");
      }
      scenario_tip_node_by_phylo_tip[member] = tip_index + 1;
    }
  }

  Rcpp::IntegerMatrix final_ids(existing_ids.nrow(), existing_ids.ncol());
  std::fill(final_ids.begin(), final_ids.end(), NA_INTEGER);
  SEXP dimnames = existing_ids.attr("dimnames");
  if (dimnames != R_NilValue) {
    final_ids.attr("dimnames") = dimnames;
  }

  // Rebuild one matrix row per phylogeny tip. Each pass finds that tip's final
  // scenario tip, constructs its root-to-tip edge chain, and advances a single
  // chain pointer across all committed time columns.
  for (int row = 0; row < phylo_tip_labels.size(); ++row) {
    std::string phylo_tip_label = Rcpp::as<std::string>(phylo_tip_labels[row]);
    std::map<std::string, int>::const_iterator tip_it =
      scenario_tip_node_by_phylo_tip.find(phylo_tip_label);
    if (tip_it == scenario_tip_node_by_phylo_tip.end()) {
      // A terminal-only phylogenetic tip need not be a final scenario leaf.
      // Resolve its endpoint-inclusive historical cells through traversal-edge
      // provenance, carrying the last canonical edge after its endpoint.
      for (int column = 0; column < time_vec.size(); ++column) {
        int raw_edge_id = existing_ids(row, column);
        if (raw_edge_id == NA_INTEGER) {
          continue;
        }
        long long time_key = column_time_key[column];
        int canonical_edge_id = -1;
        long long latest_start_key = std::numeric_limits<long long>::min();
        int carried_edge_id = -1;
        long long latest_end_key = std::numeric_limits<long long>::min();

        // Inspect only final edges that actually descend from this cell's raw
        // traversal edge. One pass preserves the former global scan's choice:
        // the latest overlapping edge owns an in-range cell, while the edge
        // ending latest is carried after a terminal-only endpoint.
        std::map<int, std::vector<int> >::const_iterator provenance_it =
          final_edges_by_raw_id.find(raw_edge_id);
        if (provenance_it != final_edges_by_raw_id.end()) {
          for (int edge_index : provenance_it->second) {
            if (edge_end_key[edge_index] > latest_end_key) {
              carried_edge_id = edge_index;
              latest_end_key = edge_end_key[edge_index];
            }
            if (time_key >= edge_start_key[edge_index] &&
                time_key <= edge_end_key[edge_index] &&
                edge_start_key[edge_index] >= latest_start_key) {
              canonical_edge_id = edge_index;
              latest_start_key = edge_start_key[edge_index];
            }
          }
        }
        if (canonical_edge_id < 0) {
          canonical_edge_id = carried_edge_id;
        }
        if (canonical_edge_id < 0) {
          Rcpp::stop("Terminal-only matrix cell has no canonical scenario edge");
        }
        final_ids(row, column) = canonical_edge_id + 1;
      }
      continue;
    }

    std::vector<int> tip_to_root;
    int child = tip_it->second;
    // Follow entering edges until the root has no entering edge. One pass adds
    // one edge above the represented scenario tip; the node depth bound rejects
    // a malformed cycle instead of looping indefinitely.
    while (child >= 1 && child <= max_node && edge_by_child[child] != -1) {
      int edge_index = edge_by_child[child];
      tip_to_root.push_back(edge_index);
      child = edge(edge_index, 0);
      if (tip_to_root.size() > static_cast<std::size_t>(nedge)) {
        Rcpp::stop("Final scenario topology contains a cycle");
      }
    }
    if (tip_to_root.empty() || child != root) {
      Rcpp::stop("A final scenario tip does not trace back to the root");
    }
    std::reverse(tip_to_root.begin(), tip_to_root.end());

    std::size_t active_position = 0;
    long long terminal_key = edge_end_key[tip_to_root.back()];
    // Walk the shared time vector for this row. The active edge pointer moves
    // only toward the tip, making the row cost linear in columns plus depth.
    for (int column = 0; column < time_vec.size(); ++column) {
      long long time_key = column_time_key[column];

      // A finalized scenario tip can be shorter than one of the biological
      // tips grouped under the same complete history. In that case the raw
      // traversal matrix remains populated after the grouped scenario edge's
      // endpoint. Carry the final public edge across exactly those populated
      // cells so every existing state/phylogeny cell retains a scenario-tree
      // coordinate; truly missing raw cells remain missing.
      if (time_key > terminal_key) {
        if (existing_ids(row, column) != NA_INTEGER) {
          final_ids(row, column) = tip_to_root.back() + 1;
        }
        continue;
      }

      if (time_key < edge_start_key[tip_to_root.front()]) {
        Rcpp::stop("A scenario matrix column precedes its root scenario edge");
      }

      // At a node-time boundary, advance through every child starting on that
      // same integer key. This assigns the post-event, deepest child edge and
      // also handles zero-length chains deterministically.
      while (active_position + 1 < tip_to_root.size() &&
             edge_start_key[tip_to_root[active_position + 1]] <= time_key) {
        ++active_position;
      }

      int active_edge = tip_to_root[active_position];
      if (time_key > edge_end_key[active_edge]) {
        Rcpp::stop("No finalized scenario edge owns one active matrix column");
      }
      final_ids(row, column) = active_edge + 1;
    }
  }

  return final_ids;
}

/**
 * Precomputed time and label information for every mapped segment on one tree.
 *
 * Matrix construction repeatedly asks which segment of a known edge owns a
 * committed time column. Computing node depths and segment boundaries inside
 * that cell loop would multiply topology work by matrix size. This compact
 * layout pays for those calculations once per tree and lets each cell resolve
 * its step with a binary search over one edge's segment starts.
 */
struct V34MatrixMapLayout {
  std::vector<long long> edge_start_key;
  std::vector<long long> edge_end_key;
  std::vector<std::vector<long long> > segment_start_key;
  std::vector<std::vector<std::string> > state_label;
  std::vector<std::vector<std::string> > path_label;
};

/**
 * Convert one state-mapped tree into edge-local integer-key segment layouts.
 *
 * `include_paths` controls whether parallel path-map labels are read and
 * cached. State maps are always required because their durations define the
 * step boundaries. Topology is traversed root-first so every edge start time is
 * inherited from its parent node rather than inferred from dense matrices.
 */
static V34MatrixMapLayout pst_v34_matrix_map_layout(
    const Rcpp::List& tree,
    double time_scale,
    bool include_paths,
    const std::string& tree_name) {
  if (!tree.containsElementNamed("edge") ||
      !tree.containsElementNamed("edge.length") ||
      !tree.containsElementNamed("maps")) {
    Rcpp::stop(tree_name + " lacks edge, edge.length, or maps");
  }

  Rcpp::IntegerMatrix edge = tree["edge"];
  Rcpp::NumericVector edge_length = tree["edge.length"];
  Rcpp::List maps = tree["maps"];
  int nedge = edge.nrow();
  if (edge.ncol() != 2 || nedge < 1 ||
      edge_length.size() != nedge || maps.size() != nedge) {
    Rcpp::stop(tree_name + " topology, lengths, and maps are not edge-parallel");
  }

  Rcpp::List path_maps;
  if (include_paths) {
    if (!tree.containsElementNamed("path_maps")) {
      Rcpp::stop(tree_name + " lacks path_maps requested for matrix output");
    }
    path_maps = Rcpp::as<Rcpp::List>(tree["path_maps"]);
    if (path_maps.size() != nedge) {
      Rcpp::stop(tree_name + " path maps are not edge-parallel");
    }
  }

  int max_node = 0;
  // Establish the node domain and reject invalid ids before allocating node-
  // indexed topology vectors.
  for (int edge_id = 0; edge_id < nedge; ++edge_id) {
    int parent = edge(edge_id, 0);
    int child = edge(edge_id, 1);
    if (parent < 1 || child < 1) {
      Rcpp::stop(tree_name + " topology contains a nonpositive node id");
    }
    max_node = std::max(max_node, std::max(parent, child));
  }

  std::vector<char> is_child(max_node + 1, 0);
  std::vector<std::vector<int> > outgoing(max_node + 1);
  // Record every outgoing edge and ensure no child receives two parents.
  std::vector<int> entering_edge(max_node + 1, -1);
  for (int edge_id = 0; edge_id < nedge; ++edge_id) {
    int parent = edge(edge_id, 0);
    int child = edge(edge_id, 1);
    if (entering_edge[child] != -1) {
      Rcpp::stop(tree_name + " topology gives one node multiple parents");
    }
    entering_edge[child] = edge_id;
    is_child[child] = 1;
    outgoing[parent].push_back(edge_id);
  }

  int root = -1;
  // A valid rooted tree has exactly one parent node without an entering edge.
  for (int node = 1; node <= max_node; ++node) {
    if (!outgoing[node].empty() && !is_child[node]) {
      if (root != -1) {
        Rcpp::stop(tree_name + " topology has multiple roots");
      }
      root = node;
    }
  }
  if (root == -1) {
    Rcpp::stop(tree_name + " topology has no root");
  }

  std::vector<double> node_time(
    max_node + 1,
    std::numeric_limits<double>::quiet_NaN()
  );
  std::vector<int> queue;
  queue.reserve(max_node);
  node_time[root] = 0.0;
  queue.push_back(root);
  int visited_edges = 0;

  // Propagate absolute node times from the root through every edge once.
  for (std::size_t queue_id = 0; queue_id < queue.size(); ++queue_id) {
    int parent = queue[queue_id];
    for (int edge_id : outgoing[parent]) {
      double length = edge_length[edge_id];
      if (!std::isfinite(length) || length < 0.0) {
        Rcpp::stop(tree_name + " contains an invalid edge length");
      }
      int child = edge(edge_id, 1);
      if (std::isfinite(node_time[child])) {
        Rcpp::stop(tree_name + " topology reaches one child more than once");
      }
      node_time[child] = node_time[parent] + length;
      queue.push_back(child);
      ++visited_edges;
    }
  }
  if (visited_edges != nedge) {
    Rcpp::stop(tree_name + " topology contains unreachable edges");
  }

  V34MatrixMapLayout layout;
  layout.edge_start_key.resize(nedge);
  layout.edge_end_key.resize(nedge);
  layout.segment_start_key.resize(nedge);
  layout.state_label.resize(nedge);
  if (include_paths) {
    layout.path_label.resize(nedge);
  }

  // Build one segment-start vector per edge. Segment durations are accumulated
  // in double precision and every boundary is quantized with V34's single
  // construction tolerance, matching event batching and final edge ownership.
  for (int edge_id = 0; edge_id < nedge; ++edge_id) {
    double edge_start = node_time[edge(edge_id, 0)];
    double edge_end = node_time[edge(edge_id, 1)];
    layout.edge_start_key[edge_id] = pst_v34_batch_time_key(
      edge_start,
      time_scale
    );
    layout.edge_end_key[edge_id] = pst_v34_batch_time_key(
      edge_end,
      time_scale
    );

    Rcpp::NumericVector edge_map = maps[edge_id];
    Rcpp::CharacterVector state_names = edge_map.names();
    if (edge_map.size() < 1 || state_names.size() != edge_map.size()) {
      Rcpp::stop(tree_name + " contains an empty or unnamed state map");
    }
    layout.segment_start_key[edge_id].reserve(edge_map.size());
    layout.state_label[edge_id].reserve(edge_map.size());

    Rcpp::NumericVector edge_path_map;
    Rcpp::CharacterVector path_names;
    if (include_paths) {
      edge_path_map = Rcpp::as<Rcpp::NumericVector>(path_maps[edge_id]);
      path_names = edge_path_map.names();
      if (edge_path_map.size() != edge_map.size() ||
          path_names.size() != edge_map.size()) {
        Rcpp::stop(tree_name + " state and path map steps are not aligned");
      }
      layout.path_label[edge_id].reserve(edge_map.size());
    }

    double segment_start = edge_start;
    for (int step_id = 0; step_id < edge_map.size(); ++step_id) {
      double duration = edge_map[step_id];
      if (!std::isfinite(duration) || duration < 0.0) {
        Rcpp::stop(tree_name + " contains an invalid map duration");
      }
      layout.segment_start_key[edge_id].push_back(
        pst_v34_batch_time_key(segment_start, time_scale)
      );
      if (state_names[step_id] == NA_STRING) {
        Rcpp::stop(tree_name + " contains a missing state-map label");
      }
      layout.state_label[edge_id].push_back(
        Rcpp::as<std::string>(state_names[step_id])
      );
      if (include_paths) {
        if (path_names[step_id] == NA_STRING) {
          Rcpp::stop(tree_name + " contains a missing path-map label");
        }
        layout.path_label[edge_id].push_back(
          Rcpp::as<std::string>(path_names[step_id])
        );
      }
      segment_start += duration;
    }
  }

  return layout;
}

/**
 * Resolve one 1-based tree edge and integer time key to a 1-based map step.
 *
 * `upper_bound` chooses the latest segment beginning on or before the matrix
 * time. Therefore a state change exactly on a committed time key selects the
 * post-event step, including several zero-length/collapsed starts on one key.
 * A terminal carry cell after a shorter tip's endpoint resolves to that edge's
 * final map step.
 */
static int pst_v34_matrix_map_step(
    const V34MatrixMapLayout& layout,
    int edge_id,
    long long time_key,
    const std::string& tree_name) {
  int edge_index = edge_id - 1;
  if (edge_index < 0 ||
      edge_index >= static_cast<int>(layout.segment_start_key.size())) {
    Rcpp::stop(tree_name + " matrix contains an invalid edge id");
  }
  const std::vector<long long>& starts = layout.segment_start_key[edge_index];
  if (starts.empty()) {
    Rcpp::stop(tree_name + " edge has no mapped steps");
  }
  if (time_key < layout.edge_start_key[edge_index]) {
    Rcpp::stop(tree_name + " matrix time precedes its addressed edge");
  }
  if (time_key >= layout.edge_end_key[edge_index]) {
    return static_cast<int>(starts.size());
  }

  std::vector<long long>::const_iterator after = std::upper_bound(
    starts.begin(),
    starts.end(),
    time_key
  );
  if (after == starts.begin()) {
    Rcpp::stop(tree_name + " matrix time precedes its first mapped step");
  }
  return static_cast<int>(after - starts.begin());
}

/**
 * Resolve an endpoint-inclusive terminal cell across a post-event scenario step.
 *
 * Scenario continuation can change state/path in the same atomic batch that an
 * unrelated member terminates. The continued edge's ordinary time lookup then
 * selects the post-event survivor step, while the terminal row still owns the
 * immediately preceding endpoint step. Rewind only across segment starts on
 * this exact time key and only to an exact phylogeny state/path match. A cell
 * that is not a terminal endpoint, or whose scenario history has no such
 * boundary, retains the ordinary post-event step.
 */
static int pst_v34_terminal_endpoint_scenario_step(
    const V34MatrixMapLayout& scenario_layout,
    int scenario_edge_id,
    int ordinary_step_id,
    long long time_key,
    const std::string& endpoint_state,
    const std::string& endpoint_path,
    bool include_paths) {
  int edge_index = scenario_edge_id - 1;
  int ordinary_index = ordinary_step_id - 1;
  if (edge_index < 0 ||
      edge_index >= static_cast<int>(scenario_layout.segment_start_key.size()) ||
      ordinary_index < 0 ||
      ordinary_index >= static_cast<int>(
        scenario_layout.segment_start_key[edge_index].size()
      )) {
    Rcpp::stop("Scenario terminal endpoint lookup received an invalid step");
  }

  const std::vector<long long>& starts =
    scenario_layout.segment_start_key[edge_index];
  const std::vector<std::string>& states =
    scenario_layout.state_label[edge_index];
  for (int candidate = ordinary_index - 1; candidate >= 0; --candidate) {
    // The candidate must end at this exact atomic boundary. Once its following
    // step begins earlier, older history cannot describe the endpoint cell.
    if (starts[candidate + 1] < time_key) {
      break;
    }
    if (starts[candidate + 1] != time_key ||
        states[candidate] != endpoint_state) {
      continue;
    }
    if (include_paths &&
        scenario_layout.path_label[edge_index][candidate] != endpoint_path) {
      continue;
    }
    return candidate + 1;
  }
  return ordinary_step_id;
}

/**
 * Build exact scenario and phylogeny map coordinates for every matrix cell.
 *
 * The existing scenario labels are used only for dimensions/dimnames when the
 * finalized edge matrix is rebuilt. Existing phylogeny edge ids already refer
 * to the unchanged input-tree edge rows. Parallel step matrices then identify
 * the final map element on each addressed edge. Optional character matrices
 * expose the path label resolved by those coordinates without changing either
 * input tree or any unrelated trajectory field.
 */
Rcpp::List pst_v34_scenario_matrix_coordinates_impl(
    Rcpp::List scenario_tree,
    Rcpp::List phylo_tree,
    Rcpp::NumericVector time_vec,
    Rcpp::IntegerMatrix existing_scenario_ids,
    Rcpp::IntegerMatrix existing_phylo_ids,
    bool include_paths,
    double time_tolerance) {
  if (existing_scenario_ids.nrow() != existing_phylo_ids.nrow() ||
      existing_scenario_ids.ncol() != existing_phylo_ids.ncol() ||
      existing_scenario_ids.ncol() != time_vec.size()) {
    Rcpp::stop("Scenario and phylogeny edge-id matrices are not aligned");
  }
  if (!phylo_tree.containsElementNamed("tip.label")) {
    Rcpp::stop("Phylogeny lacks tip labels required for scenario-matrix rows");
  }
  Rcpp::CharacterVector phylo_tip_labels = phylo_tree["tip.label"];
  if (phylo_tip_labels.size() != existing_scenario_ids.nrow()) {
    Rcpp::stop("Phylogeny tips do not match scenario-matrix rows");
  }

  // Reuse the finalized topology traversal that already implements terminal
  // merged-tip membership and post-event edge ownership.
  Rcpp::IntegerMatrix scenario_edge_ids =
    pst_v34_final_scenario_edge_ids_impl(
      scenario_tree,
      phylo_tip_labels,
      time_vec,
      existing_scenario_ids,
      time_tolerance
    );

  double time_scale = pst_v34_batch_time_scale(time_tolerance);
  V34MatrixMapLayout scenario_layout = pst_v34_matrix_map_layout(
    scenario_tree,
    time_scale,
    include_paths,
    "Scenario tree"
  );
  V34MatrixMapLayout phylo_layout = pst_v34_matrix_map_layout(
    phylo_tree,
    time_scale,
    include_paths,
    "Phylogeny"
  );
  Rcpp::IntegerMatrix phylo_edge = phylo_tree["edge"];

  Rcpp::IntegerMatrix scenario_step_ids(
    existing_scenario_ids.nrow(),
    existing_scenario_ids.ncol()
  );
  Rcpp::IntegerMatrix phylo_step_ids(
    existing_phylo_ids.nrow(),
    existing_phylo_ids.ncol()
  );
  std::fill(scenario_step_ids.begin(), scenario_step_ids.end(), NA_INTEGER);
  std::fill(phylo_step_ids.begin(), phylo_step_ids.end(), NA_INTEGER);

  Rcpp::CharacterMatrix scenario_paths;
  Rcpp::CharacterMatrix phylo_paths;
  if (include_paths) {
    scenario_paths = Rcpp::CharacterMatrix(
      existing_scenario_ids.nrow(),
      existing_scenario_ids.ncol()
    );
    phylo_paths = Rcpp::CharacterMatrix(
      existing_phylo_ids.nrow(),
      existing_phylo_ids.ncol()
    );
    std::fill(scenario_paths.begin(), scenario_paths.end(), NA_STRING);
    std::fill(phylo_paths.begin(), phylo_paths.end(), NA_STRING);
  }

  SEXP dimnames = existing_scenario_ids.attr("dimnames");
  if (dimnames != R_NilValue) {
    scenario_step_ids.attr("dimnames") = dimnames;
    phylo_step_ids.attr("dimnames") = dimnames;
    if (include_paths) {
      scenario_paths.attr("dimnames") = dimnames;
      phylo_paths.attr("dimnames") = dimnames;
    }
  }

  std::vector<long long> time_key(time_vec.size());
  for (int column = 0; column < time_vec.size(); ++column) {
    time_key[column] = pst_v34_batch_time_key(time_vec[column], time_scale);
  }

  // Resolve both tree coordinates for one rectangular cell at a time. Missing
  // edge ids remain missing step/path cells; mismatched coverage is rejected so
  // a saved matrix never exposes only half of a biological correspondence.
  for (int column = 0; column < existing_scenario_ids.ncol(); ++column) {
    for (int row = 0; row < existing_scenario_ids.nrow(); ++row) {
      int scenario_edge_id = scenario_edge_ids(row, column);
      int phylo_edge_id = existing_phylo_ids(row, column);
      bool scenario_missing = scenario_edge_id == NA_INTEGER;
      bool phylo_missing = phylo_edge_id == NA_INTEGER;
      if (scenario_missing != phylo_missing) {
        Rcpp::stop(
          "Scenario and phylogeny edge-id coverage does not match at row %d, "
          "column %d, time %.17g, tip %s: scenario_edge_id=%d, "
          "raw_scenario_edge_id=%d, phylo_edge_id=%d",
          row + 1,
          column + 1,
          time_vec[column],
          Rcpp::as<std::string>(phylo_tip_labels[row]).c_str(),
          scenario_edge_id,
          existing_scenario_ids(row, column),
          phylo_edge_id
        );
      }
      if (scenario_missing) {
        continue;
      }

      int phylo_step_id = pst_v34_matrix_map_step(
        phylo_layout,
        phylo_edge_id,
        time_key[column],
        "Phylogeny"
      );
      int scenario_step_id = pst_v34_matrix_map_step(
        scenario_layout,
        scenario_edge_id,
        time_key[column],
        "Scenario tree"
      );
      int phylo_edge_index = phylo_edge_id - 1;
      bool terminal_endpoint =
        phylo_edge_index >= 0 && phylo_edge_index < phylo_edge.nrow() &&
        phylo_edge(phylo_edge_index, 1) == row + 1 &&
        time_key[column] == phylo_layout.edge_end_key[phylo_edge_index];
      if (terminal_endpoint) {
        const std::string& endpoint_state =
          phylo_layout.state_label[phylo_edge_index][phylo_step_id - 1];
        std::string endpoint_path;
        if (include_paths) {
          endpoint_path =
            phylo_layout.path_label[phylo_edge_index][phylo_step_id - 1];
        }
        scenario_step_id = pst_v34_terminal_endpoint_scenario_step(
          scenario_layout,
          scenario_edge_id,
          scenario_step_id,
          time_key[column],
          endpoint_state,
          endpoint_path,
          include_paths
        );
      }
      scenario_step_ids(row, column) = scenario_step_id;
      phylo_step_ids(row, column) = phylo_step_id;

      if (include_paths) {
        scenario_paths(row, column) =
          scenario_layout.path_label[scenario_edge_id - 1][scenario_step_id - 1];
        phylo_paths(row, column) =
          phylo_layout.path_label[phylo_edge_id - 1][phylo_step_id - 1];
      }
    }
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("scenario_edge_ids") = scenario_edge_ids,
    Rcpp::Named("scenario_edge_step_ids") = scenario_step_ids,
    Rcpp::Named("phylo_edge_step_ids") = phylo_step_ids
  );
  if (include_paths) {
    out["scenario_paths"] = scenario_paths;
    out["phylo_paths"] = phylo_paths;
  }
  return out;
}

/**
 * Serialize the public transition tree from traversal-owned canonical paths.
 *
 * Each input edge row is one path that appeared in a committed batch. Parent
 * identity, state, lineage support, represented tips, and terminal groups are
 * already fixed by traversal. This function only assigns public node numbers,
 * state-size labels, and simmap containers; it never groups a dense history or
 * infers a transition sequence from completed matrices.
 */
static Rcpp::List pst_v34_transition_tree_from_compact_support(
    const Rcpp::List& transition_support,
    const Rcpp::CharacterVector& tip_labels,
    double time_tolerance) {
  Rcpp::IntegerVector edge_parent_path_id =
    transition_support["edge_parent_path_id"];
  Rcpp::IntegerVector edge_path_id = transition_support["edge_path_id"];
  Rcpp::IntegerVector edge_state_id = transition_support["edge_state_id"];
  Rcpp::IntegerVector edge_lineage_size =
    transition_support["edge_lineage_size"];
  Rcpp::List edge_tip_ids = transition_support["edge_tip_ids"];
  Rcpp::IntegerVector terminal_path_id =
    transition_support["terminal_path_id"];
  Rcpp::List terminal_tip_ids = transition_support["terminal_tip_ids"];
  Rcpp::IntegerVector path_parent_by_id =
    transition_support["path_parent_by_id"];
  Rcpp::CharacterVector state_labels_by_id =
    transition_support["state_labels_by_id"];
  Rcpp::CharacterVector expanded_state_space =
    transition_support["expanded_state_space"];

  int nrecord = edge_path_id.size();
  // Compact edge fields are record-parallel. A mismatch would pair one path's
  // topology with another path's state, support, or represented tips.
  if (edge_parent_path_id.size() != nrecord ||
      edge_state_id.size() != nrecord ||
      edge_lineage_size.size() != nrecord ||
      edge_tip_ids.size() != nrecord) {
    Rcpp::stop("V34 compact transition support fields are not aligned");
  }
  // Terminal path ids and represented-tip groups are one-to-one.
  if (terminal_path_id.size() != terminal_tip_ids.size()) {
    Rcpp::stop("V34 compact transition terminal groups are not aligned");
  }

  std::map<int, std::string> state_labels;
  Rcpp::CharacterVector state_label_names = state_labels_by_id.names();
  // Decode the immutable state lookup once for state-size label formatting.
  for (int index = 0; index < state_labels_by_id.size(); ++index) {
    int state_id = std::atoi(
      Rcpp::as<std::string>(state_label_names[index]).c_str()
    );
    state_labels[state_id] = Rcpp::as<std::string>(
      state_labels_by_id[index]
    );
  }

  std::set<int> emitted_paths;
  std::map<int, int> state_by_path;
  std::map<int, int> size_by_path;
  std::map<int, std::vector<int> > support_tips_by_path;
  // Index each committed path record. The path id itself is the stable old
  // child-node id used until deterministic public renumbering.
  for (int record_id = 0; record_id < nrecord; ++record_id) {
    int path_id = edge_path_id[record_id];
    // Duplicate or out-of-domain path records would create ambiguous topology.
    if (path_id < 1 ||
        path_id > path_parent_by_id.size() ||
        emitted_paths.count(path_id)) {
      Rcpp::stop("V34 compact transition support contains an invalid path");
    }
    emitted_paths.insert(path_id);
    state_by_path[path_id] = edge_state_id[record_id];
    size_by_path[path_id] = edge_lineage_size[record_id];
    support_tips_by_path[path_id] = Rcpp::as<std::vector<int> >(
      edge_tip_ids[record_id]
    );
  }

  std::map<int, std::vector<int> > terminal_tips_by_path;
  // Index traversal-owned terminal history groups by final canonical path.
  for (int terminal_id = 0;
       terminal_id < terminal_path_id.size();
       ++terminal_id) {
    int path_id = terminal_path_id[terminal_id];
    // A terminal history must end on a path represented by a committed edge.
    if (!emitted_paths.count(path_id)) {
      Rcpp::stop("V34 terminal history references an uncommitted path");
    }
    terminal_tips_by_path[path_id] = Rcpp::as<std::vector<int> >(
      terminal_tip_ids[terminal_id]
    );
  }

  std::map<int, int> effective_parent_by_path;
  std::map<int, std::vector<int> > children_by_path;
  int max_old_id = path_parent_by_id.size();
  // Resolve every emitted path to its nearest emitted ancestor. Paths absent
  // from committed rows are zero-duration intermediates and are bypassed using
  // the complete live parent registry, never by parsing public labels.
  for (int path_id : emitted_paths) {
    int parent_path_id = path_parent_by_id[path_id - 1];
    // Climb only through registered, uncommitted ancestors. Each iteration
    // moves toward the root and cannot revisit a path in the canonical tree.
    while (parent_path_id != NA_INTEGER &&
           parent_path_id > 0 &&
           !emitted_paths.count(parent_path_id)) {
      // Parent ids outside the registry indicate corrupted traversal support.
      if (parent_path_id > path_parent_by_id.size()) {
        Rcpp::stop("V34 transition parent path is outside the registry");
      }
      parent_path_id = path_parent_by_id[parent_path_id - 1];
    }
    effective_parent_by_path[path_id] = parent_path_id;
    // Only emitted parents own public transition-tree children. Initial paths
    // remain attached to the synthetic public root established below.
    if (parent_path_id != NA_INTEGER &&
        parent_path_id > 0 &&
        emitted_paths.count(parent_path_id)) {
      children_by_path[parent_path_id].push_back(path_id);
    }
  }

  int root_old_id = ++max_old_id;
  std::vector<std::pair<int, int> > old_edges;
  std::vector<V34TransitionMap> ordered_maps;
  std::map<int, std::vector<int> > represented_tips_by_leaf;
  // Emit one unit transition edge per committed canonical path. Unit length is
  // transition depth, matching the established T-tree diversity definition.
  for (int path_id : emitted_paths) {
    int parent_path_id = effective_parent_by_path[path_id];
    int parent_old_id = root_old_id;
    // A represented emitted ancestor supplies the public parent node.
    if (parent_path_id != NA_INTEGER &&
        parent_path_id > 0 &&
        emitted_paths.count(parent_path_id)) {
      parent_old_id = parent_path_id;
    }
    int state_id = state_by_path[path_id];
    int lineage_size = size_by_path[path_id];
    // Every public transition segment requires a known state and positive live
    // support; malformed support stops before a partial simmap is returned.
    if (!state_labels.count(state_id) || lineage_size < 1) {
      Rcpp::stop("V34 compact transition edge has invalid state or support");
    }
    old_edges.push_back(std::make_pair(parent_old_id, path_id));
    ordered_maps.push_back(pst_v34_make_transition_map(
      std::vector<int>(1, state_id),
      std::vector<int>(1, lineage_size),
      std::vector<double>(1, 1.0),
      state_labels
    ));

    bool has_children = !children_by_path[path_id].empty();
    std::map<int, std::vector<int> >::const_iterator terminal =
      terminal_tips_by_path.find(path_id);
    bool has_terminal_tips =
      terminal != terminal_tips_by_path.end() && !terminal->second.empty();
    // A path with no descendants is itself a terminal transition history.
    // During partial debug recovery, active support tips become deterministic
    // provisional tips when no biological terminal has yet closed.
    if (!has_children) {
      represented_tips_by_leaf[path_id] = has_terminal_tips ?
        terminal->second : support_tips_by_path[path_id];
    } else if (has_terminal_tips) {
      // A history ending at an internal path needs a terminal placeholder so
      // its tips coexist with longer child histories. STOPPED is removed by
      // the existing pure formatter, leaving a zero-length inherited-state tip.
      int terminal_old_id = ++max_old_id;
      old_edges.push_back(std::make_pair(path_id, terminal_old_id));
      V34TransitionMap stopped;
      stopped.labels.push_back("STOPPED_1");
      stopped.durations.push_back(1.0);
      ordered_maps.push_back(stopped);
      represented_tips_by_leaf[terminal_old_id] = terminal->second;
    }
  }

  std::set<int> parent_nodes;
  std::set<int> all_nodes;
  // Classify old topology nodes from compact edges before public renumbering.
  for (const std::pair<int, int>& edge : old_edges) {
    parent_nodes.insert(edge.first);
    all_nodes.insert(edge.first);
    all_nodes.insert(edge.second);
  }
  std::vector<int> terminal_nodes;
  std::vector<int> internal_nodes;
  // Nodes without outgoing compact edges are public transition-tree tips; all
  // other nodes remain internal, including the synthetic root.
  for (int node_id : all_nodes) {
    // Leaf and internal node domains must be disjoint for valid ape numbering.
    if (!parent_nodes.count(node_id)) {
      terminal_nodes.push_back(node_id);
    } else {
      internal_nodes.push_back(node_id);
    }
  }
  std::map<int, int> node_map;
  // Assign contiguous tip ids first in deterministic old-id order.
  for (int index = 0;
       index < static_cast<int>(terminal_nodes.size());
       ++index) {
    node_map[terminal_nodes[index]] = index + 1;
  }
  // Ape-compatible helpers require the root to be the first internal id,
  // immediately after all public tips. Assign it before the remaining
  // deterministic internal-node order.
  if (!internal_nodes.empty()) {
    node_map[root_old_id] = static_cast<int>(terminal_nodes.size()) + 1;
  }
  int next_internal_id = static_cast<int>(terminal_nodes.size()) + 2;
  // Assign every nonroot internal node after the fixed public root id.
  for (int old_node_id : internal_nodes) {
    // The root was assigned above and must not consume a second public id.
    if (old_node_id == root_old_id) {
      continue;
    }
    node_map[old_node_id] = next_internal_id++;
  }

  Rcpp::IntegerMatrix edge(old_edges.size(), 2);
  Rcpp::NumericVector edge_length(old_edges.size());
  Rcpp::List maps(old_edges.size());
  Rcpp::NumericMatrix mapped_edge(
    old_edges.size(),
    expanded_state_space.size()
  );
  std::map<std::string, int> mapped_column;
  // Cache each declared state-size label's output column once.
  for (int column = 0; column < expanded_state_space.size(); ++column) {
    mapped_column[Rcpp::as<std::string>(expanded_state_space[column])] = column;
  }
  Rcpp::CharacterVector mapped_rownames(old_edges.size());
  // Serialize topology, maps, branch lengths, and mapped exposure in parallel.
  for (int edge_id = 0;
       edge_id < static_cast<int>(old_edges.size());
       ++edge_id) {
    int parent = node_map[old_edges[edge_id].first];
    int child = node_map[old_edges[edge_id].second];
    edge(edge_id, 0) = parent;
    edge(edge_id, 1) = child;
    maps[edge_id] = pst_v34_transition_map_to_vector(ordered_maps[edge_id]);
    mapped_rownames[edge_id] =
      std::to_string(parent) + "," + std::to_string(child);
    double length = 0.0;
    // Sum this compact edge's already-formatted unit/placeholder segments and
    // write their exposure under matching declared labels.
    for (int segment_id = 0;
         segment_id < static_cast<int>(ordered_maps[edge_id].labels.size());
         ++segment_id) {
      const std::string& label = ordered_maps[edge_id].labels[segment_id];
      double duration = ordered_maps[edge_id].durations[segment_id];
      length += duration;
      std::map<std::string, int>::const_iterator column =
        mapped_column.find(label);
      // Every compact label should belong to the declared expanded state
      // domain; unknown labels indicate support/formatter divergence.
      if (column == mapped_column.end()) {
        Rcpp::stop("V34 compact transition label is outside mapped.edge");
      }
      mapped_edge(edge_id, column->second) += duration;
    }
    edge_length[edge_id] = length;
  }

  Rcpp::CharacterVector public_tip_labels(terminal_nodes.size());
  std::map<std::string, std::vector<int> > represented_tip_ids_by_label;
  // Format each terminal history from traversal-owned represented tip ids.
  for (int tip_id = 0;
       tip_id < static_cast<int>(terminal_nodes.size());
       ++tip_id) {
    int old_node_id = terminal_nodes[tip_id];
    std::map<int, std::vector<int> >::const_iterator represented =
      represented_tips_by_leaf.find(old_node_id);
    // Every leaf must own biological or provisional represented tips.
    if (represented == represented_tips_by_leaf.end() ||
        represented->second.empty()) {
      Rcpp::stop("V34 compact transition leaf has no represented tips");
    }
    std::string label = pst_v34_join_tip_labels(
      represented->second,
      tip_labels
    );
    if (represented_tip_ids_by_label.count(label)) {
      Rcpp::stop("V34 compact transition tips have duplicate public labels");
    }
    represented_tip_ids_by_label[label] = represented->second;
    SET_STRING_ELT(
      public_tip_labels,
      tip_id,
      Rf_mkCharCE(label.c_str(), CE_UTF8)
    );
  }

  edge.attr("dimnames") = Rcpp::List::create(
    R_NilValue,
    Rcpp::CharacterVector::create("parent", "child")
  );
  mapped_edge.attr("dimnames") = Rcpp::List::create(
    mapped_rownames,
    expanded_state_space
  );
  Rcpp::List tree = Rcpp::List::create(
    Rcpp::Named("edge") = edge,
    Rcpp::Named("edge.length") = edge_length,
    Rcpp::Named("tip.label") = public_tip_labels,
    Rcpp::Named("Nnode") = static_cast<int>(internal_nodes.size()),
    Rcpp::Named("maps") = maps,
    Rcpp::Named("mapped.edge") = mapped_edge
  );
  tree.attr("class") = Rcpp::CharacterVector::create("phylo", "simmap");
  // Canonical paths form one edge per transition step, while the public T tree
  // stores uninterrupted unary path chains as multi-segment edges. Contracting
  // those topology-only nodes is deterministic formatting and preserves every
  // traversal-owned state, support value, represented tip, and path boundary.
  V34ScenarioTree compact_tree = pst_v34_scenario_tree_from_r(
    tree,
    time_tolerance
  );
  pst_v34_contract_unary_internal_edges_cpp(compact_tree, time_tolerance);
  Rcpp::List compact_output = pst_v34_scenario_tree_to_r(compact_tree);
  Rcpp::IntegerMatrix compact_edge = compact_output["edge"];
  Rcpp::List compact_maps = compact_output["maps"];
  Rcpp::NumericMatrix compact_mapped_edge(
    compact_edge.nrow(),
    expanded_state_space.size()
  );
  std::map<std::string, int> compact_column;
  // Restore the declared state-size column domain after topology compaction;
  // the generic scenario serializer intentionally emits an empty matrix.
  for (int column = 0; column < expanded_state_space.size(); ++column) {
    compact_column[Rcpp::as<std::string>(expanded_state_space[column])] = column;
  }
  Rcpp::CharacterVector compact_rownames(compact_edge.nrow());
  // Rebuild mapped exposure from the already-compacted edge maps. This is
  // numeric serialization only and does not infer path or transition identity.
  for (int edge_id = 0; edge_id < compact_edge.nrow(); ++edge_id) {
    compact_rownames[edge_id] =
      std::to_string(compact_edge(edge_id, 0)) + "," +
      std::to_string(compact_edge(edge_id, 1));
    Rcpp::NumericVector edge_map = compact_maps[edge_id];
    Rcpp::CharacterVector edge_labels = edge_map.names();
    // Copy every compact map duration into its predeclared state-size column.
    for (int segment_id = 0; segment_id < edge_map.size(); ++segment_id) {
      std::string label = Rcpp::as<std::string>(edge_labels[segment_id]);
      std::map<std::string, int>::const_iterator column = compact_column.find(label);
      // Compaction concatenates existing labels only; a missing column proves
      // the compact output diverged from traversal-owned support.
      if (column == compact_column.end()) {
        Rcpp::stop("V34 compact transition map lost its mapped-edge column");
      }
      compact_mapped_edge(edge_id, column->second) += edge_map[segment_id];
    }
  }
  compact_mapped_edge.attr("dimnames") = Rcpp::List::create(
    compact_rownames,
    expanded_state_space
  );
  compact_output["mapped.edge"] = compact_mapped_edge;
  Rcpp::CharacterVector compact_tip_labels = compact_output["tip.label"];
  Rcpp::List compact_represented_tip_ids(compact_tip_labels.size());
  // Restore the traversal-owned biological tip ids parallel to public
  // transition tips after unary contraction. Later size formatting uses these
  // ids to resolve the latest real terminal time without parsing merged labels.
  for (int tip_id = 0; tip_id < compact_tip_labels.size(); ++tip_id) {
    std::string label = Rcpp::as<std::string>(compact_tip_labels[tip_id]);
    std::map<std::string, std::vector<int> >::const_iterator represented =
      represented_tip_ids_by_label.find(label);
    if (represented == represented_tip_ids_by_label.end()) {
      Rcpp::stop("V34 compact transition tip lost represented biological ids");
    }
    compact_represented_tip_ids[tip_id] = Rcpp::wrap(represented->second);
  }
  compact_output.attr("pst_v34_represented_tip_ids") =
    compact_represented_tip_ids;
  return compact_output;
}

/**
 * Build the public transition tree from compact transition support.
 *
 * This is an Rcpp boundary: inputs are R objects and outputs must preserve documented V34 public shapes.
 */
Rcpp::List pst_v34_transition_tree_from_support_impl(
    Rcpp::List transition_support,
    Rcpp::CharacterVector tip_labels,
    double time_tolerance) {
  // V34 accepts only traversal-owned compact path records. Dense history
  // matrices require post-traversal regrouping and are an ownership failure.
  if (!transition_support.containsElementNamed("edge_path_id")) {
    Rcpp::stop("V34 transition support omitted compact path records");
  }
  SEXP support_digest = transition_support.attr(
    "pst_v34_transition_support_digest"
  );
  if (Rf_isNull(support_digest)) {
    Rcpp::stop("V34 transition support omitted its ownership digest");
  }
  Rcpp::List out = pst_v34_transition_tree_from_compact_support(
    transition_support,
    tip_labels,
    time_tolerance
  );
  out.attr("pst_v34_transition_support_digest") = support_digest;
  return out;

}

/**
 * Remove traversal-only STOPPED segments from a public transition tree.
 *
 * This is an Rcpp boundary: inputs are R objects and outputs must preserve documented V34 public shapes.
 */
Rcpp::List pst_v34_correct_stopped_branches_impl(
    Rcpp::List trans_tree,
    std::string stop_char,
    bool shrink_bl) {
  trans_tree = Rcpp::clone(trans_tree);
  Rcpp::IntegerMatrix edge = trans_tree["edge"];
  Rcpp::NumericVector edge_length = trans_tree["edge.length"];
  Rcpp::List maps = trans_tree["maps"];
  Rcpp::NumericMatrix mapped_edge = trans_tree["mapped.edge"];
  Rcpp::List mapped_dimnames = mapped_edge.attr("dimnames");
  Rcpp::CharacterVector mapped_rownames = mapped_dimnames[0];
  Rcpp::CharacterVector mapped_colnames = mapped_dimnames[1];

  std::vector<int> kept_cols;
  std::map<std::string, int> kept_col_index;
  // Retain only biological mapped-edge columns and build their new column index.
  for (int col = 0; col < mapped_colnames.size(); ++col) {
    std::string label = Rcpp::as<std::string>(mapped_colnames[col]);
    // Labels without the STOPPED marker remain part of the public transition state space.
    if (label.find(stop_char) == std::string::npos) {
      kept_col_index[label] = static_cast<int>(kept_cols.size());
      kept_cols.push_back(col);
    }
  }
  Rcpp::NumericMatrix cleaned_mapped(mapped_edge.nrow(), kept_cols.size());
  Rcpp::CharacterVector cleaned_colnames(kept_cols.size());
  // Copy retained mapped-edge columns without changing edge-row order.
  for (int new_col = 0; new_col < static_cast<int>(kept_cols.size()); ++new_col) {
    int old_col = kept_cols[new_col];
    cleaned_colnames[new_col] = mapped_colnames[old_col];
    // Preserve each biological edge's exposure under the compacted column domain.
    for (int row = 0; row < mapped_edge.nrow(); ++row) {
      cleaned_mapped(row, new_col) = mapped_edge(row, old_col);
    }
  }

  // Clean STOPPED segments independently on every transition-tree edge.
  for (int edge_id = 0; edge_id < maps.size(); ++edge_id) {
    Rcpp::NumericVector edge_map = maps[edge_id];
    Rcpp::CharacterVector edge_names = edge_map.names();
    std::vector<int> stop_indices;
    double stop_time = 0.0;
    // Locate all STOPPED segments and sum their nonbiological duration.
    for (int i = 0; i < edge_names.size(); ++i) {
      std::string label = Rcpp::as<std::string>(edge_names[i]);
      // Mark only labels carrying the reserved STOPPED token.
      if (label.find(stop_char) != std::string::npos) {
        stop_indices.push_back(i);
        stop_time += edge_map[i];
      }
    }
    // Edges without STOPPED exposure already satisfy the public contract.
    if (stop_indices.empty()) {
      continue;
    }

    // A mixed edge has a preceding biological state that can absorb or precede STOPPED time.
    if (edge_map.size() > 1) {
      int previous_index = edge_map.size() - 2;
      std::string last_state = Rcpp::as<std::string>(edge_names[previous_index]);
      // Shrink mode removes terminal bookkeeping time from physical branch length.
      if (shrink_bl) {
        edge_length[edge_id] = edge_length[edge_id] - stop_time;
      } else {
        // Nonshrink mode transfers STOPPED duration to the immediately preceding biological run.
        for (int stop_index : stop_indices) {
          int target = stop_index - 1;
          // Only a STOPPED segment with a preceding run can transfer its duration.
          if (target >= 0) {
            edge_map[target] = edge_map[target] + edge_map[stop_index];
          }
        }
        auto col_it = kept_col_index.find(last_state);
        // Mirror transferred duration into the preceding state's mapped-edge column.
        if (col_it != kept_col_index.end()) {
          cleaned_mapped(edge_id, col_it->second) = stop_time;
        }
      }

      int keep_count = edge_map.size() - static_cast<int>(stop_indices.size());
      Rcpp::NumericVector new_map(keep_count);
      Rcpp::CharacterVector new_names(keep_count);
      std::set<int> stop_set(stop_indices.begin(), stop_indices.end());
      int pos = 0;
      // Rebuild the named map from biological segments only.
      for (int i = 0; i < edge_map.size(); ++i) {
        // Exclude every traversal-only STOPPED segment from the public map.
        if (stop_set.count(i)) {
          continue;
        }
        new_map[pos] = edge_map[i];
        new_names[pos] = edge_names[i];
        ++pos;
      }
      new_map.attr("names") = new_names;
      maps[edge_id] = new_map;
    } else {
      // An all-STOPPED edge has no local biological predecessor; inherit the
      // terminal state of its parent edge when topology provides one.
      int parent_node = edge(edge_id, 0);
      std::string last_state;
      // Locate the incoming parent edge whose terminal state precedes this branch.
      for (int parent_edge = 0; parent_edge < edge.nrow(); ++parent_edge) {
        // The edge ending at this branch's parent supplies the inherited state.
        if (edge(parent_edge, 1) == parent_node) {
          Rcpp::NumericVector parent_map = maps[parent_edge];
          Rcpp::CharacterVector parent_names = parent_map.names();
          // A nonempty parent map exposes its terminal biological state.
          if (parent_names.size() > 0) {
            last_state = Rcpp::as<std::string>(parent_names[parent_names.size() - 1]);
          }
          break;
        }
      }
      // Root-level all-STOPPED edges fall back to their own encoded predecessor label.
      if (last_state.empty()) {
        last_state = Rcpp::as<std::string>(edge_names[0]);
      }
      Rcpp::NumericVector new_map(1);
      Rcpp::CharacterVector new_names(1);
      new_map[0] = stop_time;
      new_names[0] = last_state;
      new_map.attr("names") = new_names;
      maps[edge_id] = new_map;
      auto col_it = kept_col_index.find(last_state);
      // Attribute transferred duration only when the inherited state is in the public domain.
      if (col_it != kept_col_index.end()) {
        cleaned_mapped(edge_id, col_it->second) += stop_time;
      }
      // Shrink mode removes the entire all-STOPPED branch exposure.
      if (shrink_bl) {
        edge_length[edge_id] = edge_length[edge_id] - stop_time;
        new_map[0] = 0.0;
        maps[edge_id] = new_map;
        // Clear the inherited-state mapped exposure together with the shrunken map.
        if (col_it != kept_col_index.end()) {
          cleaned_mapped(edge_id, col_it->second) = 0.0;
        }
      }
    }
  }

  cleaned_mapped.attr("dimnames") = Rcpp::List::create(mapped_rownames, cleaned_colnames);
  trans_tree["edge.length"] = edge_length;
  trans_tree["maps"] = maps;
  trans_tree["mapped.edge"] = cleaned_mapped;
  return trans_tree;
}

/**
 * Strip V34 size suffixes while preserving public simmap topology and exposure.
 *
 * `tree` is cloned before mutation. Scenario-size source labels are reduced to
 * their biological states, while already-state-only source labels are copied
 * exactly. Adjacent equal states are merged in edge-local time order, and
 * `mapped.edge` is rebuilt from the resulting maps. This formatter owns V34
 * label generation at the C++ output boundary; malformed map/topology shapes
 * fail before any partial object is returned to R.
 */
Rcpp::List pst_v34_state_only_tree_impl(
    Rcpp::List tree,
    std::string summary_suffix,
    Rcpp::Nullable<Rcpp::CharacterVector> state_levels,
    bool include_edge_boundaries,
    bool source_labels_are_state_only) {
  tree = Rcpp::clone(tree);
  Rcpp::IntegerMatrix edge = tree["edge"];
  Rcpp::List maps = tree["maps"];
  int nedge = edge.nrow();
  bool has_path_maps = tree.containsElementNamed("path_maps");
  Rcpp::List path_maps = has_path_maps ?
    Rcpp::as<Rcpp::List>(tree["path_maps"]) : Rcpp::List();
  if (has_path_maps && path_maps.size() != nedge) {
    Rcpp::stop("V34 state-only formatter received misaligned path maps");
  }

  // A simmap tree must keep one map per public edge; otherwise state exposure
  // cannot be aligned to topology rows without guessing.
  if (maps.size() != nedge) {
    Rcpp::stop("V34 state-only tree formatter received misaligned edge maps");
  }

  Rcpp::List state_maps(nedge);
  Rcpp::List compact_path_maps(nedge);
  std::set<std::string> state_domain;
  // Rewrite every edge map from size/path labels to biological state labels.
  // The loop reads the public edge-order map list and writes an edge-parallel
  // state-only map list; after each iteration, `state_maps[edge_id]` is aligned
  // with `edge[edge_id, ]` and all emitted states are present in `state_domain`.
  for (int edge_id = 0; edge_id < nedge; ++edge_id) {
    Rcpp::NumericVector edge_map = maps[edge_id];
    Rcpp::CharacterVector edge_labels = edge_map.names();
    Rcpp::NumericVector edge_path_map;
    Rcpp::CharacterVector edge_path_labels;
    if (has_path_maps) {
      edge_path_map = Rcpp::as<Rcpp::NumericVector>(path_maps[edge_id]);
      edge_path_labels = edge_path_map.names();
      if (edge_path_map.size() != edge_map.size()) {
        Rcpp::stop("V34 state-only formatter found nonparallel path segments");
      }
    }

    // Every duration must have a label because the state-only formatter cannot
    // infer biological state from a bare numeric segment.
    if (edge_labels.size() != edge_map.size()) {
      Rcpp::stop("V34 state-only tree formatter received an unnamed map segment");
    }

    std::vector<std::string> compact_labels;
    std::vector<double> compact_durations;
    std::vector<std::string> compact_path_labels;
    // Process mapped segments in chronological order on this edge. The loop
    // reads original labelled durations and writes merged biological-state
    // runs; after each iteration, adjacent equal states have one accumulated
    // duration and no segment has changed its chronological position.
    for (int segment_id = 0; segment_id < edge_map.size(); ++segment_id) {
      // Missing segment labels are malformed construction output; treating
      // them as text would silently create a bogus public state.
      if (edge_labels[segment_id] == NA_STRING) {
        Rcpp::stop("V34 state-only tree formatter received an NA map label");
      }
      std::string source_label = Rcpp::as<std::string>(edge_labels[segment_id]);
      std::string state_label = source_labels_are_state_only ?
        source_label : pst_v34_scenario_state(source_label);
      double duration = edge_map[segment_id];
      std::string path_label = has_path_maps ?
        Rcpp::as<std::string>(edge_path_labels[segment_id]) : std::string();

      // Encoded state-size sources may collapse after the size suffix is
      // removed. Already state-only sources are traversal/public segments and
      // retain every boundary, including adjacent equal labels.
      if (!source_labels_are_state_only &&
          !compact_labels.empty() &&
          compact_labels.back() == state_label &&
          (!has_path_maps || compact_path_labels.back() == path_label)) {
        compact_durations.back() += duration;
      } else {
        // A changed biological state starts a new public mapped segment.
        compact_labels.push_back(state_label);
        compact_durations.push_back(duration);
        if (has_path_maps) {
          compact_path_labels.push_back(path_label);
        }
      }
      state_domain.insert(state_label);
    }

    Rcpp::NumericVector out_map(compact_durations.size());
    Rcpp::CharacterVector out_labels(compact_labels.size());
    // Serialize the compact edge-local state runs back into an R named vector.
    // The loop writes labels and durations positionally so map names and values
    // remain a single biological segment stream.
    for (int segment_id = 0; segment_id < static_cast<int>(compact_durations.size()); ++segment_id) {
      out_map[segment_id] = compact_durations[segment_id];
      out_labels[segment_id] = compact_labels[segment_id];
    }
    out_map.attr("names") = out_labels;
    state_maps[edge_id] = out_map;
    if (has_path_maps) {
      Rcpp::NumericVector out_path_map(compact_durations.size());
      Rcpp::CharacterVector out_path_labels(compact_path_labels.size());
      for (int segment_id = 0;
           segment_id < static_cast<int>(compact_durations.size());
           ++segment_id) {
        out_path_map[segment_id] = compact_durations[segment_id];
        out_path_labels[segment_id] = compact_path_labels[segment_id];
      }
      out_path_map.attr("names") = out_path_labels;
      compact_path_maps[edge_id] = out_path_map;
    }
  }

  std::vector<std::string> states(state_domain.begin(), state_domain.end());
  std::map<std::string, int> state_col;
  // Cache the deterministic mapped-edge column for every biological state that
  // appears in the rewritten maps.
  for (int state_id = 0; state_id < static_cast<int>(states.size()); ++state_id) {
    state_col[states[state_id]] = state_id;
  }

  Rcpp::NumericMatrix mapped_edge(nedge, states.size());
  Rcpp::CharacterVector mapped_rownames(nedge);
  Rcpp::CharacterVector mapped_colnames(states.size());
  // Label mapped-edge rows by the unchanged public topology endpoints.
  for (int edge_id = 0; edge_id < nedge; ++edge_id) {
    mapped_rownames[edge_id] =
      std::to_string(edge(edge_id, 0)) + "," + std::to_string(edge(edge_id, 1));
  }
  // Emit mapped-edge columns in deterministic lexical state order.
  for (int state_id = 0; state_id < static_cast<int>(states.size()); ++state_id) {
    mapped_colnames[state_id] = states[state_id];
  }
  // Accumulate state exposure from the native state-only maps. The loop reads
  // each edge's compact map and writes exactly one mapped-edge row, preserving
  // the invariant that row sums equal edge-map sums.
  for (int edge_id = 0; edge_id < nedge; ++edge_id) {
    Rcpp::NumericVector edge_map = state_maps[edge_id];
    Rcpp::CharacterVector edge_labels = edge_map.names();
    // Visit map segments positionally so every duration is credited to its own state column.
    for (int segment_id = 0; segment_id < edge_map.size(); ++segment_id) {
      std::string state_label = Rcpp::as<std::string>(edge_labels[segment_id]);
      auto col_it = state_col.find(state_label);
      // The state domain was built from these maps; missing columns indicate
      // internal corruption and must not be repaired by dropping exposure.
      if (col_it == state_col.end()) {
        Rcpp::stop("V34 state-only tree formatter lost a mapped state column");
      }
      mapped_edge(edge_id, col_it->second) += edge_map[segment_id];
    }
  }
  mapped_edge.attr("dimnames") = Rcpp::List::create(mapped_rownames, mapped_colnames);

  tree["maps"] = state_maps;
  if (has_path_maps) {
    tree["path_maps"] = compact_path_maps;
  }
  tree["mapped.edge"] = mapped_edge;
  tree.attr("class") = Rcpp::CharacterVector::create("simmap", "phylo");
  // Materialization-time SS parts are optional because debug and intermediate
  // callers may need only the state-only tree. When requested, attach compact
  // metrics before the public tree leaves the native output boundary.
  if (!summary_suffix.empty()) {
    // Summary capture needs the same biological state domain used by the final
    // compatibility SS tables; without it, row/column ordering would become a
    // caller-side reconstruction decision.
    if (state_levels.isNull()) {
      Rcpp::stop("V34 state-only summary capture requires state levels");
    }
    Rcpp::CharacterVector summary_states(state_levels);
    tree.attr("pst_v34_tree_summary_parts") =
      pst_v34_capture_tree_summary_parts_cpp(
        tree,
        summary_suffix,
        summary_states,
        include_edge_boundaries
      );
  }
  return tree;
}

/**
 * Relabel a state tree into root-to-state path labels.
 *
 * This is an Rcpp boundary: inputs are R objects and outputs must preserve documented V34 public shapes.
 */
Rcpp::List pst_v34_pathify_state_tree_impl(
    Rcpp::List tree,
    Rcpp::List root_policy,
    std::string sep) {
  tree = Rcpp::clone(tree);
  Rcpp::IntegerMatrix edge = tree["edge"];
  Rcpp::List maps = tree["maps"];
  Rcpp::CharacterVector tip_labels = tree["tip.label"];
  int ntips = tip_labels.size();
  int nedge = edge.nrow();
  int root = ntips + 1;

  std::set<std::string> state_set;
  // Collect the complete biological state domain from all input mapped segments.
  for (int edge_id = 0; edge_id < maps.size(); ++edge_id) {
    Rcpp::NumericVector edge_map = maps[edge_id];
    Rcpp::CharacterVector names = edge_map.names();
    // Add each nonmissing state label once to the canonical domain.
    for (int i = 0; i < names.size(); ++i) {
      // Missing labels do not define a biological state or path transition.
      if (names[i] != NA_STRING) {
        state_set.insert(Rcpp::as<std::string>(names[i]));
      }
    }
  }

  std::vector<std::string> states(state_set.begin(), state_set.end());
  std::map<std::string, int> state_ids;
  // Assign deterministic compact ids to lexical state labels for path-registry lookup.
  for (int i = 0; i < static_cast<int>(states.size()); ++i) {
    state_ids[states[i]] = i + 1;
  }
  int n_state = static_cast<int>(states.size());

  Rcpp::CharacterVector anchor = root_policy["anchor"];
  std::vector<std::string> root_components;
  // Seed the root path with every explicit root-policy anchor state in order.
  for (int i = 0; i < anchor.size(); ++i) {
    // Missing anchor entries do not contribute a biological root-history step.
    if (anchor[i] != NA_STRING) {
      root_components.push_back(Rcpp::as<std::string>(anchor[i]));
    }
  }

  std::vector<std::vector<std::string> > path_components(1);
  std::vector<std::string> path_labels(1);
  std::vector<int> path_parent(1, NA_INTEGER);
  std::vector<int> path_added_state_id(1, NA_INTEGER);
  path_components.push_back(root_components);
  path_labels.push_back(pst_v34_path_label(root_components, sep));
  path_parent.push_back(NA_INTEGER);
  path_added_state_id.push_back(NA_INTEGER);
  std::vector<std::vector<int> > transition_to_path(2, std::vector<int>(n_state + 1, 0));

  int max_node = root;
  // Bound node ids before allocating per-node path and terminal-state ownership.
  for (int i = 0; i < nedge; ++i) {
    max_node = std::max(max_node, edge(i, 0));
    max_node = std::max(max_node, edge(i, 1));
  }
  std::vector<int> path_id_at_node(max_node + 1, 0);
  std::vector<int> tail_state_id_at_node(max_node + 1, NA_INTEGER);
  path_id_at_node[root] = 1;

  // A known final anchor state establishes the root path's current terminal state.
  if (anchor.size() > 0 && anchor[anchor.size() - 1] != NA_STRING) {
    std::string root_tail = Rcpp::as<std::string>(anchor[anchor.size() - 1]);
    auto tail_it = state_ids.find(root_tail);
    // Set root state ownership only when the anchor belongs to the observed state domain.
    if (tail_it != state_ids.end()) {
      tail_state_id_at_node[root] = tail_it->second;
    }
  }

  Rcpp::List new_maps(nedge);
  Rcpp::List component_maps(nedge);
  Rcpp::List path_id_maps(nedge);
  std::vector<int> ordered = pst_v34_preorder_edges(edge, root);
  // Traverse edges root-first so every parent node's path is known before its daughters.
  for (int edge_id : ordered) {
    int row = edge_id - 1;
    int parent = edge(row, 0);
    int child = edge(row, 1);
    int current_path_id = path_id_at_node[parent];
    int current_tail_state_id = tail_state_id_at_node[parent];
    // Missing parent path ownership indicates broken topology or path propagation and cannot be formatted safely.
    if (current_path_id == 0 || current_path_id == NA_INTEGER) {
      Rcpp::stop("Missing V34 path id at parent node %d", parent);
    }

    Rcpp::NumericVector edge_map = maps[row];
    Rcpp::CharacterVector edge_states = edge_map.names();
    Rcpp::CharacterVector labels(edge_map.size());
    Rcpp::IntegerVector segment_path_ids(edge_map.size());
    Rcpp::List segment_components(edge_map.size());

    // Advance through mapped states in biological time order and annotate each segment with its path.
    for (int segment_id = 0; segment_id < edge_map.size(); ++segment_id) {
      std::string state = Rcpp::as<std::string>(edge_states[segment_id]);
      int state_id = state_ids[state];
      // Only a real state change extends the transition sequence; continuations retain the current path.
      if (current_tail_state_id != state_id) {
        int next_path_id = transition_to_path[current_path_id][state_id];
        // Reuse a previously registered child path for the same parent-path/state transition.
        if (next_path_id > 0) {
          current_path_id = next_path_id;
        } else {
          // First arrival to this parent-path/state pair creates one canonical child path.
          int parent_path_id = current_path_id;
          current_path_id = static_cast<int>(path_components.size());
          std::vector<std::string> next_components = path_components[parent_path_id];
          next_components.push_back(state);
          path_components.push_back(next_components);
          path_labels.push_back(pst_v34_path_label(next_components, sep));
          path_parent.push_back(parent_path_id);
          path_added_state_id.push_back(state_id);
          transition_to_path[parent_path_id][state_id] = current_path_id;
          transition_to_path.push_back(std::vector<int>(n_state + 1, 0));
        }
        current_tail_state_id = state_id;
      }

      labels[segment_id] = path_labels[current_path_id];
      segment_path_ids[segment_id] = current_path_id;
      segment_components[segment_id] = Rcpp::wrap(path_components[current_path_id]);
    }

    Rcpp::NumericVector out_map = Rcpp::clone(edge_map);
    out_map.attr("names") = labels;
    new_maps[row] = out_map;
    path_id_maps[row] = segment_path_ids;
    component_maps[row] = segment_components;
    path_id_at_node[child] = current_path_id;
    tail_state_id_at_node[child] = current_tail_state_id;
  }

  std::set<std::string> path_set;
  // Collect every path label that occupies positive mapped-tree exposure.
  for (int edge_id = 0; edge_id < nedge; ++edge_id) {
    Rcpp::NumericVector edge_map = new_maps[edge_id];
    Rcpp::CharacterVector names = edge_map.names();
    // Add nonmissing segment paths to the mapped-edge column domain.
    for (int i = 0; i < names.size(); ++i) {
      // Missing labels do not define a public transition sequence.
      if (names[i] != NA_STRING) {
        path_set.insert(Rcpp::as<std::string>(names[i]));
      }
    }
  }
  std::vector<std::string> paths(path_set.begin(), path_set.end());
  std::map<std::string, int> path_col;
  // Cache each lexical path label's mapped-edge column.
  for (int i = 0; i < static_cast<int>(paths.size()); ++i) {
    path_col[paths[i]] = i;
  }

  Rcpp::NumericMatrix mapped_edge(nedge, paths.size());
  Rcpp::CharacterVector mapped_rownames(nedge);
  Rcpp::CharacterVector mapped_colnames(paths.size());
  // Label mapped-edge rows by the unchanged public topology endpoints.
  for (int i = 0; i < nedge; ++i) {
    mapped_rownames[i] = std::to_string(edge(i, 0)) + "," + std::to_string(edge(i, 1));
  }
  // Emit mapped-edge columns in deterministic lexical path order.
  for (int i = 0; i < static_cast<int>(paths.size()); ++i) {
    mapped_colnames[i] = paths[i];
  }
  // Accumulate each edge's duration by canonical path label.
  for (int edge_id = 0; edge_id < nedge; ++edge_id) {
    Rcpp::NumericVector edge_map = new_maps[edge_id];
    Rcpp::CharacterVector names = edge_map.names();
    // Visit segment labels and durations together so path exposure remains exact.
    for (int i = 0; i < edge_map.size(); ++i) {
      auto col_it = path_col.find(Rcpp::as<std::string>(names[i]));
      // Only paths declared in the public column domain receive mapped exposure.
      if (col_it != path_col.end()) {
        mapped_edge(edge_id, col_it->second) += edge_map[i];
      }
    }
  }
  mapped_edge.attr("dimnames") = Rcpp::List::create(mapped_rownames, mapped_colnames);

  Rcpp::List dictionary;
  // Build the path dictionary directly from each labelled segment's registered components.
  for (int edge_id = 0; edge_id < nedge; ++edge_id) {
    Rcpp::NumericVector edge_map = new_maps[edge_id];
    Rcpp::CharacterVector names = edge_map.names();
    Rcpp::List segments = component_maps[edge_id];
    // Associate every public path label with its root-to-state component vector.
    for (int segment_id = 0; segment_id < edge_map.size(); ++segment_id) {
      dictionary[Rcpp::as<std::string>(names[segment_id])] = segments[segment_id];
    }
  }

  Rcpp::IntegerVector lookup_path_id(path_components.size() - 1);
  Rcpp::CharacterVector lookup_label(path_components.size() - 1);
  Rcpp::IntegerVector lookup_parent(path_components.size() - 1);
  Rcpp::IntegerVector lookup_added(path_components.size() - 1);
  // Serialize the native path registry into the public lookup table.
  for (int path_id = 1; path_id < static_cast<int>(path_components.size()); ++path_id) {
    int pos = path_id - 1;
    lookup_path_id[pos] = path_id;
    lookup_label[pos] = path_labels[path_id];
    lookup_parent[pos] = path_parent[path_id];
    lookup_added[pos] = path_added_state_id[path_id];
  }
  Rcpp::DataFrame lookup = Rcpp::DataFrame::create(
    Rcpp::Named("path_id") = lookup_path_id,
    Rcpp::Named("label") = lookup_label,
    Rcpp::Named("parent_path_id") = lookup_parent,
    Rcpp::Named("added_state_id") = lookup_added
  );

  tree["maps"] = new_maps;
  tree["mapped.edge"] = mapped_edge;
  tree.attr("class") = Rcpp::CharacterVector::create("simmap", "phylo");
  tree.attr("pst_root_policy") = root_policy;
  tree.attr("pst_path_components") = component_maps;
  tree.attr("pst_path_id_maps") = path_id_maps;
  tree.attr("pst_path_dictionary") = dictionary;
  tree.attr("pst_v4_path_lookup") = lookup;
  return tree;
}

/**
 * Materialize transition-tree size maps from committed online path support.
 *
 * `state_tree` and `path_tree` must expose identical topology and parallel map
 * segments. Final cumulative rows own internal/through sizes, while terminal
 * contributions committed with LTT/STT batches own terminal-edge tip sizes.
 * The function returns fresh map lists without mutating either tree. Missing
 * or duplicate path labels, mismatched topology, and malformed TT matrices
 * fail before a partial size-map family crosses the Rcpp boundary.
 */
Rcpp::List pst_v34_transition_size_maps_impl(
    Rcpp::List state_tree,
    Rcpp::List path_tree,
    Rcpp::NumericMatrix cumulative_lineage_trans,
    Rcpp::NumericMatrix cumulative_scenario_trans,
    Rcpp::CharacterVector path_labels_by_id,
    Rcpp::IntegerVector terminal_lineage_counts_by_id,
    Rcpp::IntegerVector terminal_scenario_counts_by_id,
    Rcpp::IntegerVector terminal_trajectory_group_id_by_tip) {
  Rcpp::IntegerMatrix state_edge = state_tree["edge"];
  Rcpp::IntegerMatrix path_edge = path_tree["edge"];
  Rcpp::List state_maps = state_tree["maps"];
  Rcpp::List path_maps = path_tree["maps"];
  int nedge = state_edge.nrow();

  // Tree projections must remain topology-parallel because each numeric size
  // is written to the same biological segment in both state and path views.
  if (path_edge.nrow() != nedge ||
      state_maps.size() != nedge ||
      path_maps.size() != nedge) {
    Rcpp::stop("V34 transition size-map trees are not edge-parallel");
  }
  // Verify every topology endpoint before segment-level values are allocated;
  // a reordered or changed path projection cannot safely share map positions.
  for (int edge_id = 0; edge_id < nedge; ++edge_id) {
    // Different endpoints prove the two projections no longer identify the
    // same public transition-tree edge.
    if (state_edge(edge_id, 0) != path_edge(edge_id, 0) ||
        state_edge(edge_id, 1) != path_edge(edge_id, 1)) {
      Rcpp::stop("V34 transition size-map trees have different topology");
    }
  }

  {
    if (cumulative_lineage_trans.nrow() < 1 ||
        cumulative_scenario_trans.nrow() < 1) {
      Rcpp::stop("V34 cumulative size-map matrices have no committed rows");
    }
    Rcpp::CharacterVector lineage_names =
      Rcpp::colnames(cumulative_lineage_trans);
    Rcpp::CharacterVector scenario_names =
      Rcpp::colnames(cumulative_scenario_trans);
    std::map<std::string, double> lineage_by_path;
    std::map<std::string, double> scenario_by_path;
    int lineage_final_row = cumulative_lineage_trans.nrow() - 1;
    int scenario_final_row = cumulative_scenario_trans.nrow() - 1;
    // Index the final cumulative value of every canonical lineage path.
    for (int column = 0; column < cumulative_lineage_trans.ncol(); ++column) {
      if (lineage_names[column] == NA_STRING) {
        Rcpp::stop("V34 cumulative lineage trans contains an unnamed column");
      }
      std::string label = Rcpp::as<std::string>(lineage_names[column]);
      if (label == "total" || label == "summary_time_vec") {
        continue;
      }
      if (lineage_by_path.count(label)) {
        Rcpp::stop("V34 cumulative lineage trans contains duplicate paths");
      }
      lineage_by_path[label] = cumulative_lineage_trans(
        lineage_final_row,
        column
      );
    }
    // Index the final cumulative value of every canonical scenario path.
    for (int column = 0; column < cumulative_scenario_trans.ncol(); ++column) {
      if (scenario_names[column] == NA_STRING) {
        Rcpp::stop("V34 cumulative scenario trans contains an unnamed column");
      }
      std::string label = Rcpp::as<std::string>(scenario_names[column]);
      if (label == "total" || label == "summary_time_vec") {
        continue;
      }
      if (scenario_by_path.count(label)) {
        Rcpp::stop("V34 cumulative scenario trans contains duplicate paths");
      }
      scenario_by_path[label] = cumulative_scenario_trans(
        scenario_final_row,
        column
      );
    }
    if (lineage_by_path.size() != scenario_by_path.size()) {
      Rcpp::stop("V34 cumulative lineage/scenario trans paths are not parallel");
    }
    for (std::map<std::string, double>::const_iterator path =
           lineage_by_path.begin();
         path != lineage_by_path.end();
         ++path) {
      if (!scenario_by_path.count(path->first)) {
        Rcpp::stop("V34 cumulative lineage/scenario trans paths are unmatched");
      }
    }

    if (path_labels_by_id.size() < 1) {
      Rcpp::stop("V34 terminal size support has no canonical path labels");
    }
    std::map<std::string, double> terminal_lineage_by_path;
    std::map<std::string, double> terminal_scenario_by_path;
    auto index_terminal_counts = [&] (
        const Rcpp::IntegerVector& counts,
        const std::string& family,
        std::map<std::string, double>& values_by_path) {
      Rcpp::CharacterVector count_names = counts.names();
      if (counts.size() != count_names.size()) {
        Rcpp::stop("V34 terminal " + family + " counts lack path ids");
      }
      // Resolve each sparse positive terminal contribution through the shared
      // C++ path registry. No public tree or active path-wide TT row is scanned.
      for (int index = 0; index < counts.size(); ++index) {
        if (count_names[index] == NA_STRING || counts[index] < 0) {
          Rcpp::stop("V34 terminal " + family + " count is invalid");
        }
        std::string path_id_text = Rcpp::as<std::string>(count_names[index]);
        char* end = nullptr;
        long path_id = std::strtol(path_id_text.c_str(), &end, 10);
        if (end == path_id_text.c_str() || *end != '\0' ||
            path_id < 1 || path_id > path_labels_by_id.size()) {
          Rcpp::stop("V34 terminal " + family + " count has an unknown path id");
        }
        Rcpp::String label_value = path_labels_by_id[path_id - 1];
        if (label_value == NA_STRING) {
          Rcpp::stop("V34 terminal " + family + " path has no label");
        }
        std::string label = static_cast<std::string>(label_value);
        if (values_by_path.count(label)) {
          Rcpp::stop("V34 terminal " + family + " counts duplicate a path");
        }
        values_by_path[label] = counts[index];
      }
    };
    index_terminal_counts(
      terminal_lineage_counts_by_id,
      "lineage",
      terminal_lineage_by_path
    );
    index_terminal_counts(
      terminal_scenario_counts_by_id,
      "scenario",
      terminal_scenario_by_path
    );

    Rcpp::List scenario_size_maps(nedge);
    Rcpp::List lineage_size_maps(nedge);
    Rcpp::List scenario_size_path_maps(nedge);
    Rcpp::List lineage_size_path_maps(nedge);
    int ntip = Rcpp::CharacterVector(state_tree["tip.label"]).size();
    Rcpp::List represented_tip_ids = state_tree.attr(
      "pst_v34_represented_tip_ids"
    );
    if (represented_tip_ids.size() != ntip ||
        terminal_trajectory_group_id_by_tip.size() < 1) {
      Rcpp::stop("V34 terminal transition support is not tip-aligned");
    }
    // Assign cumulative support to internal segments and terminal counts to the
    // final segment of every tip edge. State/path projections share values.
    for (int edge_id = 0; edge_id < nedge; ++edge_id) {
      Rcpp::NumericVector state_map = state_maps[edge_id];
      Rcpp::NumericVector path_map = path_maps[edge_id];
      Rcpp::CharacterVector state_labels = state_map.names();
      Rcpp::CharacterVector path_labels = path_map.names();
      if (state_map.size() != path_map.size() ||
          state_labels.size() != state_map.size() ||
          path_labels.size() != path_map.size()) {
        Rcpp::stop("V34 transition cumulative size-map segments are not parallel");
      }
      Rcpp::NumericVector lineage_state(state_map.size());
      Rcpp::NumericVector scenario_state(state_map.size());
      Rcpp::NumericVector lineage_path(path_map.size());
      Rcpp::NumericVector scenario_path(path_map.size());
      lineage_state.attr("names") = Rcpp::clone(state_labels);
      scenario_state.attr("names") = Rcpp::clone(state_labels);
      lineage_path.attr("names") = Rcpp::clone(path_labels);
      scenario_path.attr("names") = Rcpp::clone(path_labels);
      bool terminal_edge = state_edge(edge_id, 1) <= ntip;
      for (int segment_id = 0; segment_id < path_map.size(); ++segment_id) {
        if (path_labels[segment_id] == NA_STRING) {
          Rcpp::stop("V34 transition cumulative size map has a missing path");
        }
        std::string path_label =
          Rcpp::as<std::string>(path_labels[segment_id]);
        std::map<std::string, double>::const_iterator lineage_value =
          lineage_by_path.find(path_label);
        std::map<std::string, double>::const_iterator scenario_value =
          scenario_by_path.find(path_label);
        if (lineage_value == lineage_by_path.end() ||
            scenario_value == scenario_by_path.end()) {
          Rcpp::stop("V34 transition cumulative size map has an unmatched path");
        }
        bool terminal_segment = terminal_edge &&
          segment_id == path_map.size() - 1;
        double lineage_size = lineage_value->second;
        double scenario_size = scenario_value->second;
        // The last segment of a tip edge records only lineage/scenario support
        // that terminates on this exact canonical path. Copying a path-wide
        // active LTT/STT cell would include unrelated histories sharing the
        // path at the same batch timestamp.
        if (terminal_segment) {
          std::map<std::string, double>::const_iterator terminal_lineage =
            terminal_lineage_by_path.find(path_label);
          std::map<std::string, double>::const_iterator terminal_scenario =
            terminal_scenario_by_path.find(path_label);
          lineage_size = terminal_lineage == terminal_lineage_by_path.end() ?
            0.0 : terminal_lineage->second;
          scenario_size = terminal_scenario == terminal_scenario_by_path.end() ?
            0.0 : terminal_scenario->second;

          int terminal_tip_node = state_edge(edge_id, 1);
          if (terminal_tip_node < 1 || terminal_tip_node > ntip) {
            Rcpp::stop("V34 terminal transition edge has an invalid tip node");
          }
          Rcpp::IntegerVector represented = represented_tip_ids[
            terminal_tip_node - 1
          ];
          if (represented.size() < 1) {
            Rcpp::stop("V34 terminal transition edge represents no phylogeny tips");
          }
          std::set<int> terminal_groups;
          for (int represented_id : represented) {
            if (represented_id < 1 ||
                represented_id > terminal_trajectory_group_id_by_tip.size()) {
              Rcpp::stop("V34 represented terminal tip id is out of range");
            }
            int group_id = terminal_trajectory_group_id_by_tip[
              represented_id - 1
            ];
            if (group_id < 1) {
              Rcpp::stop("V34 represented terminal tip has no trajectory group");
            }
            terminal_groups.insert(group_id);
          }
          // Terminal transition sizes belong to this exact public tip group,
          // not every terminal history that happens to share its final path.
          lineage_size = represented.size();
          scenario_size = terminal_groups.size();
        }
        lineage_state[segment_id] = lineage_size;
        lineage_path[segment_id] = lineage_size;
        scenario_state[segment_id] = scenario_size;
        scenario_path[segment_id] = scenario_size;
      }
      scenario_size_maps[edge_id] = scenario_state;
      lineage_size_maps[edge_id] = lineage_state;
      scenario_size_path_maps[edge_id] = scenario_path;
      lineage_size_path_maps[edge_id] = lineage_path;
    }
    return Rcpp::List::create(
      Rcpp::Named("scenario_size_maps") = scenario_size_maps,
      Rcpp::Named("lineage_size_maps") = lineage_size_maps,
      Rcpp::Named("scenario_size_path_maps") = scenario_size_path_maps,
      Rcpp::Named("lineage_size_path_maps") = lineage_size_path_maps,
      Rcpp::Named("construction_phase") =
        "native_tree_output_from_cumulative_through_and_terminal_tt"
    );
  }

}
