#include "pst_v35_tree_output.hpp"

#include "pst_v35_batch.hpp"
#include "pst_v35_ss_output.hpp"

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

struct V35TransitionMap {
  std::vector<std::string> labels;
  std::vector<double> durations;
};

struct V35ScenarioEdge {
  int parent;
  int child;
  V35TransitionMap map;
  V35TransitionMap path_map;
  std::set<int> raw_edge_ids;
};

struct V35ScenarioTree {
  std::vector<V35ScenarioEdge> edges;
  std::vector<std::string> tip_labels;
  int nnode;
  bool has_path_maps = false;
};

/**
 * Join represented tip labels in stable public tip order after scenario-tip merging.
 *
 * Inputs are read only and failures stop before a partially formatted public tree is returned.
 */
static std::string pst_v35_join_tip_labels(
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
static V35TransitionMap pst_v35_make_transition_map(
    const std::vector<int>& state_ids,
    const std::vector<int>& lineage_counts,
    const std::vector<double>& durations,
    const std::map<int, std::string>& state_labels) {
  V35TransitionMap out;
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
static Rcpp::NumericVector pst_v35_transition_map_to_vector(const V35TransitionMap& edge_map) {
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
static Rcpp::NumericVector pst_v35_scenario_edge_map_to_vector(const V35TransitionMap& edge_map) {
  return pst_v35_transition_map_to_vector(edge_map);
}

/**
 * Create an internal scenario-map segment vector from states, sizes, and durations.
 *
 * Inputs are read only and failures stop before a partially formatted public tree is returned.
 */
static std::pair<V35TransitionMap, V35TransitionMap> pst_v35_make_scenario_maps(
    const Rcpp::IntegerVector& state_ids,
    const Rcpp::IntegerVector& path_ids,
    const Rcpp::IntegerVector& sizes,
    const Rcpp::NumericVector& durations,
    const std::map<int, std::string>& state_labels,
    const std::map<int, std::string>& path_labels) {
  V35TransitionMap state_out;
  V35TransitionMap path_out;
  int n = durations.size();
  if (state_ids.size() != n ||
      path_ids.size() != n ||
      sizes.size() != n ||
      n < 1) {
    Rcpp::stop("V35 scenario record vectors are not nonempty and parallel");
  }
  // Convert the committed state, size, and duration vectors in their original
  // segment order so topology and scenario support remain aligned.
  for (int i = 0; i < n; ++i) {
    double duration = durations[i];
    if (!R_finite(duration) || duration < 0.0) {
      Rcpp::stop("V35 scenario segment duration is nonfinite or negative");
    }
    int state_id = state_ids[i];
    int path_id = path_ids[i];
    int size = sizes[i];
    auto label_it = state_labels.find(state_id);
    std::string state = label_it == state_labels.end() ? std::string("NA") : label_it->second;
    std::string label = state + "_" + std::to_string(size);
    std::map<int, std::string>::const_iterator path_it = path_labels.find(path_id);
    if (path_it == path_labels.end()) {
      Rcpp::stop("V35 scenario segment contains an unregistered path id");
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
static std::string pst_v35_escape_path_state(const std::string& state) {
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
static std::string pst_v35_path_label(
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
    out += pst_v35_escape_path_state(components[i]);
  }
  return out;
}

/**
 * Return phylogeny edge ids in stable root-first order for tree serialization.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static std::vector<int> pst_v35_preorder_edges(const Rcpp::IntegerMatrix& edge, int root) {
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
static std::vector<std::string> pst_v35_split_csv(const std::string& label) {
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
static std::string pst_v35_join_csv(const std::vector<std::string>& labels) {
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
static std::string pst_v35_scenario_state(const std::string& label) {
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
static int pst_v35_scenario_size(const std::string& label) {
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
static std::string pst_v35_scenario_label(const std::string& state, int size) {
  return state + "_" + std::to_string(size);
}

/**
 * Drop zero-duration entries and merge adjacent identical labels in a transition map.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static V35TransitionMap pst_v35_compact_transition_map(const V35TransitionMap& edge_map) {
  V35TransitionMap out;
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
static void pst_v35_compact_scenario_edge_maps(
    V35ScenarioEdge& edge,
    double time_tolerance) {
  if (edge.map.labels.size() != edge.path_map.labels.size() ||
      edge.map.durations.size() != edge.path_map.durations.size() ||
      edge.map.labels.size() != edge.map.durations.size()) {
    Rcpp::stop("V35 scenario state/path maps lost segment alignment");
  }
  V35TransitionMap state_out;
  V35TransitionMap path_out;
  std::string final_zero_state;
  std::string final_zero_path;
  bool saw_zero_segment = false;
  // Retain the shared chronological boundary unless both identities match.
  for (int segment_id = 0;
       segment_id < static_cast<int>(edge.map.durations.size());
       ++segment_id) {
    double state_duration = edge.map.durations[segment_id];
    double path_duration = edge.path_map.durations[segment_id];
    if (!R_finite(state_duration) || !R_finite(path_duration) ||
        std::fabs(state_duration - path_duration) > time_tolerance) {
      Rcpp::stop("V35 scenario state/path durations diverged");
    }
    if (state_duration == 0.0) {
      // Retain the final committed identity in case this entire biological
      // edge is a zero-length terminal history. Zero intermediates on an edge
      // with positive exposure remain topology-only and are omitted below.
      final_zero_state = edge.map.labels[segment_id];
      final_zero_path = edge.path_map.labels[segment_id];
      saw_zero_segment = true;
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
  if (state_out.labels.empty() && saw_zero_segment) {
    state_out.labels.push_back(final_zero_state);
    state_out.durations.push_back(0.0);
    path_out.labels.push_back(final_zero_path);
    path_out.durations.push_back(0.0);
  }
  edge.map = state_out;
  edge.path_map = path_out;
}

/** Concatenate maps without compacting boundaries owned by a parallel map. */
static V35TransitionMap pst_v35_concat_maps_raw(
    const std::vector<V35TransitionMap>& maps) {
  V35TransitionMap out;
  for (const V35TransitionMap& map : maps) {
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
static double pst_v35_map_sum(const V35TransitionMap& edge_map) {
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
static V35TransitionMap pst_v35_concat_maps(const std::vector<V35TransitionMap>& maps) {
  V35TransitionMap out;
  // Append source maps in edge-path order before compacting their shared boundaries.
  for (const V35TransitionMap& edge_map : maps) {
    // Copy each label/duration pair together so map identity cannot drift from exposure.
    for (int i = 0; i < static_cast<int>(edge_map.durations.size()); ++i) {
      out.labels.push_back(edge_map.labels[i]);
      out.durations.push_back(edge_map.durations[i]);
    }
  }
  return pst_v35_compact_transition_map(out);
}

/**
 * Find the largest node id before topology compaction allocates replacements.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static int pst_v35_max_scenario_node(const V35ScenarioTree& tree) {
  int max_node = 0;
  // Inspect both endpoints of every scenario edge to bound future node allocation.
  for (const V35ScenarioEdge& edge : tree.edges) {
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
static void pst_v35_renumber_internal_nodes(V35ScenarioTree& tree) {
  int ntips = static_cast<int>(tree.tip_labels.size());
  std::set<int> all_nodes;
  std::vector<int> parents_order;
  std::set<int> parents_seen;
  std::set<int> children;
  // Collect the node domain and first parent order from every retained edge;
  // this order makes root selection deterministic after compaction.
  for (const V35ScenarioEdge& edge : tree.edges) {
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
  for (V35ScenarioEdge& edge : tree.edges) {
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
static void pst_v35_renumber_tree_nodes(V35ScenarioTree& tree) {
  std::vector<int> parents_order;
  std::set<int> parents_seen;
  std::vector<int> children_order;
  std::set<int> children_seen;
  std::set<int> all_nodes_set;
  // Discover first-occurrence parent/child order and the complete node domain
  // after topology edits that may have converted former internals into leaves.
  for (const V35ScenarioEdge& edge : tree.edges) {
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
  for (V35ScenarioEdge& edge : tree.edges) {
    edge.parent = node_map[edge.parent];
    edge.child = node_map[edge.child];
  }
}

/**
 * Convert an R simmap list into the internal C++ scenario-tree representation.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static V35ScenarioTree pst_v35_scenario_tree_from_r(
    Rcpp::List tree,
    double time_tolerance) {
  Rcpp::IntegerMatrix edge = tree["edge"];
  Rcpp::List maps = tree["maps"];
  bool has_path_maps = tree.containsElementNamed("path_maps");
  Rcpp::List path_maps = has_path_maps ?
    Rcpp::as<Rcpp::List>(tree["path_maps"]) : Rcpp::List();
  if (has_path_maps && path_maps.size() != maps.size()) {
    Rcpp::stop("V35 scenario state/path maps are not edge-parallel");
  }
  Rcpp::CharacterVector tip_labels = tree["tip.label"];
  Rcpp::IntegerVector raw_scenario_edge_ids =
    tree.containsElementNamed("raw_scenario_edge_ids") ?
      Rcpp::as<Rcpp::IntegerVector>(tree["raw_scenario_edge_ids"]) :
      Rcpp::IntegerVector();
  if (raw_scenario_edge_ids.size() > 0 &&
      raw_scenario_edge_ids.size() != edge.nrow()) {
    Rcpp::stop("V35 raw scenario-edge ids are not edge-parallel");
  }
  V35ScenarioTree out;
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
      Rcpp::stop("V35 scenario state/path maps are not segment-parallel");
    }
    V35TransitionMap map;
    V35TransitionMap path_map;
    // Copy each named duration as one internal map segment.
    for (int j = 0; j < edge_map.size(); ++j) {
      map.labels.push_back(Rcpp::as<std::string>(names[j]));
      map.durations.push_back(edge_map[j]);
      path_map.labels.push_back(Rcpp::as<std::string>(path_names[j]));
      path_map.durations.push_back(edge_path_map[j]);
    }
    out.edges.push_back(V35ScenarioEdge{
      edge(i, 0),
      edge(i, 1),
      map,
      path_map,
      std::set<int>{
        raw_scenario_edge_ids.size() > 0 ?
          raw_scenario_edge_ids[i] : i + 1
      }
    });
    pst_v35_compact_scenario_edge_maps(
      out.edges.back(),
      time_tolerance
    );
  }
  return out;
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
static Rcpp::List pst_v35_scenario_lineage_size_maps(
    const V35ScenarioTree& tree) {
  Rcpp::List lineage_size_maps(tree.edges.size());
  // Convert every canonical scenario edge without changing edge or segment order.
  for (int edge_id = 0; edge_id < static_cast<int>(tree.edges.size()); ++edge_id) {
    const V35TransitionMap& edge_map = tree.edges[edge_id].map;
    Rcpp::NumericVector sizes(edge_map.labels.size());
    Rcpp::CharacterVector states(edge_map.labels.size());
    // Split each committed scenario label into the state key and displayed lineage size.
    for (int segment_id = 0;
         segment_id < static_cast<int>(edge_map.labels.size());
         ++segment_id) {
      const std::string& label = edge_map.labels[segment_id];
      states[segment_id] = pst_v35_scenario_state(label);
      sizes[segment_id] = pst_v35_scenario_size(label);
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
static Rcpp::List pst_v35_scenario_tree_to_r(const V35ScenarioTree& tree) {
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
    edge_length[i] = pst_v35_map_sum(tree.edges[i].map);
    maps[i] = pst_v35_transition_map_to_vector(tree.edges[i].map);
    path_maps[i] = pst_v35_transition_map_to_vector(tree.edges[i].path_map);
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
      const V35TransitionMap& parent_paths =
        tree.edges[parent_edge->second].path_map;
      const V35TransitionMap& child_paths = tree.edges[edge_id].path_map;
      if (parent_paths.labels.empty() || child_paths.labels.empty()) {
        Rcpp::stop("V35 final scenario edge lacks path identity for classification");
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
    Rcpp::Named("lineage_size_maps") = pst_v35_scenario_lineage_size_maps(tree)
  );
  out["raw_scenario_edge_ids"] = raw_scenario_edge_ids;
  out.attr("pst_v35_raw_scenario_edge_ids") = raw_scenario_edge_ids;
  out.attr("class") = Rcpp::CharacterVector::create("simmap", "phylo");
  return out;
}

/**
 * Remove zero-length scenario-tree edges created by cleanup while preserving mapped durations.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static void pst_v35_contract_zero_length_edges_cpp(
    V35ScenarioTree& tree,
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
      if (pst_v35_map_sum(tree.edges[edge_index].map) > time_tolerance ||
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
      V35TransitionMap bridge_map = tree.edges[edge_index].map;
      // A unary bridge contributes its map prefix to its sole continuation.
      if (child_edges.size() == 1) {
        int child_edge_id = child_edges[0];
        tree.edges[child_edge_id].map = pst_v35_concat_maps_raw(
          std::vector<V35TransitionMap>{bridge_map, tree.edges[child_edge_id].map}
        );
        tree.edges[child_edge_id].path_map = pst_v35_concat_maps_raw(
          std::vector<V35TransitionMap>{
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
  for (V35ScenarioEdge& edge : tree.edges) {
    pst_v35_compact_scenario_edge_maps(edge, time_tolerance);
  }
  pst_v35_renumber_internal_nodes(tree);
}

/**
 * Remove unary internal nodes created by compaction and repair topology metadata.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static void pst_v35_contract_unary_internal_edges_cpp(
    V35ScenarioTree& tree,
    double time_tolerance) {
  int ntips = static_cast<int>(tree.tip_labels.size());
  int max_node = pst_v35_max_scenario_node(tree);
  std::vector<int> parents;
  std::vector<int> children;
  std::set<int> parent_set;
  std::set<int> child_set;
  // Build parent/child sets and degree inputs from the current scenario topology.
  for (const V35ScenarioEdge& edge : tree.edges) {
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
  std::vector<V35ScenarioEdge> new_edges;
  // Start output chains only at non-unary parents; unary-parent edges are consumed by their ancestor chain.
  for (int edge_id = 0; edge_id < static_cast<int>(tree.edges.size()); ++edge_id) {
    int parent = tree.edges[edge_id].parent;
    // Skip an edge whose parent will be absorbed into an earlier concatenated chain.
    if (parent >= 0 && parent < static_cast<int>(is_unary.size()) && is_unary[parent]) {
      continue;
    }
    int child = tree.edges[edge_id].child;
    std::vector<V35TransitionMap> edge_maps;
    std::vector<V35TransitionMap> edge_path_maps;
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
      const V35ScenarioEdge& child_edge = tree.edges[child_edge_id];
      bool duplicate_endpoint_contribution = false;
      if (child_edge.child <= ntips &&
          pst_v35_map_sum(child_edge.map) <= time_tolerance &&
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
          pst_v35_scenario_state(edge_maps.back().labels.back()) ==
            pst_v35_scenario_state(child_edge.map.labels.front()) &&
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
    new_edges.push_back(V35ScenarioEdge{
      parent,
      child,
      pst_v35_concat_maps_raw(edge_maps),
      pst_v35_concat_maps_raw(edge_path_maps),
      raw_edge_ids
    });
    pst_v35_compact_scenario_edge_maps(
      new_edges.back(),
      time_tolerance
    );
  }
  tree.edges = new_edges;
  // Compact concatenated chains and restore valid public node numbering.
  for (V35ScenarioEdge& edge : tree.edges) {
    pst_v35_compact_scenario_edge_maps(edge, time_tolerance);
  }
  pst_v35_renumber_internal_nodes(tree);
}

/**
 * Test exact mapped-segment identity before sibling scenario histories merge.
 *
 * Inputs are read only; outputs preserve canonical R-facing ids, topology, and map alignment.
 */
static bool pst_v35_maps_identical(const V35TransitionMap& left, const V35TransitionMap& right) {
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
static V35TransitionMap pst_v35_multiply_support_map(const V35TransitionMap& edge_map, int multiplicity) {
  // A single represented lineage or an empty trajectory needs no support scaling.
  if (multiplicity <= 1 || edge_map.labels.empty()) {
    return edge_map;
  }
  V35TransitionMap out = edge_map;
  // Multiply each segment's lineage count because all merged tips share the full mapped history.
  for (std::string& label : out.labels) {
    std::string state = pst_v35_scenario_state(label);
    int size = pst_v35_scenario_size(label);
    label = pst_v35_scenario_label(state, size * multiplicity);
  }
  return out;
}

/**
 * Restore public tip-label order after C++ scenario-tree compaction merges tips.
 *
 * Inputs are read only and failures stop before a partially formatted public tree is returned.
 */
static void pst_v35_reorder_tip_labels_cpp(
    V35ScenarioTree& tree,
    const Rcpp::CharacterVector& original_tip_labels) {
  std::map<std::string, int> order;
  // Rank original biological labels once so composite scenario tips can be restored to input order.
  for (int i = 0; i < original_tip_labels.size(); ++i) {
    order[Rcpp::as<std::string>(original_tip_labels[i])] = i;
  }
  // Reorder every possibly composite scenario-tip label without changing its represented-tip set.
  for (std::string& label : tree.tip_labels) {
    std::vector<std::string> tips = pst_v35_split_csv(label);
    std::stable_sort(tips.begin(), tips.end(), [&](const std::string& a, const std::string& b) {
      auto ai = order.find(a);
      auto bi = order.find(b);
      int av = ai == order.end() ? std::numeric_limits<int>::max() : ai->second;
      int bv = bi == order.end() ? std::numeric_limits<int>::max() : bi->second;
      return av < bv;
    });
    label = pst_v35_join_csv(tips);
  }
}

struct V35TerminalEventRow {
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
static Rcpp::DataFrame pst_v35_terminal_events_from_tree(
    const V35ScenarioTree& tree,
    const Rcpp::List& summary,
    double time_tolerance) {
  Rcpp::DataFrame terminal_records = Rcpp::as<Rcpp::DataFrame>(
    summary["terminal_records"]
  );
  Rcpp::CharacterVector input_tip_labels = summary["tip_labels"];
  if (terminal_records.nrows() != input_tip_labels.size()) {
    Rcpp::stop("V35 terminal records are not aligned one-to-one with input tips");
  }
  Rcpp::IntegerVector tip_id = terminal_records["tip_id"];
  Rcpp::CharacterVector tip_label = terminal_records["tip_label"];
  Rcpp::IntegerVector phylo_edge_id = terminal_records["final_phylo_edge_id"];
  if (!terminal_records.containsElementNamed("stable_scenario_id")) {
    Rcpp::stop("V35 terminal records lack stable scenario identity");
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
        "V35 direct scenario edge lacks one exact stable identity"
      );
    }
    int stable_id = *tree.edges[edge_id].raw_edge_ids.begin();
    if (stable_id < 1 ||
        !public_edge_id_by_stable_id
          .insert(std::make_pair(stable_id, edge_id))
          .second) {
      Rcpp::stop("V35 direct scenario tree has duplicate stable identity");
    }
  }
  int max_node = pst_v35_max_scenario_node(tree);
  std::vector<int> edge_by_child(max_node + 1, -1);
  std::vector<std::vector<int> > outgoing(max_node + 1);
  for (int edge_id = 0; edge_id < nedge; ++edge_id) {
    int parent = tree.edges[edge_id].parent;
    int child = tree.edges[edge_id].child;
    if (parent < 1 || child < 1 || parent > max_node || child > max_node ||
        edge_by_child[child] != -1) {
      Rcpp::stop("V35 terminal-event topology is malformed");
    }
    edge_by_child[child] = edge_id;
    outgoing[parent].push_back(edge_id);
  }
  int root = -1;
  for (int node = 1; node <= max_node; ++node) {
    if (!outgoing[node].empty() && edge_by_child[node] == -1) {
      if (root != -1) {
        Rcpp::stop("V35 terminal-event topology has multiple roots");
      }
      root = node;
    }
  }
  if (root < 1) {
    Rcpp::stop("V35 terminal-event topology has no root");
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
      double length = pst_v35_map_sum(tree.edges[edge_id].map);
      if (!std::isfinite(length) || length < 0.0) {
        Rcpp::stop("V35 terminal-event scenario edge has invalid length");
      }
      node_time[child] = node_time[parent] + length;
      queue.push_back(child);
    }
  }

  double time_scale = pst_v35_batch_time_scale(time_tolerance);
  std::vector<long long> edge_start_key(nedge);
  std::vector<long long> edge_end_key(nedge);
  for (int edge_id = 0; edge_id < nedge; ++edge_id) {
    edge_start_key[edge_id] = pst_v35_batch_time_key(
      node_time[tree.edges[edge_id].parent],
      time_scale
    );
    edge_end_key[edge_id] = pst_v35_batch_time_key(
      node_time[tree.edges[edge_id].child],
      time_scale
    );
  }

  std::vector<V35TerminalEventRow> events;
  events.reserve(terminal_records.nrows());
  long long maximum_terminal_key = std::numeric_limits<long long>::min();
  for (int row = 0; row < terminal_records.nrows(); ++row) {
    if (tip_id[row] < 1 || tip_id[row] > input_tip_labels.size() ||
        phylo_edge_id[row] < 1 || stable_scenario_id[row] < 1 ||
        path_id[row] < 1 || !R_finite(terminal_time[row]) ||
        !state_id_by_path.count(path_id[row]) ||
        !state_by_id.count(state_id_by_path[path_id[row]])) {
      Rcpp::stop("V35 terminal record contains invalid endpoint identity");
    }
    long long time_key = pst_v35_batch_time_key(
      terminal_time[row],
      time_scale
    );
    std::map<int, int>::const_iterator public_edge =
      public_edge_id_by_stable_id.find(stable_scenario_id[row]);
    if (public_edge == public_edge_id_by_stable_id.end()) {
      std::ostringstream diagnostic;
      diagnostic
        << "V35 terminal record has no public edge for its stable identity"
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
        << "V35 terminal time lies outside its exact stable scenario edge"
        << " [tip=" << Rcpp::as<std::string>(tip_label[row])
        << ", stable_scenario_id=" << stable_scenario_id[row]
        << ", terminal_time_key=" << time_key
        << ", edge_time_key=" << edge_start_key[canonical_edge_id]
        << ".." << edge_end_key[canonical_edge_id]
        << "]";
      Rcpp::stop(diagnostic.str());
    }
    events.push_back(V35TerminalEventRow{
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
      const V35TerminalEventRow& left,
      const V35TerminalEventRow& right) {
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
  out.attr("pst_v35_terminal_time_keys") = out_terminal_time_key;
  return out;
}

} // namespace

/**
 * Serialize the final traversal-owned scenario registry without changing its
 * topology, stable identities, or segment boundaries.
 */
Rcpp::List pst_v35_scenario_tree_initial_impl(Rcpp::List summary) {
  if (!summary.containsElementNamed("scenario_edge_records") ||
      !summary.containsElementNamed("state_lookup") ||
      !summary.containsElementNamed("path_lookup") ||
      !summary.containsElementNamed("tip_labels")) {
    Rcpp::stop("V35 direct scenario serialization lacks required summary fields");
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
        std::string("V35 direct scenario serialization lacks ") + field
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
    Rcpp::stop("V35 direct scenario records are not nonempty and edge-parallel");
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
          << "V35 direct scenario record has an invalid parent stable ID"
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
        Rcpp::stop("V35 direct scenario registry contains a parent cycle");
      }
      cursor = parent;
    }
  }
  if (root_record_count < 1) {
    Rcpp::stop("V35 direct scenario registry contains no root record");
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
    Rcpp::stop("V35 direct scenario registry contains no terminal topology edge");
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
    std::pair<V35TransitionMap, V35TransitionMap> serialized =
      pst_v35_make_scenario_maps(
        states,
        paths,
        sizes,
        durations,
        state_labels,
        path_labels
      );
    maps[edge_index] =
      pst_v35_scenario_edge_map_to_vector(serialized.first);
    path_maps[edge_index] =
      pst_v35_transition_map_to_vector(serialized.second);
    edge_lengths[edge_index] = pst_v35_map_sum(serialized.first);
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
        "V35 direct scenario leaf lacks traversal-owned tip membership"
      );
    }
    std::ostringstream label;
    for (int member_index = 0;
         member_index < represented_tip_ids.size();
         ++member_index) {
      int tip_id = represented_tip_ids[member_index];
      if (tip_id < 1 || tip_id > input_tip_labels.size()) {
        Rcpp::stop("V35 direct scenario leaf contains an invalid tip ID");
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
 * Finalize direct scenario output by validating stable topology and exposing
 * traversal-owned state, path, support, and edge-type records unchanged.
 */
Rcpp::List pst_v35_canonicalize_scenario_tree_impl(
    Rcpp::List tree,
    Rcpp::List summary,
    double time_tolerance) {
  if (!std::isfinite(time_tolerance) || time_tolerance <= 0.0) {
    Rcpp::stop("V35 scenario-tree time tolerance must be finite and positive");
  }
  if (!tree.containsElementNamed("edge") ||
      !tree.containsElementNamed("edge.length") ||
      !tree.containsElementNamed("maps") ||
      !tree.containsElementNamed("path_maps") ||
      !tree.containsElementNamed("raw_scenario_edge_ids") ||
      !summary.containsElementNamed("scenario_edge_records")) {
    Rcpp::stop("V35 direct scenario finalization lacks required fields");
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
        std::string("V35 direct scenario finalization lacks ") + field
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
    Rcpp::stop("V35 direct scenario output is not edge-parallel");
  }

  std::set<int> seen_stable_ids;
  std::map<int, int> public_edge_by_stable_id;
  std::map<int, int> public_edge_by_child_node;
  for (int edge_index = 0; edge_index < nedge; ++edge_index) {
    int stable_id = stable_ids[edge_index];
    if (stable_id < 1 || stable_id > nedge ||
        !seen_stable_ids.insert(stable_id).second) {
      Rcpp::stop("V35 direct scenario output has invalid stable IDs");
    }
    public_edge_by_stable_id[stable_id] = edge_index;
    if (!public_edge_by_child_node
          .insert(std::make_pair(edge(edge_index, 1), edge_index))
          .second) {
      Rcpp::stop("V35 direct scenario topology has duplicate child nodes");
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
        << "V35 direct scenario parent mismatch [stable_id="
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
      Rcpp::stop("V35 direct scenario record has an invalid edge type");
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
      Rcpp::stop("V35 direct scenario record has misaligned segments");
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
        Rcpp::stop("V35 direct scenario record has an invalid segment time");
      }
      cursor += duration;
      if (std::abs(cursor - end_times[segment_index]) > time_tolerance) {
        std::ostringstream diagnostic;
        diagnostic
          << "V35 direct scenario duration/end-time mismatch [stable_id="
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
          "V35 direct scenario public map changed a committed duration"
        );
      }
      state_map[segment_index] = duration;
      state_names[segment_index] = pst_v35_scenario_state(
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
      Rcpp::stop("V35 direct scenario edge length changed committed duration");
    }
  }

  Rcpp::List out = Rcpp::clone(tree);
  out["maps"] = public_maps;
  out["edge_type"] = Rcpp::clone(edge_types);
  out["size_maps"] = Rcpp::List::create(
    Rcpp::Named("lineage_size_maps") = lineage_size_maps
  );
  out["raw_scenario_edge_ids"] = singleton_provenance;
  out.attr("pst_v35_raw_scenario_edge_ids") = singleton_provenance;
  V35ScenarioTree terminal_event_tree =
    pst_v35_scenario_tree_from_r(tree, time_tolerance);
  out["terminal_events"] = pst_v35_terminal_events_from_tree(
    terminal_event_tree,
    summary,
    time_tolerance
  );
  out.attr("pst_v35_terminal_trajectory_merged_labels") =
    Rcpp::CharacterVector(0);
  out.attr("class") = Rcpp::CharacterVector::create("simmap", "phylo");
  return out;
}

/**
 * Preserve the traversal-owned summary exactly. Stable scenario identities and
 * LDIF are committed together before tree serialization, so output formatting
 * has no summary rewrite authority.
 */
Rcpp::List pst_v35_canonicalize_summary_scenario_ids_impl(
    Rcpp::List summary,
    Rcpp::List scenario_tree) {
  if (!scenario_tree.containsElementNamed("raw_scenario_edge_ids")) {
    Rcpp::stop("V35 direct scenario tree lacks stable-ID provenance");
  }
  return Rcpp::clone(summary);
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
struct V35MatrixMapLayout {
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
static V35MatrixMapLayout pst_v35_matrix_map_layout(
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

  V35MatrixMapLayout layout;
  layout.edge_start_key.resize(nedge);
  layout.edge_end_key.resize(nedge);
  layout.segment_start_key.resize(nedge);
  layout.state_label.resize(nedge);
  if (include_paths) {
    layout.path_label.resize(nedge);
  }

  // Build one segment-start vector per edge. Segment durations are accumulated
  // in double precision and every boundary is quantized with V35's single
  // construction tolerance, matching event batching and final edge ownership.
  for (int edge_id = 0; edge_id < nedge; ++edge_id) {
    double edge_start = node_time[edge(edge_id, 0)];
    double edge_end = node_time[edge(edge_id, 1)];
    layout.edge_start_key[edge_id] = pst_v35_batch_time_key(
      edge_start,
      time_scale
    );
    layout.edge_end_key[edge_id] = pst_v35_batch_time_key(
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
        pst_v35_batch_time_key(segment_start, time_scale)
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
static int pst_v35_matrix_map_step(
    const V35MatrixMapLayout& layout,
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
static int pst_v35_terminal_endpoint_scenario_step(
    const V35MatrixMapLayout& scenario_layout,
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
Rcpp::List pst_v35_scenario_matrix_coordinates_impl(
    Rcpp::List scenario_tree,
    Rcpp::List phylo_tree,
    Rcpp::NumericVector time_vec,
    Rcpp::IntegerMatrix existing_scenario_ids,
    Rcpp::IntegerMatrix existing_phylo_ids,
    Rcpp::IntegerMatrix existing_path_ids,
    Rcpp::DataFrame path_lookup,
    bool include_paths,
    double time_tolerance) {
  if (existing_scenario_ids.nrow() != existing_phylo_ids.nrow() ||
      existing_scenario_ids.ncol() != existing_phylo_ids.ncol() ||
      existing_scenario_ids.nrow() != existing_path_ids.nrow() ||
      existing_scenario_ids.ncol() != existing_path_ids.ncol() ||
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

  if (!scenario_tree.containsElementNamed("edge") ||
      scenario_tree.attr("pst_v35_raw_scenario_edge_ids") == R_NilValue) {
    Rcpp::stop(
      "V35 scenario matrix materialization lacks traversal-owned edge IDs"
    );
  }
  Rcpp::IntegerMatrix scenario_edge = scenario_tree["edge"];
  Rcpp::List raw_scenario_edge_ids = Rcpp::as<Rcpp::List>(
    scenario_tree.attr("pst_v35_raw_scenario_edge_ids")
  );
  int nscenario_edge = scenario_edge.nrow();
  if (scenario_edge.ncol() != 2 ||
      nscenario_edge < 1 ||
      raw_scenario_edge_ids.size() != nscenario_edge) {
    Rcpp::stop(
      "V35 scenario matrix materialization received misaligned edge IDs"
    );
  }

  // Direct serialization preserves one traversal-owned stable ID for every
  // public edge. Invert that edge-parallel mapping once, then translate the
  // already committed matrix IDs without walking topology or reconstructing
  // biological ownership.
  std::vector<int> public_edge_by_stable_id(nscenario_edge + 1, 0);
  bool stable_ids_are_public_rows = true;
  for (int edge_index = 0; edge_index < nscenario_edge; ++edge_index) {
    Rcpp::IntegerVector stable_ids = raw_scenario_edge_ids[edge_index];
    if (stable_ids.size() != 1) {
      Rcpp::stop(
        "V35 scenario matrix materialization requires singleton stable IDs"
      );
    }
    int stable_id = stable_ids[0];
    if (stable_id < 1 ||
        stable_id > nscenario_edge ||
        public_edge_by_stable_id[stable_id] != 0) {
      Rcpp::stop(
        "V35 scenario matrix materialization received invalid stable IDs"
      );
    }
    public_edge_by_stable_id[stable_id] = edge_index + 1;
    stable_ids_are_public_rows =
      stable_ids_are_public_rows && stable_id == edge_index + 1;
  }

  Rcpp::IntegerMatrix scenario_edge_ids = stable_ids_are_public_rows ?
    existing_scenario_ids :
    Rcpp::clone(existing_scenario_ids);
  if (!stable_ids_are_public_rows) {
    for (R_xlen_t cell = 0; cell < scenario_edge_ids.size(); ++cell) {
      int stable_id = scenario_edge_ids[cell];
      if (stable_id == NA_INTEGER) {
        continue;
      }
      if (stable_id < 1 ||
          stable_id > nscenario_edge ||
          public_edge_by_stable_id[stable_id] == 0) {
        Rcpp::stop(
          "V35 scenario matrix contains an unknown traversal-owned edge ID"
        );
      }
      scenario_edge_ids[cell] = public_edge_by_stable_id[stable_id];
    }
  }

  double time_scale = pst_v35_batch_time_scale(time_tolerance);
  V35MatrixMapLayout scenario_layout = pst_v35_matrix_map_layout(
    scenario_tree,
    time_scale,
    include_paths,
    "Scenario tree"
  );
  V35MatrixMapLayout phylo_layout = pst_v35_matrix_map_layout(
    phylo_tree,
    time_scale,
    false,
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
  std::vector<SEXP> path_label_by_id;
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
    Rcpp::IntegerVector lookup_path_id = path_lookup["path_id"];
    Rcpp::CharacterVector lookup_label = path_lookup["label"];
    int maximum_path_id = 0;
    for (int path_id : lookup_path_id) {
      maximum_path_id = std::max(maximum_path_id, path_id);
    }
    path_label_by_id.assign(maximum_path_id + 1, NA_STRING);
    for (int index = 0; index < lookup_path_id.size(); ++index) {
      int path_id = lookup_path_id[index];
      if (path_id < 1 ||
          path_id >= static_cast<int>(path_label_by_id.size()) ||
          path_label_by_id[path_id] != NA_STRING) {
        Rcpp::stop("V35 matrix path lookup contains invalid path IDs");
      }
      path_label_by_id[path_id] = lookup_label[index];
    }
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
    time_key[column] = pst_v35_batch_time_key(time_vec[column], time_scale);
  }

  std::vector<int> last_scenario_edge(existing_scenario_ids.nrow(), 0);
  std::vector<int> last_scenario_step(existing_scenario_ids.nrow(), 0);
  std::vector<int> last_phylo_edge(existing_scenario_ids.nrow(), 0);
  std::vector<int> last_phylo_step(existing_scenario_ids.nrow(), 0);
  const auto step_from_monotone_cursor = [](
      const V35MatrixMapLayout& layout,
      int edge_id,
      long long key,
      int& last_edge_id,
      int& last_step_id,
      const std::string& tree_name) {
    if (edge_id != last_edge_id || last_step_id < 1) {
      last_edge_id = edge_id;
      last_step_id = pst_v35_matrix_map_step(
        layout,
        edge_id,
        key,
        tree_name
      );
      return last_step_id;
    }
    int edge_index = edge_id - 1;
    const std::vector<long long>& starts =
      layout.segment_start_key[edge_index];
    if (key < layout.edge_start_key[edge_index]) {
      Rcpp::stop(tree_name + " matrix time precedes its addressed edge");
    }
    if (key >= layout.edge_end_key[edge_index]) {
      last_step_id = static_cast<int>(starts.size());
      return last_step_id;
    }
    while (last_step_id < static_cast<int>(starts.size()) &&
           starts[last_step_id] <= key) {
      ++last_step_id;
    }
    return last_step_id;
  };

  // Resolve both tree coordinates for one rectangular cell at a time. Missing
  // edge ids remain missing step/path cells; mismatched coverage is rejected so
  // a saved matrix never exposes only half of a biological correspondence.
  for (int column = 0; column < existing_scenario_ids.ncol(); ++column) {
    for (int row = 0; row < existing_scenario_ids.nrow(); ++row) {
      int scenario_edge_id = scenario_edge_ids(row, column);
      int phylo_edge_id = existing_phylo_ids(row, column);
      int path_id = existing_path_ids(row, column);
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
        if (path_id != NA_INTEGER) {
          Rcpp::stop("V35 missing matrix edge retains a path ID");
        }
        continue;
      }
      if (path_id < 1 ||
          (include_paths &&
           (path_id >= static_cast<int>(path_label_by_id.size()) ||
            path_label_by_id[path_id] == NA_STRING))) {
        Rcpp::stop("V35 scenario matrix contains an unknown path ID");
      }

      int phylo_step_id = step_from_monotone_cursor(
        phylo_layout,
        phylo_edge_id,
        time_key[column],
        last_phylo_edge[row],
        last_phylo_step[row],
        "Phylogeny"
      );
      int scenario_step_id = step_from_monotone_cursor(
        scenario_layout,
        scenario_edge_id,
        time_key[column],
        last_scenario_edge[row],
        last_scenario_step[row],
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
          endpoint_path = Rcpp::as<std::string>(path_label_by_id[path_id]);
        }
        scenario_step_id = pst_v35_terminal_endpoint_scenario_step(
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
        scenario_paths(row, column) = path_label_by_id[path_id];
        phylo_paths(row, column) = path_label_by_id[path_id];
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
static Rcpp::List pst_v35_transition_tree_from_compact_support(
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
    Rcpp::stop("V35 compact transition support fields are not aligned");
  }
  // Terminal path ids and represented-tip groups are one-to-one.
  if (terminal_path_id.size() != terminal_tip_ids.size()) {
    Rcpp::stop("V35 compact transition terminal groups are not aligned");
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
      Rcpp::stop("V35 compact transition support contains an invalid path");
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
      Rcpp::stop("V35 terminal history references an uncommitted path");
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
        Rcpp::stop("V35 transition parent path is outside the registry");
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
  std::vector<V35TransitionMap> ordered_maps;
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
      Rcpp::stop("V35 compact transition edge has invalid state or support");
    }
    old_edges.push_back(std::make_pair(parent_old_id, path_id));
    ordered_maps.push_back(pst_v35_make_transition_map(
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
      // its tips coexist with longer child histories. Emit that placeholder
      // directly in the committed final biological state with zero duration.
      int terminal_old_id = ++max_old_id;
      old_edges.push_back(std::make_pair(path_id, terminal_old_id));
      ordered_maps.push_back(pst_v35_make_transition_map(
        std::vector<int>(1, state_by_path[path_id]),
        std::vector<int>(1, size_by_path[path_id]),
        std::vector<double>(1, 0.0),
        state_labels
      ));
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
    maps[edge_id] = pst_v35_transition_map_to_vector(ordered_maps[edge_id]);
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
        Rcpp::stop("V35 compact transition label is outside mapped.edge");
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
      Rcpp::stop("V35 compact transition leaf has no represented tips");
    }
    std::string label = pst_v35_join_tip_labels(
      represented->second,
      tip_labels
    );
    if (represented_tip_ids_by_label.count(label)) {
      Rcpp::stop("V35 compact transition tips have duplicate public labels");
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
  V35ScenarioTree compact_tree = pst_v35_scenario_tree_from_r(
    tree,
    time_tolerance
  );
  pst_v35_contract_unary_internal_edges_cpp(compact_tree, time_tolerance);
  Rcpp::List compact_output = pst_v35_scenario_tree_to_r(compact_tree);
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
        Rcpp::stop("V35 compact transition map lost its mapped-edge column");
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
      Rcpp::stop("V35 compact transition tip lost represented biological ids");
    }
    compact_represented_tip_ids[tip_id] = Rcpp::wrap(represented->second);
  }
  compact_output.attr("pst_v35_represented_tip_ids") =
    compact_represented_tip_ids;
  return compact_output;
}

/**
 * Build the public transition tree from compact transition support.
 *
 * This is an Rcpp boundary: inputs are R objects and outputs must preserve documented V35 public shapes.
 */
Rcpp::List pst_v35_transition_tree_from_support_impl(
    Rcpp::List transition_support,
    Rcpp::CharacterVector tip_labels,
    double time_tolerance) {
  // V35 accepts only traversal-owned compact path records. Dense history
  // matrices require post-traversal regrouping and are an ownership failure.
  if (!transition_support.containsElementNamed("edge_path_id")) {
    Rcpp::stop("V35 transition support omitted compact path records");
  }
  SEXP support_digest = transition_support.attr(
    "pst_v35_transition_support_digest"
  );
  if (Rf_isNull(support_digest)) {
    Rcpp::stop("V35 transition support omitted its ownership digest");
  }
  Rcpp::List out = pst_v35_transition_tree_from_compact_support(
    transition_support,
    tip_labels,
    time_tolerance
  );
  out.attr("pst_v35_transition_support_digest") = support_digest;
  return out;

}

/**
 * Strip V35 size suffixes while preserving public simmap topology and exposure.
 *
 * `tree` is cloned before mutation. Scenario-size source labels are reduced to
 * their biological states, while already-state-only source labels are copied
 * exactly. Adjacent equal states are merged in edge-local time order, and
 * `mapped.edge` is rebuilt from the resulting maps. This formatter owns V35
 * label generation at the C++ output boundary; malformed map/topology shapes
 * fail before any partial object is returned to R.
 */
Rcpp::List pst_v35_state_only_tree_impl(
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
    Rcpp::stop("V35 state-only formatter received misaligned path maps");
  }

  // A simmap tree must keep one map per public edge; otherwise state exposure
  // cannot be aligned to topology rows without guessing.
  if (maps.size() != nedge) {
    Rcpp::stop("V35 state-only tree formatter received misaligned edge maps");
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
        Rcpp::stop("V35 state-only formatter found nonparallel path segments");
      }
    }

    // Every duration must have a label because the state-only formatter cannot
    // infer biological state from a bare numeric segment.
    if (edge_labels.size() != edge_map.size()) {
      Rcpp::stop("V35 state-only tree formatter received an unnamed map segment");
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
        Rcpp::stop("V35 state-only tree formatter received an NA map label");
      }
      std::string source_label = Rcpp::as<std::string>(edge_labels[segment_id]);
      std::string state_label = source_labels_are_state_only ?
        source_label : pst_v35_scenario_state(source_label);
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
        Rcpp::stop("V35 state-only tree formatter lost a mapped state column");
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
      Rcpp::stop("V35 state-only summary capture requires state levels");
    }
    Rcpp::CharacterVector summary_states(state_levels);
    tree.attr("pst_v35_tree_summary_parts") =
      pst_v35_capture_tree_summary_parts_cpp(
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
 * This is an Rcpp boundary: inputs are R objects and outputs must preserve documented V35 public shapes.
 */
Rcpp::List pst_v35_pathify_state_tree_impl(
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
  path_labels.push_back(pst_v35_path_label(root_components, sep));
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
  std::vector<int> ordered = pst_v35_preorder_edges(edge, root);
  // Traverse edges root-first so every parent node's path is known before its daughters.
  for (int edge_id : ordered) {
    int row = edge_id - 1;
    int parent = edge(row, 0);
    int child = edge(row, 1);
    int current_path_id = path_id_at_node[parent];
    int current_tail_state_id = tail_state_id_at_node[parent];
    // Missing parent path ownership indicates broken topology or path propagation and cannot be formatted safely.
    if (current_path_id == 0 || current_path_id == NA_INTEGER) {
      Rcpp::stop("Missing V35 path id at parent node %d", parent);
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
          path_labels.push_back(pst_v35_path_label(next_components, sep));
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
Rcpp::List pst_v35_transition_size_maps_impl(
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
    Rcpp::stop("V35 transition size-map trees are not edge-parallel");
  }
  // Verify every topology endpoint before segment-level values are allocated;
  // a reordered or changed path projection cannot safely share map positions.
  for (int edge_id = 0; edge_id < nedge; ++edge_id) {
    // Different endpoints prove the two projections no longer identify the
    // same public transition-tree edge.
    if (state_edge(edge_id, 0) != path_edge(edge_id, 0) ||
        state_edge(edge_id, 1) != path_edge(edge_id, 1)) {
      Rcpp::stop("V35 transition size-map trees have different topology");
    }
  }

  {
    if (cumulative_lineage_trans.nrow() < 1 ||
        cumulative_scenario_trans.nrow() < 1) {
      Rcpp::stop("V35 cumulative size-map matrices have no committed rows");
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
        Rcpp::stop("V35 cumulative lineage trans contains an unnamed column");
      }
      std::string label = Rcpp::as<std::string>(lineage_names[column]);
      if (label == "total" || label == "summary_time_vec") {
        continue;
      }
      if (lineage_by_path.count(label)) {
        Rcpp::stop("V35 cumulative lineage trans contains duplicate paths");
      }
      lineage_by_path[label] = cumulative_lineage_trans(
        lineage_final_row,
        column
      );
    }
    // Index the final cumulative value of every canonical scenario path.
    for (int column = 0; column < cumulative_scenario_trans.ncol(); ++column) {
      if (scenario_names[column] == NA_STRING) {
        Rcpp::stop("V35 cumulative scenario trans contains an unnamed column");
      }
      std::string label = Rcpp::as<std::string>(scenario_names[column]);
      if (label == "total" || label == "summary_time_vec") {
        continue;
      }
      if (scenario_by_path.count(label)) {
        Rcpp::stop("V35 cumulative scenario trans contains duplicate paths");
      }
      scenario_by_path[label] = cumulative_scenario_trans(
        scenario_final_row,
        column
      );
    }
    if (lineage_by_path.size() != scenario_by_path.size()) {
      Rcpp::stop("V35 cumulative lineage/scenario trans paths are not parallel");
    }
    for (std::map<std::string, double>::const_iterator path =
           lineage_by_path.begin();
         path != lineage_by_path.end();
         ++path) {
      if (!scenario_by_path.count(path->first)) {
        Rcpp::stop("V35 cumulative lineage/scenario trans paths are unmatched");
      }
    }

    if (path_labels_by_id.size() < 1) {
      Rcpp::stop("V35 terminal size support has no canonical path labels");
    }
    std::map<std::string, double> terminal_lineage_by_path;
    std::map<std::string, double> terminal_scenario_by_path;
    auto index_terminal_counts = [&] (
        const Rcpp::IntegerVector& counts,
        const std::string& family,
        std::map<std::string, double>& values_by_path) {
      Rcpp::CharacterVector count_names = counts.names();
      if (counts.size() != count_names.size()) {
        Rcpp::stop("V35 terminal " + family + " counts lack path ids");
      }
      // Resolve each sparse positive terminal contribution through the shared
      // C++ path registry. No public tree or active path-wide TT row is scanned.
      for (int index = 0; index < counts.size(); ++index) {
        if (count_names[index] == NA_STRING || counts[index] < 0) {
          Rcpp::stop("V35 terminal " + family + " count is invalid");
        }
        std::string path_id_text = Rcpp::as<std::string>(count_names[index]);
        char* end = nullptr;
        long path_id = std::strtol(path_id_text.c_str(), &end, 10);
        if (end == path_id_text.c_str() || *end != '\0' ||
            path_id < 1 || path_id > path_labels_by_id.size()) {
          Rcpp::stop("V35 terminal " + family + " count has an unknown path id");
        }
        Rcpp::String label_value = path_labels_by_id[path_id - 1];
        if (label_value == NA_STRING) {
          Rcpp::stop("V35 terminal " + family + " path has no label");
        }
        std::string label = static_cast<std::string>(label_value);
        if (values_by_path.count(label)) {
          Rcpp::stop("V35 terminal " + family + " counts duplicate a path");
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
      "pst_v35_represented_tip_ids"
    );
    if (represented_tip_ids.size() != ntip ||
        terminal_trajectory_group_id_by_tip.size() < 1) {
      Rcpp::stop("V35 terminal transition support is not tip-aligned");
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
        Rcpp::stop("V35 transition cumulative size-map segments are not parallel");
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
          Rcpp::stop("V35 transition cumulative size map has a missing path");
        }
        std::string path_label =
          Rcpp::as<std::string>(path_labels[segment_id]);
        std::map<std::string, double>::const_iterator lineage_value =
          lineage_by_path.find(path_label);
        std::map<std::string, double>::const_iterator scenario_value =
          scenario_by_path.find(path_label);
        if (lineage_value == lineage_by_path.end() ||
            scenario_value == scenario_by_path.end()) {
          Rcpp::stop("V35 transition cumulative size map has an unmatched path");
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
            Rcpp::stop("V35 terminal transition edge has an invalid tip node");
          }
          Rcpp::IntegerVector represented = represented_tip_ids[
            terminal_tip_node - 1
          ];
          if (represented.size() < 1) {
            Rcpp::stop("V35 terminal transition edge represents no phylogeny tips");
          }
          std::set<int> terminal_groups;
          for (int represented_id : represented) {
            if (represented_id < 1 ||
                represented_id > terminal_trajectory_group_id_by_tip.size()) {
              Rcpp::stop("V35 represented terminal tip id is out of range");
            }
            int group_id = terminal_trajectory_group_id_by_tip[
              represented_id - 1
            ];
            if (group_id < 1) {
              Rcpp::stop("V35 represented terminal tip has no trajectory group");
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
