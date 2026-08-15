#include "pst_v34_debug.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {

/** One retained edge and its clipped map prefix before node renumbering. */
struct V34DebugEdge {
  int original_edge_index;
  int parent_old_id;
  int child_old_id;
  bool child_is_provisional;
  int provisional_id;
  double length;
  Rcpp::NumericVector map;
  int retained_segments;
};

/**
 * Clip a named duration map to one retained edge length.
 *
 * At least one zero-duration state is retained for a frontier edge clipped at
 * its parent. This keeps simmap labels available for root-initialized snapshots.
 */
Rcpp::NumericVector v34_debug_clip_duration_map(
    const Rcpp::NumericVector& source,
    double retained_length,
    double time_tolerance) {
  Rcpp::CharacterVector source_names = source.names();
  std::vector<double> values;
  std::vector<std::string> labels;
  double remaining = std::max(retained_length, 0.0);

  // Consume source segments in biological order until the retained edge length
  // is represented exactly; each iteration appends at most one clipped segment.
  for (int segment_id = 0; segment_id < source.size(); ++segment_id) {
    // Once positive retained duration is exhausted, later biological states
    // lie beyond the snapshot and must not enter the clipped map.
    if (remaining <= time_tolerance && !values.empty()) {
      break;
    }
    double source_duration = source[segment_id];
    double kept_duration = std::min(source_duration, remaining);
    // A zero-length frontier still needs its first state label. Later
    // zero-length segments are omitted because they are not publicly present.
    if (kept_duration > time_tolerance || values.empty()) {
      values.push_back(std::max(kept_duration, 0.0));
      labels.push_back(Rcpp::as<std::string>(source_names[segment_id]));
    }
    remaining -= kept_duration;
  }

  // Source maps are required for simmap clipping; an empty map cannot identify
  // the active frontier state.
  if (values.empty()) {
    Rcpp::stop("V34 debug clipping received an empty edge map");
  }
  // A retained length beyond the source map indicates inconsistent tree/map
  // lengths and must fail instead of silently stretching the last state.
  if (remaining > time_tolerance) {
    Rcpp::stop("V34 debug cutoff exceeds an edge's mapped duration");
  }

  Rcpp::NumericVector out = Rcpp::wrap(values);
  out.names() = Rcpp::wrap(labels);
  return out;
}

/**
 * Clip edge-parallel size-map families to retained map segment counts.
 *
 * Numeric size values are support annotations, not durations; clipping keeps
 * the same leading number of values as the retained duration map.
 */
Rcpp::List v34_debug_clip_size_map_container(
    const Rcpp::List& container,
    const std::vector<V34DebugEdge>& edges) {
  Rcpp::List out = Rcpp::clone(container);
  Rcpp::CharacterVector family_names = container.names();

  // Visit every named size-map family, preserving its public family name and
  // constructing one edge map per retained clipped edge.
  for (int family_id = 0; family_id < container.size(); ++family_id) {
    Rcpp::List source_family = container[family_id];
    Rcpp::List clipped_family(edges.size());
    // Copy the retained prefix from the corresponding original edge map.
    for (int edge_id = 0; edge_id < static_cast<int>(edges.size()); ++edge_id) {
      const V34DebugEdge& edge = edges[edge_id];
      Rcpp::NumericVector source_map =
        source_family[edge.original_edge_index];
      int keep = std::min(
        edge.retained_segments,
        static_cast<int>(source_map.size())
      );
      // Every retained duration segment must have a parallel size annotation.
      if (keep < edge.retained_segments) {
        Rcpp::stop("V34 debug size map is shorter than its duration map");
      }
      Rcpp::NumericVector clipped_map = source_map[Rcpp::Range(0, keep - 1)];
      Rcpp::CharacterVector source_names = source_map.names();
      // Rcpp numeric range slicing does not preserve names automatically; copy
      // the parallel state/path labels required by summary materialization.
      if (source_names.size() == source_map.size()) {
        clipped_map.names() = source_names[Rcpp::Range(0, keep - 1)];
      } else {
        // A size map without parallel labels cannot be assigned biologically.
        Rcpp::stop("V34 debug size map is missing state/path names");
      }
      clipped_family[edge_id] = clipped_map;
    }
    out[family_id] = clipped_family;
  }
  out.attr("names") = family_names;
  return out;
}

}  // namespace

/**
 * Clip a finalized or committed-prefix simmap tree at one debug snapshot time.
 *
 * The input tree is read only. Returned topology contains all exposure through
 * `cutoff_time`; branches crossing the cutoff end in deterministic provisional
 * frontier tips, and attached size maps are clipped segment-parallel with maps.
 * Node ids are rebuilt under the `phylo` contract. Invalid time, missing tree
 * fields, inconsistent map labels, or an unreachable root raises an R error
 * without mutating the source tree.
 */
Rcpp::List pst_v34_debug_clip_tree_cpp(
    const Rcpp::List& tree,
    double cutoff_time,
    double time_tolerance) {
  // Snapshot time is root-relative and cannot precede root initialization.
  if (!std::isfinite(cutoff_time) || cutoff_time < -time_tolerance) {
    Rcpp::stop("V34 debug cutoff must be a finite nonnegative time");
  }
  // Public tree topology and map fields are mandatory for a clipped trajectory.
  if (!tree.containsElementNamed("edge") ||
      !tree.containsElementNamed("edge.length") ||
      !tree.containsElementNamed("maps") ||
      !tree.containsElementNamed("tip.label")) {
    Rcpp::stop("V34 debug clipping requires a mapped rooted tree");
  }

  Rcpp::IntegerMatrix edge = tree["edge"];
  Rcpp::NumericVector edge_length = tree["edge.length"];
  Rcpp::List maps = tree["maps"];
  Rcpp::CharacterVector tip_labels = tree["tip.label"];
  int original_tip_count = tip_labels.size();
  int original_edge_count = edge.nrow();
  // Edge rows, lengths, and maps are one parallel topology contract.
  if (edge.ncol() != 2 ||
      edge_length.size() != original_edge_count ||
      maps.size() != original_edge_count) {
    Rcpp::stop("V34 debug tree edge fields are not parallel");
  }

  std::set<int> parent_nodes;
  std::set<int> child_nodes;
  int maximum_node_id = 0;
  // Index original topology and determine the root as the sole parent that is
  // never a child. Every iteration records one edge's endpoint ownership.
  for (int edge_id = 0; edge_id < original_edge_count; ++edge_id) {
    int parent = edge(edge_id, 0);
    int child = edge(edge_id, 1);
    parent_nodes.insert(parent);
    child_nodes.insert(child);
    maximum_node_id = std::max(maximum_node_id, std::max(parent, child));
  }
  std::vector<int> roots;
  // Inspect parent ids for the unique root candidate.
  for (int parent : parent_nodes) {
    // A parent absent from child ids has no incoming edge and is the root.
    if (child_nodes.find(parent) == child_nodes.end()) {
      roots.push_back(parent);
    }
  }
  // A rooted phylogeny must expose exactly one root.
  if (roots.size() != 1U) {
    Rcpp::stop("V34 debug clipping requires exactly one rooted component");
  }
  int root_old_id = roots[0];

  std::vector<std::vector<int> > outgoing(maximum_node_id + 1);
  std::vector<double> node_depth(maximum_node_id + 1, R_NaReal);
  // Build parent-to-edge adjacency used for one forward depth traversal.
  for (int edge_id = 0; edge_id < original_edge_count; ++edge_id) {
    outgoing[edge(edge_id, 0)].push_back(edge_id);
  }
  node_depth[root_old_id] = 0.0;
  std::vector<int> node_queue(1, root_old_id);
  // Traverse reached topology from root to descendants, assigning each child
  // its absolute root-relative endpoint time.
  for (int queue_id = 0; queue_id < static_cast<int>(node_queue.size()); ++queue_id) {
    int parent = node_queue[queue_id];
    // Advance every outgoing edge and enqueue internal children for expansion.
    for (int edge_id : outgoing[parent]) {
      int child = edge(edge_id, 1);
      node_depth[child] = node_depth[parent] + edge_length[edge_id];
      // Only internal children own outgoing edges that need later traversal.
      if (child > original_tip_count) {
        node_queue.push_back(child);
      }
    }
  }

  double tree_height = 0.0;
  // Find the latest terminal endpoint to recognize an exact final snapshot.
  for (int tip_id = 1; tip_id <= original_tip_count; ++tip_id) {
    tree_height = std::max(tree_height, node_depth[tip_id]);
  }
  // The final snapshot must preserve the exact public tree, including all
  // attributes and edge ordering, rather than rebuilding equivalent topology.
  if (cutoff_time >= tree_height - time_tolerance) {
    return Rcpp::clone(tree);
  }

  std::vector<V34DebugEdge> retained_edges;
  int next_provisional_id = 1;
  // Examine every original edge in public order. Reached edges are retained in
  // that order; crossing edges terminate at deterministic provisional tips.
  for (int edge_id = 0; edge_id < original_edge_count; ++edge_id) {
    int parent = edge(edge_id, 0);
    int child = edge(edge_id, 1);
    double start = node_depth[parent];
    double end = node_depth[child];
    // Edges whose parent lies after the cutoff are biologically unreached.
    if (start > cutoff_time + time_tolerance) {
      continue;
    }
    bool fully_reached = end <= cutoff_time + time_tolerance;
    double retained_length = fully_reached ?
      edge_length[edge_id] : std::max(cutoff_time - start, 0.0);
    Rcpp::NumericVector clipped_map = v34_debug_clip_duration_map(
      Rcpp::as<Rcpp::NumericVector>(maps[edge_id]),
      retained_length,
      time_tolerance
    );
    retained_edges.push_back(V34DebugEdge{
      edge_id,
      parent,
      child,
      !fully_reached,
      fully_reached ? 0 : next_provisional_id++,
      retained_length,
      clipped_map,
      static_cast<int>(clipped_map.size())
    });
  }
  // Root snapshots require at least one zero-length frontier edge.
  if (retained_edges.empty()) {
    Rcpp::stop("V34 debug cutoff produced no root frontier edges");
  }

  std::vector<std::string> clipped_tip_labels;
  std::map<int, int> old_tip_to_new;
  std::map<int, int> provisional_to_new;
  // Assign leaf ids in retained edge order for deterministic snapshots.
  for (const V34DebugEdge& retained : retained_edges) {
    // Crossing edges create unique provisional frontier tips.
    if (retained.child_is_provisional) {
      int new_tip_id = static_cast<int>(clipped_tip_labels.size()) + 1;
      provisional_to_new[retained.provisional_id] = new_tip_id;
      clipped_tip_labels.push_back(
        ".pst_v34_frontier_" +
        std::to_string(retained.original_edge_index + 1)
      );
    } else if (retained.child_old_id <= original_tip_count &&
               old_tip_to_new.find(retained.child_old_id) == old_tip_to_new.end()) {
      // Fully reached biological tips preserve their original labels.
      int new_tip_id = static_cast<int>(clipped_tip_labels.size()) + 1;
      old_tip_to_new[retained.child_old_id] = new_tip_id;
      clipped_tip_labels.push_back(
        Rcpp::as<std::string>(tip_labels[retained.child_old_id - 1])
      );
    }
  }

  std::set<int> retained_internal_nodes;
  // Every retained edge parent is an internal node in the clipped topology;
  // fully reached internal children are also retained for outgoing edges.
  for (const V34DebugEdge& retained : retained_edges) {
    retained_internal_nodes.insert(retained.parent_old_id);
    // Internal child endpoints reached by the cutoff remain topology nodes.
    if (!retained.child_is_provisional &&
        retained.child_old_id > original_tip_count) {
      retained_internal_nodes.insert(retained.child_old_id);
    }
  }
  std::map<int, int> old_internal_to_new;
  int next_internal_id = static_cast<int>(clipped_tip_labels.size()) + 1;
  // Renumber internal nodes after all tips, as required by `phylo`.
  for (int old_node_id : retained_internal_nodes) {
    old_internal_to_new[old_node_id] = next_internal_id++;
  }

  Rcpp::IntegerMatrix clipped_edge(retained_edges.size(), 2);
  Rcpp::NumericVector clipped_edge_length(retained_edges.size());
  Rcpp::List clipped_maps(retained_edges.size());
  // Materialize each retained edge under the new tip/internal id spaces.
  for (int edge_id = 0; edge_id < static_cast<int>(retained_edges.size()); ++edge_id) {
    const V34DebugEdge& retained = retained_edges[edge_id];
    clipped_edge(edge_id, 0) = old_internal_to_new[retained.parent_old_id];
    // Provisional, biological-tip, and internal children use distinct id maps.
    if (retained.child_is_provisional) {
      clipped_edge(edge_id, 1) = provisional_to_new[retained.provisional_id];
    } else if (retained.child_old_id <= original_tip_count) {
      // A reached biological tip uses the compact id assigned from its original label.
      clipped_edge(edge_id, 1) = old_tip_to_new[retained.child_old_id];
    } else {
      // A reached internal child uses the separately compacted internal-node id domain.
      clipped_edge(edge_id, 1) = old_internal_to_new[retained.child_old_id];
    }
    clipped_edge_length[edge_id] = retained.length;
    clipped_maps[edge_id] = retained.map;
  }
  clipped_edge.attr("dimnames") = Rcpp::List::create(
    R_NilValue,
    Rcpp::CharacterVector::create("parent", "child")
  );

  std::set<std::string> mapped_states;
  // Collect every retained map label into deterministic lexical column order.
  for (const V34DebugEdge& retained : retained_edges) {
    Rcpp::CharacterVector labels = retained.map.names();
    // Each retained map segment contributes its state/path column label.
    for (int segment_id = 0; segment_id < labels.size(); ++segment_id) {
      mapped_states.insert(Rcpp::as<std::string>(labels[segment_id]));
    }
  }
  std::vector<std::string> state_names(mapped_states.begin(), mapped_states.end());
  std::map<std::string, int> state_column;
  // Index mapped-edge columns by retained label.
  for (int state_id = 0; state_id < static_cast<int>(state_names.size()); ++state_id) {
    state_column[state_names[state_id]] = state_id;
  }
  Rcpp::NumericMatrix mapped_edge(retained_edges.size(), state_names.size());
  // Sum clipped map durations into the standard mapped-edge matrix.
  for (int edge_id = 0; edge_id < static_cast<int>(retained_edges.size()); ++edge_id) {
    Rcpp::NumericVector edge_map = retained_edges[edge_id].map;
    Rcpp::CharacterVector labels = edge_map.names();
    // Each segment adds duration to its state/path column on this edge.
    for (int segment_id = 0; segment_id < edge_map.size(); ++segment_id) {
      std::string label = Rcpp::as<std::string>(labels[segment_id]);
      mapped_edge(edge_id, state_column[label]) += edge_map[segment_id];
    }
  }
  mapped_edge.attr("dimnames") = Rcpp::List::create(
    R_NilValue,
    Rcpp::wrap(state_names)
  );

  Rcpp::List out = Rcpp::clone(tree);
  out["edge"] = clipped_edge;
  out["edge.length"] = clipped_edge_length;
  out["maps"] = clipped_maps;
  out["mapped.edge"] = mapped_edge;
  out["tip.label"] = Rcpp::wrap(clipped_tip_labels);
  out["Nnode"] = static_cast<int>(retained_internal_nodes.size());
  // Size-map containers are edge-parallel and retain the same clipped segment
  // prefixes as duration maps.
  if (tree.containsElementNamed("size_maps") && !Rf_isNull(tree["size_maps"])) {
    out["size_maps"] = v34_debug_clip_size_map_container(
      Rcpp::as<Rcpp::List>(tree["size_maps"]),
      retained_edges
    );
  }
  // Transition path-size maps follow the same edge/segment clipping contract.
  if (tree.containsElementNamed("path_size_maps") &&
      !Rf_isNull(tree["path_size_maps"])) {
    out["path_size_maps"] = v34_debug_clip_size_map_container(
      Rcpp::as<Rcpp::List>(tree["path_size_maps"]),
      retained_edges
    );
  }
  return out;
}
