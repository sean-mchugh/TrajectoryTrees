#include "pst_v35_parse.hpp"
#include <map>
#include <string>

/**
 * Parse the compact R list produced by `pst_v35_encode_input()` into a typed C++ payload without changing id bases.
 *
 * Input: the compact R payload emitted by `pst_v35_encode_input()`.
 * Output: an owning `V35Input` whose vectors retain R's 1-based biological
 * ids. The parser validates vector lengths before traversal can index them and
 * throws an R error without returning a partially initialized input.
 */
V35Input pst_v35_parse_input(const Rcpp::List& payload) {
  // Pull scalar and vector fields from the R boundary list without mutation.
  V35Input input;
  input.ntips = Rcpp::as<int>(payload["ntips"]);
  input.nedge = Rcpp::as<int>(payload["nedge"]);
  input.edge_parent = Rcpp::as<Rcpp::IntegerVector>(payload["edge_parent"]);
  input.edge_child = Rcpp::as<Rcpp::IntegerVector>(payload["edge_child"]);
  input.edge_length = Rcpp::as<Rcpp::NumericVector>(payload["edge_length"]);
  input.map_edge_id = Rcpp::as<Rcpp::IntegerVector>(payload["map_edge_id"]);
  input.map_state_id = Rcpp::as<Rcpp::IntegerVector>(payload["map_state_id"]);
  input.map_duration = Rcpp::as<Rcpp::NumericVector>(payload["map_duration"]);
  input.map_edge_offset = Rcpp::as<Rcpp::IntegerVector>(payload["map_edge_offset"]);
  input.tip_label = Rcpp::as<Rcpp::CharacterVector>(payload["tip_label"]);
  input.state_lookup = Rcpp::as<Rcpp::DataFrame>(payload["state_lookup"]);
  input.root_anchor_state_ids = Rcpp::as<Rcpp::IntegerVector>(payload["root_anchor_state_ids"]);
  // The full payload always carries the resolved root policy. Retain only the
  // semantic flag needed by traversal-time path scoring; C++ does not need to
  // interpret the policy's descriptive provenance strings.
  Rcpp::List root_policy = Rcpp::as<Rcpp::List>(payload["root_policy"]);
  input.root_is_synthetic = Rcpp::as<bool>(root_policy["synthetic"]);
  return input;
}

/**
 * Parse a simmap tree directly into the compact C++ input shape while preserving R-facing ids.
 *
 * Inputs: a mapped `phylo`, the stable state lookup, and the resolved root
 * anchor policy. Output: the same compact, validated `V35Input` used by the
 * payload parser. Ownership of copied vectors belongs to the returned input;
 * malformed map labels or edge-parallel fields cause an immediate R error.
 */
V35Input pst_v35_parse_tree_input(
    const Rcpp::List& tree,
    const Rcpp::DataFrame& state_lookup,
    const Rcpp::IntegerVector& root_anchor_state_ids,
    bool root_is_synthetic) {
  V35Input input;
  Rcpp::IntegerMatrix edge = tree["edge"];
  Rcpp::NumericVector edge_length = tree["edge.length"];
  Rcpp::CharacterVector tip_label = tree["tip.label"];
  Rcpp::List maps = tree["maps"];

  input.ntips = tip_label.size();
  input.nedge = edge.nrow();
  input.edge_parent = Rcpp::IntegerVector(input.nedge);
  input.edge_child = Rcpp::IntegerVector(input.nedge);
  // Copy each topology row into edge-parallel storage. After each iteration
  // the parent, child, and branch-length entries for that phylogeny edge share
  // the same zero-based C++ slot and retain their 1-based biological ids.
  for (int i = 0; i < input.nedge; ++i) {
    input.edge_parent[i] = edge(i, 0);
    input.edge_child[i] = edge(i, 1);
  }
  input.edge_length = edge_length;
  input.tip_label = tip_label;
  input.state_lookup = state_lookup;
  input.root_anchor_state_ids = root_anchor_state_ids;
  input.root_is_synthetic = root_is_synthetic;

  Rcpp::IntegerVector state_id_column = state_lookup["state_id"];
  Rcpp::CharacterVector state_column = state_lookup["state"];
  std::map<std::string, int> state_ids;
  // Index every declared state label once so map parsing can resolve labels to
  // stable integer ids without repeated data-frame searches.
  for (int i = 0; i < state_column.size(); ++i) {
    state_ids[Rcpp::as<std::string>(state_column[i])] = state_id_column[i];
  }

  int nsegment = 0;
  // Count mapped segments edge by edge before allocation. The cumulative
  // offset written for each edge is the exclusive start of its segment slice.
  for (int edge_id = 0; edge_id < maps.size(); ++edge_id) {
    Rcpp::NumericVector edge_map = maps[edge_id];
    nsegment += edge_map.size();
  }

  input.map_edge_id = Rcpp::IntegerVector(nsegment);
  input.map_state_id = Rcpp::IntegerVector(nsegment);
  input.map_duration = Rcpp::NumericVector(nsegment);
  input.map_edge_offset = Rcpp::IntegerVector(input.nedge + 1);
  input.map_edge_offset[0] = 1;

  int pos = 0;
  // Flatten every edge-local simmap into chronological segment arrays while
  // preserving edge ownership and the edge-local map order.
  for (int edge_id = 0; edge_id < maps.size(); ++edge_id) {
    Rcpp::NumericVector edge_map = maps[edge_id];
    Rcpp::CharacterVector names = edge_map.names();
    // Resolve each contiguous segment to its state id and duration. After each
    // iteration the flat arrays contain one complete biological segment.
    for (int segment_id = 0; segment_id < edge_map.size(); ++segment_id) {
      std::string state_name = Rcpp::as<std::string>(names[segment_id]);
      auto state_it = state_ids.find(state_name);
      // A map label absent from the supplied lookup would make every later
      // state/path counter ambiguous, so reject it before scheduling events.
      if (state_it == state_ids.end()) {
        Rcpp::stop("V35 C++ tree parser saw unknown state '%s'", state_name);
      }
      input.map_edge_id[pos] = edge_id + 1;
      input.map_state_id[pos] = state_it->second;
      input.map_duration[pos] = edge_map[segment_id];
      ++pos;
    }
    input.map_edge_offset[edge_id + 1] = pos + 1;
  }

  return input;
}
