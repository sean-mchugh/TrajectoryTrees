// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct TextValue {
  bool is_missing;
  std::string value;
};

struct TreeCoordinates {
  std::vector<double> node_times;
  std::vector<double> tip_end_times;
  std::vector<std::vector<int> > root_to_tip_edge_chains;
};

struct PreparedScenarioData {
  int tip_count;
  int interval_count;
  std::vector<std::string> tip_labels;
  SEXP tip_label_source;
  std::vector<double> interval_start_times;
  std::vector<double> interval_end_times;
  std::vector<double> interval_durations;
  std::vector<double> tip_end_times;
  std::vector<int> tip_end_columns;
  std::vector<TextValue> states;
  std::vector<TextValue> full_paths;
  std::vector<std::vector<TextValue> > paths_by_requested_depth;
  std::vector<std::vector<int> > path_ids_by_requested_depth;
  int projected_path_count;
  std::vector<int> available_transition_counts;
  std::vector<int> phylogenetic_edge_identifiers;
  std::vector<int> scenario_edge_identifiers;
  std::vector<int> scenario_edge_step_identifiers;
  std::vector<bool> lineage_is_active;
  std::vector<int> maximum_transition_counts;
  std::vector<std::string> state_levels;
  Rcpp::CharacterVector state_level_source;
  std::string root_anchor;
  double time_tolerance;
  SEXP root_policy;
  SEXP state_source;
  SEXP path_source;
  SEXP phylogenetic_edge_identifier_source;
  SEXP scenario_edge_identifier_source;
  SEXP scenario_edge_step_identifier_source;
  SEXP scenario_column_names;
  SEXP path_column_names;
  SEXP phylogenetic_edge_identifier_column_names;
  SEXP scenario_edge_identifier_column_names;
  SEXP scenario_edge_step_identifier_column_names;
};

struct ComponentMetadataRow {
  std::string matrix_name;
  TextValue state;
  TextValue path;
  std::string phase;
  std::string kind;
  int requested_maximum_transition_count;
  int actual_transition_count_used;
  std::string view;
};

struct ComponentDefinitions {
  std::vector<std::string> synchronous_component_names;
  std::vector<std::vector<std::string> > asynchronous_component_names;
  std::vector<ComponentMetadataRow> metadata_rows;
};

struct IntervalClassificationRecord {
  int interval_identifier;
  double start_time;
  double end_time;
  double duration;
  TextValue first_state;
  TextValue second_state;
  bool first_lineage_active;
  bool second_lineage_active;
  TextValue first_path;
  TextValue second_path;
  TextValue first_post_separation_path;
  TextValue second_post_separation_path;
  int first_phylogenetic_edge_identifier;
  int second_phylogenetic_edge_identifier;
  int first_scenario_edge_identifier;
  int second_scenario_edge_identifier;
  int first_scenario_step_identifier;
  int second_scenario_step_identifier;
  TextValue pair_mrca_state;
  bool root_is_synthetic;
  std::string phase;
  TextValue synchronous_class;
  std::string first_asynchronous_class;
  std::string second_asynchronous_class;
  int first_available_transition_count;
  int second_available_transition_count;
  int requested_maximum_transition_count;
  int actual_transition_count_used;
  TextValue matching_path;
};

struct SynchronousParallelTimeRecord {
  int requested_maximum_transition_count;
  int actual_transition_count_used;
  std::string path;
  std::string terminal_state;
  double duration;
};

struct AsynchronousIndependentPathTimeRecord {
  std::string lineage;
  std::string tip;
  int requested_maximum_transition_count;
  int available_transition_count;
  int actual_transition_count_used;
  std::string path;
  std::string terminal_state;
  double duration;
  double proportion;
};

struct AsynchronousParallelSimilarityRecord {
  int requested_maximum_transition_count;
  int actual_transition_count_used;
  std::string path;
  std::string terminal_state;
  double similarity;
  std::string metric;
};

struct PairCalculationResult {
  int first_tip_index;
  int second_tip_index;
  double total_comparable_time;
  double first_total_history_time;
  double second_total_history_time;
  std::vector<IntervalClassificationRecord> interval_classifications;
  std::map<std::string, double> synchronous_similarity;
  std::vector<std::map<std::string, double> > asynchronous_similarity;
  std::vector<SynchronousParallelTimeRecord> synchronous_parallel_times_by_path;
  std::vector<AsynchronousIndependentPathTimeRecord>
    asynchronous_independent_times_by_path;
  std::vector<AsynchronousParallelSimilarityRecord>
    asynchronous_parallel_similarities_by_path;
};

struct CompactLineageRun {
  int first_interval;
  int interval_after_last;
  int state;
  int phylogenetic_edge;
  std::vector<int> projected_path_ids;
};

struct CompactScenarioData {
  std::vector<std::vector<CompactLineageRun> > runs_by_tip;
  std::vector<double> cumulative_duration;
  int path_count;
};

struct CompactWorkspace {
  std::vector<double> first_shared;
  std::vector<double> second_shared;
  std::vector<double> first_conserved;
  std::vector<double> second_conserved;
  std::vector<double> first_independent;
  std::vector<double> second_independent;
  std::vector<double> shared_similarity;
  std::vector<double> conserved_similarity;
  std::vector<double> independent_similarity;
  std::vector<double> homoplasy_similarity;
  std::vector<double> parallel_similarity;
  std::vector<std::vector<double> > first_path_durations;
  std::vector<std::vector<double> > second_path_durations;
  std::vector<std::vector<std::uint64_t> > path_stamps;
  std::vector<std::vector<int> > touched_paths;
  std::vector<std::vector<int> > path_transition_counts;
  std::vector<std::vector<int> > path_terminal_states;
  std::uint64_t pair_stamp;
};

struct CompactMatrixStorage {
  std::vector<Rcpp::NumericMatrix> synchronous;
  std::vector<std::vector<Rcpp::NumericMatrix> > asynchronous;
  std::vector<double> synchronous_totals;
  std::vector<std::vector<double> > asynchronous_totals;
  Rcpp::NumericMatrix synchronous_available;
  Rcpp::NumericMatrix asynchronous_available;
  double synchronous_available_total;
  double asynchronous_available_total;
};

std::string translated_utf8(const Rcpp::String& value) {
  return std::string(Rf_translateCharUTF8(value.get_sexp()));
}

SEXP utf8_character(const std::string& value) {
  return Rf_mkCharCE(value.c_str(), CE_UTF8);
}

Rcpp::CharacterVector utf8_character_vector(
    const std::vector<std::string>& values) {
  Rcpp::CharacterVector result(values.size());
  for (std::size_t index = 0; index < values.size(); ++index) {
    result[index] = utf8_character(values[index]);
  }
  return result;
}

std::size_t scenario_cell_index(
    int tip_index,
    int interval_index,
    int interval_count) {
  return static_cast<std::size_t>(tip_index) *
    static_cast<std::size_t>(interval_count) +
    static_cast<std::size_t>(interval_index);
}

std::string replace_all(
    std::string value,
    const std::string& encoded,
    const std::string& decoded) {
  std::size_t position = 0;
  while ((position = value.find(encoded, position)) != std::string::npos) {
    value.replace(position, encoded.size(), decoded);
    position += decoded.size();
  }
  return value;
}

std::string decode_state_name(const std::string& value) {
  return replace_all(replace_all(value, "%7C", "|"), "%25", "%");
}

std::vector<std::string> path_components(const std::string& path) {
  if (path.empty() || path[0] != '|') {
    Rcpp::stop("Canonical paths must be nonmissing strings beginning with '|'");
  }

  std::vector<std::string> components;
  std::size_t component_start = 1;
  while (true) {
    const std::size_t separator = path.find('|', component_start);
    if (separator == std::string::npos) {
      components.push_back(path.substr(component_start));
      break;
    }
    components.push_back(path.substr(component_start, separator - component_start));
    component_start = separator + 1;
  }
  return components;
}

std::string terminal_state(const std::string& path) {
  const std::vector<std::string> components = path_components(path);
  return decode_state_name(components.back());
}

std::string project_transition_path(
    const std::vector<std::string>& full_path,
    int maximum_transition_count) {
  if (maximum_transition_count == 0) {
    return decode_state_name(full_path.back());
  }

  const int retained_component_count = std::min(
    static_cast<int>(full_path.size()),
    maximum_transition_count + 1
  );
  const int first_retained_component =
    static_cast<int>(full_path.size()) - retained_component_count;
  std::string projected_path;
  for (int component_index = first_retained_component;
       component_index < static_cast<int>(full_path.size());
       ++component_index) {
    projected_path += "|" + full_path[component_index];
  }
  return projected_path;
}

std::string transition_view_name(int maximum_transition_count) {
  if (maximum_transition_count == 0) return "state_only";
  return "using_up_to_" + std::to_string(maximum_transition_count) + "_" +
    (maximum_transition_count == 1 ? "transition" : "transitions");
}

std::string transition_count_label(int transition_count) {
  return std::to_string(transition_count) +
    (transition_count == 1 ? "Transition" : "Transitions");
}

std::string transition_depth_partition_label(
    int actual_transition_count,
    int requested_maximum_transition_count) {
  if (actual_transition_count == requested_maximum_transition_count) {
    return std::to_string(actual_transition_count) + "OrMoreTransitions";
  }
  return transition_count_label(actual_transition_count);
}

std::string maximum_transition_count_label(int maximum_transition_count) {
  return "UsingUpTo" + transition_count_label(maximum_transition_count);
}

TextValue text_matrix_value(
    const Rcpp::CharacterMatrix& matrix,
    int row_index,
    int column_index) {
  const Rcpp::String value = matrix(row_index, column_index);
  if (value == NA_STRING) return TextValue{true, std::string()};
  return TextValue{false, translated_utf8(value)};
}

SEXP matrix_column_names(SEXP matrix) {
  Rcpp::RObject matrix_object(matrix);
  if (!matrix_object.hasAttribute("dimnames")) return R_NilValue;
  Rcpp::List dimensions_names = matrix_object.attr("dimnames");
  if (dimensions_names.size() < 2) return R_NilValue;
  return dimensions_names[1];
}

void set_matrix_names(
    Rcpp::RObject matrix,
    SEXP tip_labels,
    SEXP column_names) {
  matrix.attr("dimnames") = Rcpp::List::create(
    tip_labels,
    column_names
  );
}

Rcpp::String source_tip_label(
    const PreparedScenarioData& prepared_scenario_data,
    int tip_index) {
  const Rcpp::CharacterVector labels(prepared_scenario_data.tip_label_source);
  return labels[tip_index];
}

Rcpp::String source_scenario_text(
    SEXP source_matrix,
    const PreparedScenarioData& prepared_scenario_data,
    int tip_index,
    int interval_index) {
  const Rcpp::CharacterMatrix matrix(source_matrix);
  return matrix(tip_index, interval_index);
}

TreeCoordinates calculate_tree_coordinates(const Rcpp::List& tree) {
  if (!tree.containsElementNamed("edge") ||
      !Rf_isMatrix(tree["edge"])) {
    Rcpp::stop("phylo$edge must be a finite positive two-column integer matrix");
  }
  const Rcpp::NumericMatrix edge = Rcpp::as<Rcpp::NumericMatrix>(tree["edge"]);
  if (edge.ncol() != 2) {
    Rcpp::stop("phylo$edge must be a finite positive two-column integer matrix");
  }
  for (int index = 0; index < edge.size(); ++index) {
    const double value = edge[index];
    if (!std::isfinite(value) || value < 1 ||
        value > std::numeric_limits<int>::max() ||
        value != std::floor(value)) {
      Rcpp::stop(
        "phylo$edge must contain representable positive integer node identifiers"
      );
    }
  }

  if (!tree.containsElementNamed("edge.length")) {
    Rcpp::stop("phylo$edge.length must contain one finite nonnegative length per edge");
  }
  const Rcpp::NumericVector edge_lengths =
    Rcpp::as<Rcpp::NumericVector>(tree["edge.length"]);
  if (edge_lengths.size() != edge.nrow()) {
    Rcpp::stop("phylo$edge.length must contain one finite nonnegative length per edge");
  }
  for (int edge_index = 0; edge_index < edge_lengths.size(); ++edge_index) {
    if (!std::isfinite(edge_lengths[edge_index]) || edge_lengths[edge_index] < 0) {
      Rcpp::stop("phylo$edge.length must contain one finite nonnegative length per edge");
    }
  }

  std::set<int> parent_nodes;
  std::set<int> child_nodes;
  int node_count = 0;
  for (int edge_index = 0; edge_index < edge.nrow(); ++edge_index) {
    const int parent_node = static_cast<int>(edge(edge_index, 0));
    const int child_node = static_cast<int>(edge(edge_index, 1));
    parent_nodes.insert(parent_node);
    child_nodes.insert(child_node);
    node_count = std::max(node_count, std::max(parent_node, child_node));
  }
  std::vector<int> root_nodes;
  std::set_difference(
    parent_nodes.begin(), parent_nodes.end(),
    child_nodes.begin(), child_nodes.end(),
    std::back_inserter(root_nodes)
  );
  if (root_nodes.size() != 1) {
    Rcpp::stop("phylogeny must have exactly one root");
  }
  const int root_node = root_nodes[0];

  std::vector<double> node_times(node_count + 1, NA_REAL);
  node_times[root_node] = 0;
  std::vector<int> unresolved_edges(edge.nrow());
  for (int edge_index = 0; edge_index < edge.nrow(); ++edge_index) {
    unresolved_edges[edge_index] = edge_index;
  }
  while (!unresolved_edges.empty()) {
    std::vector<int> still_unresolved;
    bool resolved_any_edge = false;
    for (std::size_t unresolved_index = 0;
         unresolved_index < unresolved_edges.size();
         ++unresolved_index) {
      const int edge_index = unresolved_edges[unresolved_index];
      const int parent_node = static_cast<int>(edge(edge_index, 0));
      const int child_node = static_cast<int>(edge(edge_index, 1));
      if (Rcpp::NumericVector::is_na(node_times[parent_node])) {
        still_unresolved.push_back(edge_index);
        continue;
      }

      const double candidate_time =
        node_times[parent_node] + edge_lengths[edge_index];
      if (!Rcpp::NumericVector::is_na(node_times[child_node]) &&
          node_times[child_node] != candidate_time) {
        Rcpp::stop("phylogeny child has multiple inconsistent parents");
      }
      node_times[child_node] = candidate_time;
      resolved_any_edge = true;
    }
    if (!resolved_any_edge) {
      Rcpp::stop("phylogeny edges do not form a rooted tree");
    }
    unresolved_edges.swap(still_unresolved);
  }

  if (!tree.containsElementNamed("tip.label")) {
    Rcpp::stop("trajectory_obj must contain at least two uniquely labelled tips");
  }
  const Rcpp::CharacterVector tip_labels = tree["tip.label"];
  std::vector<int> child_to_edge(node_count + 1, -1);
  for (int edge_index = 0; edge_index < edge.nrow(); ++edge_index) {
    child_to_edge[static_cast<int>(edge(edge_index, 1))] = edge_index;
  }

  std::vector<double> tip_end_times(tip_labels.size());
  std::vector<std::vector<int> > root_to_tip_edge_chains(tip_labels.size());
  for (int tip_index = 0; tip_index < tip_labels.size(); ++tip_index) {
    const int tip_node = tip_index + 1;
    tip_end_times[tip_index] = node_times[tip_node];
    int current_node = tip_node;
    std::vector<int> reversed_edge_chain;
    while (current_node != root_node) {
      if (current_node < 0 || current_node > node_count ||
          child_to_edge[current_node] < 0) {
        Rcpp::stop("a tip lacks a root-to-tip edge chain");
      }
      const int edge_index = child_to_edge[current_node];
      reversed_edge_chain.push_back(edge_index + 1);
      current_node = static_cast<int>(edge(edge_index, 0));
    }
    std::reverse(reversed_edge_chain.begin(), reversed_edge_chain.end());
    root_to_tip_edge_chains[tip_index] = reversed_edge_chain;
  }

  return TreeCoordinates{
    node_times,
    tip_end_times,
    root_to_tip_edge_chains
  };
}

bool integer_is_in_chain(int value, const std::vector<int>& edge_chain) {
  return std::find(edge_chain.begin(), edge_chain.end(), value) != edge_chain.end();
}

PreparedScenarioData prepare_scenario_data(
    const Rcpp::List& trajectory_obj,
    const std::vector<int>& maximum_transition_counts,
    bool record_complete_output,
    double time_tolerance) {
  const Rcpp::List tree = trajectory_obj["phylo"];
  const Rcpp::List scenario_matrices = trajectory_obj["scenario_mats"];
  const Rcpp::List root_policy = trajectory_obj["root_policy"];
  const Rcpp::CharacterVector tip_label_values = tree["tip.label"];
  const int tip_count = tip_label_values.size();
  if (tip_count < 2) {
    Rcpp::stop("trajectory_obj must contain at least two uniquely labelled tips");
  }

  std::vector<std::string> tip_labels(tip_count);
  std::set<std::string> unique_tip_labels;
  for (int tip_index = 0; tip_index < tip_count; ++tip_index) {
    const Rcpp::String tip_label = tip_label_values[tip_index];
    if (tip_label == NA_STRING) {
      Rcpp::stop("trajectory_obj must contain at least two uniquely labelled tips");
    }
    tip_labels[tip_index] = translated_utf8(tip_label);
    if (tip_labels[tip_index].empty() ||
        !unique_tip_labels.insert(tip_labels[tip_index]).second) {
      Rcpp::stop("trajectory_obj must contain at least two uniquely labelled tips");
    }
  }

  if (!std::isfinite(time_tolerance) || time_tolerance <= 0) {
    Rcpp::stop("time_tolerance must be one finite positive number");
  }

  if (!root_policy.containsElementNamed("anchor") ||
      root_policy["anchor"] == R_NilValue) {
    Rcpp::stop("root_policy$anchor must contain one biological state");
  }
  const Rcpp::CharacterVector root_anchor_value = root_policy["anchor"];
  if (root_anchor_value.size() != 1) {
    Rcpp::stop("root_policy$anchor must contain one biological state");
  }
  const Rcpp::String root_anchor_string = root_anchor_value[0];
  if (root_anchor_string == NA_STRING) {
    Rcpp::stop("root_policy$anchor must contain one biological state");
  }
  const std::string root_anchor(translated_utf8(root_anchor_string));
  if (root_anchor.empty()) {
    Rcpp::stop("root_policy$anchor must contain one biological state");
  }

  const std::vector<std::string> required_matrix_fields{
    "scenarios",
    "phylo_paths",
    "phylo_edge_ids",
    "scenario_edge_ids",
    "scenario_edge_step_ids"
  };
  for (std::size_t field_index = 0;
       field_index < required_matrix_fields.size();
       ++field_index) {
    const std::string& field_name = required_matrix_fields[field_index];
    if (!scenario_matrices.containsElementNamed(field_name.c_str()) ||
        !Rf_isMatrix(scenario_matrices[field_name])) {
      Rcpp::stop("scenario matrices are missing required matrix fields");
    }
  }

  const Rcpp::CharacterMatrix state_matrix = scenario_matrices["scenarios"];
  const Rcpp::CharacterMatrix path_matrix = scenario_matrices["phylo_paths"];
  const Rcpp::NumericMatrix phylogenetic_edge_identifier_matrix =
    Rcpp::as<Rcpp::NumericMatrix>(scenario_matrices["phylo_edge_ids"]);
  const Rcpp::NumericMatrix scenario_edge_identifier_matrix =
    Rcpp::as<Rcpp::NumericMatrix>(scenario_matrices["scenario_edge_ids"]);
  const Rcpp::NumericMatrix scenario_edge_step_identifier_matrix =
    Rcpp::as<Rcpp::NumericMatrix>(scenario_matrices["scenario_edge_step_ids"]);
  const int scenario_column_count = state_matrix.ncol();
  const int interval_count = scenario_column_count - 1;
  if (state_matrix.nrow() != tip_count ||
      path_matrix.nrow() != tip_count ||
      phylogenetic_edge_identifier_matrix.nrow() != tip_count ||
      scenario_edge_identifier_matrix.nrow() != tip_count ||
      scenario_edge_step_identifier_matrix.nrow() != tip_count ||
      path_matrix.ncol() != scenario_column_count ||
      phylogenetic_edge_identifier_matrix.ncol() != scenario_column_count ||
      scenario_edge_identifier_matrix.ncol() != scenario_column_count ||
      scenario_edge_step_identifier_matrix.ncol() != scenario_column_count) {
    Rcpp::stop("scenario matrices must have identical tip-by-time dimensions");
  }

  const std::vector<SEXP> matrix_sources{
    scenario_matrices["scenarios"],
    scenario_matrices["phylo_paths"],
    scenario_matrices["phylo_edge_ids"],
    scenario_matrices["scenario_edge_ids"],
    scenario_matrices["scenario_edge_step_ids"]
  };
  for (std::size_t matrix_index = 0;
       matrix_index < matrix_sources.size();
       ++matrix_index) {
    Rcpp::RObject matrix_object(matrix_sources[matrix_index]);
    if (matrix_object.hasAttribute("dimnames")) {
      Rcpp::List dimension_names = matrix_object.attr("dimnames");
      if (dimension_names.size() > 0 && dimension_names[0] != R_NilValue) {
        const Rcpp::CharacterVector row_names = dimension_names[0];
        if (row_names.size() != tip_count) {
          Rcpp::stop("scenario matrix rows must exactly match phylo tip labels");
        }
        for (int tip_index = 0; tip_index < tip_count; ++tip_index) {
          const Rcpp::String row_name = row_names[tip_index];
          if (row_name == NA_STRING ||
              translated_utf8(row_name) != tip_labels[tip_index]) {
            Rcpp::stop("scenario matrix rows must exactly match phylo tip labels");
          }
        }
      }
    }
  }

  if (!scenario_matrices.containsElementNamed("time_vec")) {
    Rcpp::stop("scenario_mats$time_vec must be finite, column-aligned, and strictly increasing");
  }
  const Rcpp::NumericVector time_values = scenario_matrices["time_vec"];
  if (time_values.size() != scenario_column_count || scenario_column_count < 2) {
    Rcpp::stop("scenario_mats$time_vec must be finite, column-aligned, and strictly increasing");
  }
  std::vector<double> time_vector(scenario_column_count);
  for (int column_index = 0;
       column_index < scenario_column_count;
       ++column_index) {
    time_vector[column_index] = time_values[column_index];
    if (!std::isfinite(time_vector[column_index]) ||
        (column_index > 0 &&
         time_vector[column_index] <= time_vector[column_index - 1])) {
      Rcpp::stop("scenario_mats$time_vec must be finite, column-aligned, and strictly increasing");
    }
  }

  const TreeCoordinates tree_coordinates = calculate_tree_coordinates(tree);
  const double time_scale = 1.0 / time_tolerance;
  const auto time_key = [time_scale](double time) {
    const long double scaled_time = static_cast<long double>(time) *
      static_cast<long double>(time_scale);
    if (scaled_time <
          static_cast<long double>(std::numeric_limits<long long>::min()) ||
        scaled_time >
          static_cast<long double>(std::numeric_limits<long long>::max())) {
      Rcpp::stop("trajectory time exceeds the integer-key range at its declared tolerance");
    }
    return static_cast<long long>(std::llrint(time * time_scale));
  };
  std::vector<long long> time_vector_keys(scenario_column_count);
  for (int column_index = 0;
       column_index < scenario_column_count;
       ++column_index) {
    time_vector_keys[column_index] = time_key(time_vector[column_index]);
  }

  std::vector<int> tip_end_columns(tip_count);
  for (int tip_index = 0; tip_index < tip_count; ++tip_index) {
    const long long tip_end_time_key = time_key(
      tree_coordinates.tip_end_times[tip_index]
    );
    int matching_column_count = 0;
    int matching_column_index = -1;
    double closest_time = time_vector[0];
    double closest_distance = std::fabs(
      time_vector[0] - tree_coordinates.tip_end_times[tip_index]
    );
    double second_closest_time = NA_REAL;
    double second_closest_distance = std::numeric_limits<double>::infinity();
    for (int column_index = 1;
         column_index < scenario_column_count;
         ++column_index) {
      const double distance = std::fabs(
        time_vector[column_index] - tree_coordinates.tip_end_times[tip_index]
      );
      if (distance < closest_distance) {
        second_closest_distance = closest_distance;
        second_closest_time = closest_time;
        closest_distance = distance;
        closest_time = time_vector[column_index];
      } else if (distance < second_closest_distance &&
                 time_vector[column_index] != closest_time) {
        second_closest_distance = distance;
        second_closest_time = time_vector[column_index];
      }
    }
    for (int column_index = 0;
         column_index < scenario_column_count;
         ++column_index) {
      if (time_vector_keys[column_index] == tip_end_time_key) {
        matching_column_index = column_index;
        ++matching_column_count;
      }
    }
    if (matching_column_count == 0) {
      std::ostringstream message;
      message.precision(17);
      message << "tip '" << tip_labels[tip_index] << "' termination time "
              << tree_coordinates.tip_end_times[tip_index]
              << " has 0 matching time_vec columns; closest time "
              << closest_time
              << ", absolute distance " << closest_distance
              << ", tolerance " << time_tolerance
              << "; second closest time " << second_closest_time
              << ", absolute distance " << second_closest_distance;
      Rcpp::stop(message.str());
    }
    if (matching_column_count > 1) {
      std::ostringstream message;
      message.precision(17);
      message << "tip '" << tip_labels[tip_index] << "' termination time "
              << tree_coordinates.tip_end_times[tip_index]
              << " maps to multiple time_vec columns at the declared time tolerance; "
              << "closest time "
              << closest_time << ", absolute distance " << closest_distance
              << ", tolerance " << time_tolerance
              << "; second closest time " << second_closest_time
              << ", absolute distance " << second_closest_distance;
      Rcpp::stop(message.str());
    }
    tip_end_columns[tip_index] = matching_column_index + 1;
  }

  std::vector<double> interval_start_times(interval_count);
  std::vector<double> interval_end_times(interval_count);
  std::vector<double> interval_durations(interval_count);
  for (int interval_index = 0; interval_index < interval_count; ++interval_index) {
    interval_start_times[interval_index] = time_vector[interval_index];
    interval_end_times[interval_index] = time_vector[interval_index + 1];
    interval_durations[interval_index] =
      interval_end_times[interval_index] - interval_start_times[interval_index];
  }

  const std::size_t scenario_cell_count =
    static_cast<std::size_t>(tip_count) *
    static_cast<std::size_t>(interval_count);
  std::vector<TextValue> states(scenario_cell_count);
  std::vector<TextValue> full_paths(scenario_cell_count);
  std::vector<int> phylogenetic_edge_identifiers(scenario_cell_count, NA_INTEGER);
  std::vector<int> scenario_edge_identifiers(scenario_cell_count, NA_INTEGER);
  std::vector<int> scenario_edge_step_identifiers(scenario_cell_count, NA_INTEGER);
  std::vector<bool> lineage_is_active(scenario_cell_count, false);

  for (int tip_index = 0; tip_index < tip_count; ++tip_index) {
    for (int interval_index = 0;
         interval_index < interval_count;
         ++interval_index) {
      const std::size_t cell_index = scenario_cell_index(
        tip_index,
        interval_index,
        interval_count
      );
      const bool is_active = interval_index < tip_end_columns[tip_index] - 1;
      lineage_is_active[cell_index] = is_active;
      states[cell_index] = text_matrix_value(state_matrix, tip_index, interval_index);
      full_paths[cell_index] = text_matrix_value(path_matrix, tip_index, interval_index);

      const double phylogenetic_edge_identifier =
        phylogenetic_edge_identifier_matrix(tip_index, interval_index);
      const double scenario_edge_identifier =
        scenario_edge_identifier_matrix(tip_index, interval_index);
      const double scenario_edge_step_identifier =
        scenario_edge_step_identifier_matrix(tip_index, interval_index);

      if (is_active &&
          (states[cell_index].is_missing || states[cell_index].value.empty() ||
           full_paths[cell_index].is_missing || full_paths[cell_index].value.empty())) {
        Rcpp::stop("active intervals require nonmissing states and canonical paths");
      }
      if (is_active &&
          (!std::isfinite(phylogenetic_edge_identifier) ||
           !std::isfinite(scenario_edge_identifier) ||
           !std::isfinite(scenario_edge_step_identifier) ||
           phylogenetic_edge_identifier < 1 ||
           scenario_edge_identifier < 1 ||
           scenario_edge_step_identifier < 1 ||
           phylogenetic_edge_identifier > std::numeric_limits<int>::max() ||
           scenario_edge_identifier > std::numeric_limits<int>::max() ||
           scenario_edge_step_identifier > std::numeric_limits<int>::max() ||
           phylogenetic_edge_identifier != std::floor(phylogenetic_edge_identifier) ||
           scenario_edge_identifier != std::floor(scenario_edge_identifier) ||
           scenario_edge_step_identifier != std::floor(scenario_edge_step_identifier))) {
        Rcpp::stop(
          "active edge and step ids must be representable positive integers"
        );
      }

      if (std::isfinite(phylogenetic_edge_identifier) &&
          phylogenetic_edge_identifier >= std::numeric_limits<int>::min() &&
          phylogenetic_edge_identifier <= std::numeric_limits<int>::max()) {
        phylogenetic_edge_identifiers[cell_index] =
          static_cast<int>(phylogenetic_edge_identifier);
      }
      if (std::isfinite(scenario_edge_identifier) &&
          scenario_edge_identifier >= std::numeric_limits<int>::min() &&
          scenario_edge_identifier <= std::numeric_limits<int>::max()) {
        scenario_edge_identifiers[cell_index] =
          static_cast<int>(scenario_edge_identifier);
      }
      if (std::isfinite(scenario_edge_step_identifier) &&
          scenario_edge_step_identifier >= std::numeric_limits<int>::min() &&
          scenario_edge_step_identifier <= std::numeric_limits<int>::max()) {
        scenario_edge_step_identifiers[cell_index] =
          static_cast<int>(scenario_edge_step_identifier);
      }

      if (is_active) {
        if (!integer_is_in_chain(
              phylogenetic_edge_identifiers[cell_index],
              tree_coordinates.root_to_tip_edge_chains[tip_index])) {
          Rcpp::stop("an active phylo edge id is outside its tip's root-to-tip chain");
        }
        const std::string decoded_terminal_state =
          terminal_state(full_paths[cell_index].value);
        if (decoded_terminal_state != states[cell_index].value) {
          std::ostringstream message;
          message.precision(17);
          message << "active path/state mismatch for tip '"
                  << tip_labels[tip_index]
                  << "', interval " << interval_index + 1
                  << " [" << interval_start_times[interval_index]
                  << ", " << interval_end_times[interval_index]
                  << "): state '" << states[cell_index].value
                  << "', path '" << full_paths[cell_index].value
                  << "', decoded terminal state '" << decoded_terminal_state
                  << "'";
          Rcpp::stop(message.str());
        }
      }
    }

    for (int interval_index = 0;
         interval_index < interval_count;
         ++interval_index) {
      if (interval_start_times[interval_index] <
            tree_coordinates.tip_end_times[tip_index] - time_tolerance &&
          interval_end_times[interval_index] >
            tree_coordinates.tip_end_times[tip_index] + time_tolerance) {
        Rcpp::stop("a tip termination clips a scenario interval");
      }
    }
  }

  if (tree.containsElementNamed("maps") && tree["maps"] != R_NilValue) {
    const Rcpp::List maps = tree["maps"];
    const Rcpp::NumericMatrix edge = Rcpp::as<Rcpp::NumericMatrix>(tree["edge"]);
    if (maps.size() == edge.nrow()) {
      for (int tip_index = 0; tip_index < tip_count; ++tip_index) {
        const int terminal_edge_identifier =
          tree_coordinates.root_to_tip_edge_chains[tip_index].back();
        const Rcpp::NumericVector terminal_map = maps[terminal_edge_identifier - 1];
        const int endpoint_column = tip_end_columns[tip_index];
        if (terminal_map.size() > 0 &&
            terminal_map[terminal_map.size() - 1] > time_tolerance &&
            endpoint_column > 1) {
          const int previous_column_index = endpoint_column - 2;
          const int endpoint_column_index = endpoint_column - 1;
          const TextValue previous_state = text_matrix_value(
            state_matrix,
            tip_index,
            previous_column_index
          );
          const TextValue endpoint_state = text_matrix_value(
            state_matrix,
            tip_index,
            endpoint_column_index
          );
          const TextValue previous_path = text_matrix_value(
            path_matrix,
            tip_index,
            previous_column_index
          );
          const TextValue endpoint_path = text_matrix_value(
            path_matrix,
            tip_index,
            endpoint_column_index
          );
          if (previous_state.is_missing != endpoint_state.is_missing ||
              previous_state.value != endpoint_state.value ||
              previous_path.is_missing != endpoint_path.is_missing ||
              previous_path.value != endpoint_path.value ||
              phylogenetic_edge_identifier_matrix(tip_index, previous_column_index) !=
                phylogenetic_edge_identifier_matrix(tip_index, endpoint_column_index)) {
            Rcpp::stop("positive terminal history does not match its endpoint snapshot");
          }
        }
      }
    }
  }

  std::vector<int> available_transition_counts(scenario_cell_count, NA_INTEGER);
  std::vector<std::vector<TextValue> > paths_by_requested_depth(
    maximum_transition_counts.size()
  );
  if (record_complete_output) {
    for (std::size_t depth = 0;
         depth < paths_by_requested_depth.size();
         ++depth) {
      paths_by_requested_depth[depth].resize(scenario_cell_count);
    }
  }
  std::vector<std::vector<int> > path_ids_by_requested_depth(
    maximum_transition_counts.size(),
    std::vector<int>(scenario_cell_count, -1)
  );
  std::unordered_map<std::string, int> projected_path_ids;
  std::unordered_map<
    std::string,
    std::pair<int, std::vector<int> >
  > compact_path_cache;
  compact_path_cache.reserve(tip_count);
  for (std::size_t cell_index = 0;
       cell_index < scenario_cell_count;
       ++cell_index) {
    if (full_paths[cell_index].is_missing) {
      if (record_complete_output) {
        for (std::size_t depth_index = 0;
             depth_index < maximum_transition_counts.size();
             ++depth_index) {
          paths_by_requested_depth[depth_index][cell_index] =
            TextValue{true, std::string()};
        }
      }
      continue;
    }

    if (!record_complete_output) {
      const std::unordered_map<
        std::string,
        std::pair<int, std::vector<int> >
      >::const_iterator cached =
        compact_path_cache.find(full_paths[cell_index].value);
      if (cached != compact_path_cache.end()) {
        available_transition_counts[cell_index] = cached->second.first;
        for (std::size_t depth = 0;
             depth < maximum_transition_counts.size();
             ++depth) {
          path_ids_by_requested_depth[depth][cell_index] =
            cached->second.second[depth];
        }
        continue;
      }
    }

    const std::vector<std::string> components =
      path_components(full_paths[cell_index].value);
    available_transition_counts[cell_index] =
      static_cast<int>(components.size()) - 1;
    for (std::size_t depth_index = 0;
         depth_index < maximum_transition_counts.size();
         ++depth_index) {
      const std::string path = project_transition_path(
        components, maximum_transition_counts[depth_index]
      );
      if (record_complete_output) {
        paths_by_requested_depth[depth_index][cell_index] =
          TextValue{false, path};
      }
      std::unordered_map<std::string, int>::iterator position =
        projected_path_ids.find(path);
      if (position == projected_path_ids.end()) {
        position = projected_path_ids.insert(
          std::make_pair(path, static_cast<int>(projected_path_ids.size()))
        ).first;
      }
      path_ids_by_requested_depth[depth_index][cell_index] = position->second;
    }
    if (!record_complete_output) {
      std::vector<int> ids(maximum_transition_counts.size());
      for (std::size_t depth = 0;
           depth < maximum_transition_counts.size();
           ++depth) {
        ids[depth] = path_ids_by_requested_depth[depth][cell_index];
      }
      compact_path_cache[full_paths[cell_index].value] =
        std::make_pair(available_transition_counts[cell_index], ids);
    }
  }

  std::unordered_map<std::string, bool> unique_states;
  for (std::size_t cell_index = 0;
       cell_index < lineage_is_active.size();
       ++cell_index) {
    if (lineage_is_active[cell_index]) {
      unique_states[states[cell_index].value] = true;
    }
  }
  Rcpp::CharacterVector unique_state_values(unique_states.size());
  R_xlen_t state_index = 0;
  for (std::unordered_map<std::string, bool>::const_iterator position =
         unique_states.begin();
       position != unique_states.end();
       ++position, ++state_index) {
    unique_state_values[state_index] = utf8_character(position->first);
  }
  const Rcpp::Environment base_environment = Rcpp::Environment::base_env();
  const Rcpp::Function sort_function = base_environment["sort"];
  const Rcpp::CharacterVector state_level_source = sort_function(
    unique_state_values
  );
  std::vector<std::string> state_levels(state_level_source.size());
  for (R_xlen_t state_index = 0;
       state_index < state_level_source.size();
       ++state_index) {
    state_levels[state_index] = translated_utf8(state_level_source[state_index]);
  }
  return PreparedScenarioData{
    tip_count,
    interval_count,
    tip_labels,
    tip_label_values,
    interval_start_times,
    interval_end_times,
    interval_durations,
    tree_coordinates.tip_end_times,
    tip_end_columns,
    states,
    full_paths,
    paths_by_requested_depth,
    path_ids_by_requested_depth,
    static_cast<int>(projected_path_ids.size()),
    available_transition_counts,
    phylogenetic_edge_identifiers,
    scenario_edge_identifiers,
    scenario_edge_step_identifiers,
    lineage_is_active,
    maximum_transition_counts,
    state_levels,
    state_level_source,
    root_anchor,
    time_tolerance,
    root_policy,
    scenario_matrices["scenarios"],
    scenario_matrices["phylo_paths"],
    scenario_matrices["phylo_edge_ids"],
    scenario_matrices["scenario_edge_ids"],
    scenario_matrices["scenario_edge_step_ids"],
    matrix_column_names(scenario_matrices["scenarios"]),
    matrix_column_names(scenario_matrices["phylo_paths"]),
    matrix_column_names(scenario_matrices["phylo_edge_ids"]),
    matrix_column_names(scenario_matrices["scenario_edge_ids"]),
    matrix_column_names(scenario_matrices["scenario_edge_step_ids"])
  };
}

void add_component_metadata_row(
    ComponentDefinitions& component_definitions,
    const std::string& matrix_name,
    const TextValue& state,
    const std::string& phase,
    const std::string& kind,
    int requested_maximum_transition_count,
    int actual_transition_count_used,
    const std::string& view) {
  component_definitions.metadata_rows.push_back(ComponentMetadataRow{
    matrix_name,
    state,
    TextValue{true, std::string()},
    phase,
    kind,
    requested_maximum_transition_count,
    actual_transition_count_used,
    view
  });
}

ComponentDefinitions build_component_definitions(
    const PreparedScenarioData& prepared_scenario_data) {
  ComponentDefinitions component_definitions;

  for (std::size_t state_index = 0;
       state_index < prepared_scenario_data.state_levels.size();
       ++state_index) {
    const std::string& state = prepared_scenario_data.state_levels[state_index];
    component_definitions.synchronous_component_names.push_back(
      state + "_Shared"
    );
    component_definitions.synchronous_component_names.push_back(
      state + "_Conserved"
    );
    component_definitions.synchronous_component_names.push_back(
      state + "_Conserved_Homoplasy"
    );
    component_definitions.synchronous_component_names.push_back(
      state + "_Convergent"
    );
    add_component_metadata_row(
      component_definitions,
      state + "_Shared",
      TextValue{false, state},
      "Shared",
      "state",
      NA_INTEGER,
      NA_INTEGER,
      "sync"
    );
    add_component_metadata_row(
      component_definitions,
      state + "_Conserved",
      TextValue{false, state},
      "Conserved",
      "state",
      NA_INTEGER,
      NA_INTEGER,
      "sync"
    );
    add_component_metadata_row(
      component_definitions,
      state + "_Conserved_Homoplasy",
      TextValue{false, state},
      "Conserved_Homoplasy",
      "state",
      NA_INTEGER,
      NA_INTEGER,
      "sync"
    );
    add_component_metadata_row(
      component_definitions,
      state + "_Convergent",
      TextValue{false, state},
      "Independent",
      "state",
      NA_INTEGER,
      0,
      "sync"
    );
  }

  const int largest_requested_maximum = *std::max_element(
    prepared_scenario_data.maximum_transition_counts.begin(),
    prepared_scenario_data.maximum_transition_counts.end()
  );
  for (int actual_transition_count = 1;
       actual_transition_count <= largest_requested_maximum;
       ++actual_transition_count) {
    for (std::size_t state_index = 0;
         state_index < prepared_scenario_data.state_levels.size();
         ++state_index) {
      const std::string& state = prepared_scenario_data.state_levels[state_index];
      const std::string matrix_name = state + "_Parallel_" +
        transition_count_label(actual_transition_count);
      component_definitions.synchronous_component_names.push_back(matrix_name);
      add_component_metadata_row(
        component_definitions,
        matrix_name,
        TextValue{false, state},
        "Independent",
        "parallel",
        NA_INTEGER,
        actual_transition_count,
        "sync"
      );
    }
  }
  component_definitions.synchronous_component_names.push_back("Divergent");
  add_component_metadata_row(
    component_definitions,
    "Divergent",
    TextValue{true, std::string()},
    "Independent",
    "residual",
    NA_INTEGER,
    NA_INTEGER,
    "sync"
  );

  component_definitions.asynchronous_component_names.resize(
    prepared_scenario_data.maximum_transition_counts.size()
  );
  for (std::size_t depth_index = 0;
       depth_index < prepared_scenario_data.maximum_transition_counts.size();
       ++depth_index) {
    const int maximum_transition_count =
      prepared_scenario_data.maximum_transition_counts[depth_index];
    const std::string view_name = transition_view_name(maximum_transition_count);
    std::vector<std::string>& component_names =
      component_definitions.asynchronous_component_names[depth_index];
    for (std::size_t state_index = 0;
         state_index < prepared_scenario_data.state_levels.size();
         ++state_index) {
      const std::string& state = prepared_scenario_data.state_levels[state_index];
      component_names.push_back(state + "_Shared");
      component_names.push_back(state + "_Conserved");
      component_names.push_back(state + "_Conserved_Homoplasy");
      add_component_metadata_row(
        component_definitions,
        state + "_Shared",
        TextValue{false, state},
        "Shared",
        "state",
        maximum_transition_count,
        NA_INTEGER,
        view_name
      );
      add_component_metadata_row(
        component_definitions,
        state + "_Conserved",
        TextValue{false, state},
        "Conserved",
        "state",
        maximum_transition_count,
        NA_INTEGER,
        view_name
      );
      add_component_metadata_row(
        component_definitions,
        state + "_Conserved_Homoplasy",
        TextValue{false, state},
        "Conserved_Homoplasy",
        "state",
        maximum_transition_count,
        NA_INTEGER,
        view_name
      );
      if (maximum_transition_count == 0) {
        component_names.push_back(state + "_Convergent");
        add_component_metadata_row(
          component_definitions,
          state + "_Convergent",
          TextValue{false, state},
          "Independent",
          "convergent",
          0,
          0,
          view_name
        );
      } else {
        for (int actual_transition_count = 1;
             actual_transition_count <= maximum_transition_count;
             ++actual_transition_count) {
          const std::string parallel_name = state + "_Parallel_" +
            transition_depth_partition_label(
              actual_transition_count,
              maximum_transition_count
            );
          component_names.push_back(parallel_name);
          add_component_metadata_row(
            component_definitions,
            parallel_name,
            TextValue{false, state},
            "Independent",
            "parallel",
            maximum_transition_count,
            actual_transition_count,
            view_name
          );
        }
        const std::string convergent_name = state + "_Convergent";
        component_names.push_back(convergent_name);
        add_component_metadata_row(
          component_definitions,
          convergent_name,
          TextValue{false, state},
          "Independent",
          "convergent",
          maximum_transition_count,
          0,
          view_name
        );
      }
    }
    component_names.push_back("Divergent");
    add_component_metadata_row(
      component_definitions,
      "Divergent",
      TextValue{true, std::string()},
      "Independent",
      "residual",
      maximum_transition_count,
      NA_INTEGER,
      view_name
    );
  }

  return component_definitions;
}

void add_named_amount(
    std::map<std::string, double>& values,
    const std::string& name,
    double amount) {
  values[name] += amount;
}

double named_amount(
    const std::map<std::string, double>& values,
    const std::string& name) {
  const std::map<std::string, double>::const_iterator position = values.find(name);
  return position == values.end() ? 0 : position->second;
}

double sum_named_amounts(const std::map<std::string, double>& values) {
  double total = 0;
  for (std::map<std::string, double>::const_iterator position = values.begin();
       position != values.end();
       ++position) {
    total += position->second;
  }
  return total;
}

double calculate_similarity_value(
    double first_proportion,
    double second_proportion,
    const std::string& async_metric) {
  if (async_metric == "bhattacharyya") {
    return std::sqrt(first_proportion * second_proportion);
  }
  return std::min(first_proportion, second_proportion);
}

std::map<std::string, double> calculate_common_similarity(
    const std::map<std::string, double>& first_durations,
    const std::map<std::string, double>& second_durations,
    double total_comparable_time,
    const std::string& async_metric) {
  std::set<std::string> names;
  for (std::map<std::string, double>::const_iterator position =
         first_durations.begin();
       position != first_durations.end();
       ++position) {
    names.insert(position->first);
  }
  for (std::map<std::string, double>::const_iterator position =
         second_durations.begin();
       position != second_durations.end();
       ++position) {
    names.insert(position->first);
  }

  std::map<std::string, double> similarity;
  for (std::set<std::string>::const_iterator name = names.begin();
       name != names.end();
       ++name) {
    similarity[*name] = calculate_similarity_value(
      named_amount(first_durations, *name) / total_comparable_time,
      named_amount(second_durations, *name) / total_comparable_time,
      async_metric
    );
  }
  return similarity;
}

void guard_similarity_composition(
    std::map<std::string, double>& values,
    double time_tolerance) {
  double total = 0;
  for (std::map<std::string, double>::iterator position = values.begin();
       position != values.end();
       ++position) {
    if (position->second < -time_tolerance ||
        position->second > 1 + time_tolerance) {
      Rcpp::stop("similarity component fell outside [0, 1]");
    }
    if (position->second < 0) position->second = 0;
    total += position->second;
  }
  if (total > 1 + time_tolerance) {
    Rcpp::stop("similarity components exceed unit mass");
  }
  if (std::fabs(total - 1) > time_tolerance) {
    Rcpp::stop("similarity components do not sum to unit mass");
  }
  if (total != 1) {
    for (std::map<std::string, double>::iterator position = values.begin();
         position != values.end();
         ++position) {
      position->second /= total;
    }
  }
}

void guard_similarity_composition(
    std::vector<double>& values,
    double time_tolerance) {
  double total = 0.0;
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (values[index] < -time_tolerance ||
        values[index] > 1.0 + time_tolerance) {
      Rcpp::stop("similarity component fell outside [0, 1]");
    }
    values[index] = std::max(0.0, values[index]);
    total += values[index];
  }
  if (std::fabs(total - 1.0) > time_tolerance) {
    Rcpp::stop("similarity components do not sum to unit mass");
  }
  if (total != 1.0) {
    for (std::size_t index = 0; index < values.size(); ++index) {
      values[index] /= total;
    }
  }
}

double duration_similarity(
    double first_duration,
    double second_duration,
    double first_history_time,
    double second_history_time,
    const std::string& metric) {
  if (metric == "bhattacharyya") {
    return std::sqrt(first_duration * second_duration) /
      std::sqrt(first_history_time * second_history_time);
  }
  return std::min(first_duration, second_duration) /
    std::min(first_history_time, second_history_time);
}

double compact_duration_similarity(
    double first_duration,
    double second_duration,
    double denominator,
    bool bhattacharyya) {
  return (bhattacharyya ?
    std::sqrt(first_duration * second_duration) :
    std::min(first_duration, second_duration)) / denominator;
}

bool root_is_synthetic(const PreparedScenarioData& data) {
  const Rcpp::List policy(data.root_policy);
  if (policy.containsElementNamed("synthetic") &&
      policy["synthetic"] != R_NilValue &&
      Rcpp::as<bool>(policy["synthetic"])) {
    return true;
  }
  if (!policy.containsElementNamed("provenance") ||
      policy["provenance"] == R_NilValue) {
    return false;
  }
  return Rcpp::as<std::string>(policy["provenance"]) == "synthetic-root";
}

CompactScenarioData build_compact_scenario_data(
    const PreparedScenarioData& data) {
  std::unordered_map<std::string, int> state_ids;
  for (std::size_t state = 0; state < data.state_levels.size(); ++state)
    state_ids[data.state_levels[state]] = state;
  const int maximum_depth = *std::max_element(
    data.maximum_transition_counts.begin(), data.maximum_transition_counts.end()
  );
  std::vector<int> prepared_depth_indices(maximum_depth + 1, -1);
  for (std::size_t index = 0; index < data.maximum_transition_counts.size();
       ++index) {
    prepared_depth_indices[data.maximum_transition_counts[index]] = index;
  }
  const bool all_depths_prepared = std::find(
    prepared_depth_indices.begin(), prepared_depth_indices.end(), -1) ==
    prepared_depth_indices.end();
  std::unordered_map<std::string, int> path_ids;
  std::unordered_map<std::string, std::vector<int> > projected_ids_by_path;
  path_ids.reserve(data.tip_count);
  projected_ids_by_path.reserve(data.tip_count);

  CompactScenarioData compact;
  compact.runs_by_tip.resize(data.tip_count);
  for (int tip = 0; tip < data.tip_count; ++tip) {
    std::vector<CompactLineageRun>& runs = compact.runs_by_tip[tip];
    std::vector<int> direct_path_ids(maximum_depth + 1);
    for (int interval = 0; interval < data.interval_count; ++interval) {
      const std::size_t cell =
        scenario_cell_index(tip, interval, data.interval_count);
      if (!data.lineage_is_active[cell]) continue;
      const int state = state_ids[data.states[cell].value];
      const int edge = data.phylogenetic_edge_identifiers[cell];
      const std::vector<int>* projected_path_ids;
      if (all_depths_prepared) {
        for (int depth = 0; depth <= maximum_depth; ++depth) {
          direct_path_ids[depth] = data.path_ids_by_requested_depth[
            prepared_depth_indices[depth]
          ][cell];
        }
        projected_path_ids = &direct_path_ids;
      } else {
        const std::string& full_path = data.full_paths[cell].value;
        std::unordered_map<std::string, std::vector<int> >::iterator cached =
          projected_ids_by_path.find(full_path);
        if (cached == projected_ids_by_path.end()) {
          const std::vector<std::string> components = path_components(full_path);
          std::vector<int> ids(maximum_depth + 1);
          for (int depth = 0; depth <= maximum_depth; ++depth) {
            const std::string path = project_transition_path(components, depth);
            std::unordered_map<std::string, int>::iterator position =
              path_ids.find(path);
            if (position == path_ids.end()) {
              position = path_ids.insert(std::make_pair(
                path, static_cast<int>(path_ids.size()))).first;
            }
            ids[depth] = position->second;
          }
          cached = projected_ids_by_path.insert(
            std::make_pair(full_path, ids)).first;
        }
        projected_path_ids = &cached->second;
      }
      if (!runs.empty() && runs.back().interval_after_last == interval &&
          runs.back().state == state &&
          runs.back().phylogenetic_edge == edge &&
          runs.back().projected_path_ids == *projected_path_ids) {
        runs.back().interval_after_last = interval + 1;
      } else {
        runs.push_back(CompactLineageRun{
          interval, interval + 1, state, edge, *projected_path_ids});
      }
    }
    const int expected_end = data.tip_end_columns[tip] - 1;
    if (runs.empty() || runs.front().first_interval != 0 ||
        runs.back().interval_after_last != expected_end)
      Rcpp::stop("compact lineage runs do not cover the active history");
  }
  compact.cumulative_duration.assign(data.interval_count + 1, 0.0);
  for (int interval = 0; interval < data.interval_count; ++interval) {
    compact.cumulative_duration[interval + 1] = compact.cumulative_duration[
      interval] + data.interval_durations[interval];
  }
  compact.path_count = all_depths_prepared ? data.projected_path_count : path_ids.size();
  return compact;
}

void add_compact_path_duration(
    CompactWorkspace& workspace,
    std::size_t depth_index,
    int path_id,
    int transition_count,
    int terminal_state,
    double duration,
    bool first_lineage) {
  if (workspace.path_stamps[depth_index][path_id] != workspace.pair_stamp) {
    workspace.path_stamps[depth_index][path_id] = workspace.pair_stamp;
    workspace.first_path_durations[depth_index][path_id] = 0.0;
    workspace.second_path_durations[depth_index][path_id] = 0.0;
    workspace.path_transition_counts[depth_index][path_id] = transition_count;
    workspace.path_terminal_states[depth_index][path_id] = terminal_state;
    workspace.touched_paths[depth_index].push_back(path_id);
  }
  std::vector<double>& durations = first_lineage ?
    workspace.first_path_durations[depth_index] :
    workspace.second_path_durations[depth_index];
  durations[path_id] += duration;
}

void calculate_compact_pair_similarities(
    const PreparedScenarioData& data,
    const CompactScenarioData& compact,
    int first_tip,
    int second_tip,
    bool synthetic_root,
    int root_anchor_state,
    const std::string& metric,
    CompactWorkspace& workspace,
    std::vector<double>& synchronous,
    std::vector<std::vector<double> >& asynchronous) {
  std::fill(synchronous.begin(), synchronous.end(), 0.0);
  std::fill(workspace.first_shared.begin(), workspace.first_shared.end(), 0.0);
  std::fill(workspace.second_shared.begin(), workspace.second_shared.end(), 0.0);
  std::fill(workspace.first_conserved.begin(),
            workspace.first_conserved.end(), 0.0);
  std::fill(workspace.second_conserved.begin(),
            workspace.second_conserved.end(), 0.0);
  std::fill(workspace.first_independent.begin(),
            workspace.first_independent.end(), 0.0);
  std::fill(workspace.second_independent.begin(),
            workspace.second_independent.end(), 0.0);
  ++workspace.pair_stamp;
  for (std::size_t depth = 0; depth < asynchronous.size(); ++depth) {
    workspace.touched_paths[depth].clear();
  }

  const int first_end = data.tip_end_columns[first_tip] - 1;
  const int second_end = data.tip_end_columns[second_tip] - 1;
  const int end = std::max(first_end, second_end);
  const std::vector<CompactLineageRun>& first_runs =
    compact.runs_by_tip[first_tip];
  const std::vector<CompactLineageRun>& second_runs =
    compact.runs_by_tip[second_tip];
  std::size_t first_index = 0;
  std::size_t second_index = 0;
  int pair_mrca_state = root_anchor_state;
  bool pair_mrca_missing = true;
  bool left_shared = false;
  int first_transition_count = -1;
  int second_transition_count = -1;
  int first_last_state = -1;
  int second_last_state = -1;
  double comparable_time = 0.0;

  for (int interval = 0; interval < end;) {
    const CompactLineageRun* first =
      interval < first_end ? &first_runs[first_index] : NULL;
    const CompactLineageRun* second =
      interval < second_end ? &second_runs[second_index] : NULL;
    int next = end;
    if (first != NULL) next = std::min(next, first->interval_after_last);
    if (second != NULL) next = std::min(next, second->interval_after_last);
    if (next <= interval) Rcpp::stop("compact lineage runs are unordered");
    const double duration = compact.cumulative_duration[next] -
      compact.cumulative_duration[interval];
    const bool shared = first != NULL && second != NULL &&
      first->phylogenetic_edge == second->phylogenetic_edge;

    if (shared) {
      if (left_shared) {
        Rcpp::stop("a tip pair cannot return to a shared phylogeny edge");
      }
      if (first->state != second->state) {
        Rcpp::stop("shared pair history has inconsistent states");
      }
      pair_mrca_state = first->state;
      pair_mrca_missing = false;
      synchronous[first->state * 4] += duration;
      workspace.first_shared[first->state] += duration;
      workspace.second_shared[second->state] += duration;
    } else {
      left_shared = left_shared || !pair_mrca_missing;
      if (pair_mrca_missing && !synthetic_root) {
        pair_mrca_state = root_anchor_state;
        pair_mrca_missing = false;
      }
      int first_count = -1;
      int second_count = -1;
      bool first_conserved = false;
      bool second_conserved = false;
      if (first != NULL) {
        if (first_transition_count < 0) {
          first_transition_count =
            (synthetic_root && pair_mrca_missing) ||
            pair_mrca_state == first->state ? 0 : 1;
        } else if (first_last_state != first->state) {
          ++first_transition_count;
        }
        first_last_state = first->state;
        first_count = first_transition_count;
        first_conserved = !synthetic_root && first_count == 0 &&
          first->state == pair_mrca_state;
        (first_conserved ?
          workspace.first_conserved : workspace.first_independent)
            [first->state] += duration;
        if (!first_conserved) {
          for (std::size_t depth = 0; depth < asynchronous.size(); ++depth) {
            if (data.maximum_transition_counts[depth] > 0) {
              const int actual = std::min(
                data.maximum_transition_counts[depth], first_count
              );
              add_compact_path_duration(
                workspace,
                depth,
                first->projected_path_ids[actual],
                actual,
                first->state,
                duration,
                true
              );
            }
          }
        }
      }
      if (second != NULL) {
        if (second_transition_count < 0) {
          second_transition_count =
            (synthetic_root && pair_mrca_missing) ||
            pair_mrca_state == second->state ? 0 : 1;
        } else if (second_last_state != second->state) {
          ++second_transition_count;
        }
        second_last_state = second->state;
        second_count = second_transition_count;
        second_conserved = !synthetic_root && second_count == 0 &&
          second->state == pair_mrca_state;
        (second_conserved ?
          workspace.second_conserved : workspace.second_independent)
            [second->state] += duration;
        if (!second_conserved) {
          for (std::size_t depth = 0; depth < asynchronous.size(); ++depth) {
            if (data.maximum_transition_counts[depth] > 0) {
              const int actual = std::min(
                data.maximum_transition_counts[depth], second_count
              );
              add_compact_path_duration(
                workspace,
                depth,
                second->projected_path_ids[actual],
                actual,
                second->state,
                duration,
                false
              );
            }
          }
        }
      }

      if (first != NULL && second != NULL) {
        const int state_count = data.state_levels.size();
        if (first_conserved && second_conserved) {
          synchronous[first->state * 4 + 1] += duration;
        } else if (first->state != second->state) {
          synchronous.back() += duration;
        } else if (first_conserved != second_conserved) {
          synchronous[first->state * 4 + 2] += duration;
        } else {
          int deepest_match = 0;
          for (std::size_t depth = 0;
               depth < data.maximum_transition_counts.size();
               ++depth) {
            const int requested = data.maximum_transition_counts[depth];
            const int first_actual = std::min(requested, first_count);
            const int second_actual = std::min(requested, second_count);
            if (first_actual > deepest_match &&
                first_actual == second_actual &&
                requested > 0 &&
                first->projected_path_ids[first_actual] ==
                  second->projected_path_ids[second_actual]) {
              deepest_match = first_actual;
            }
          }
          const int component = deepest_match == 0 ?
            first->state * 4 + 3 :
            state_count * 4 + (deepest_match - 1) * state_count + first->state;
          synchronous[component] += duration;
        }
      }
    }

    if (first != NULL && second != NULL) comparable_time += duration;
    interval = next;
    if (first != NULL && interval == first->interval_after_last) ++first_index;
    if (second != NULL && interval == second->interval_after_last) ++second_index;
  }

  for (std::size_t component = 0; component < synchronous.size(); ++component) {
    synchronous[component] /= comparable_time;
  }
  guard_similarity_composition(synchronous, data.time_tolerance);

  const double first_history = compact.cumulative_duration[first_end];
  const double second_history = compact.cumulative_duration[second_end];
  const bool bhattacharyya = metric == "bhattacharyya";
  const double denominator = bhattacharyya ?
    std::sqrt(first_history * second_history) :
    std::min(first_history, second_history);
  double baseline = 0.0;
  for (int state = 0; state < static_cast<int>(data.state_levels.size()); ++state) {
    workspace.shared_similarity[state] = compact_duration_similarity(
      workspace.first_shared[state], workspace.second_shared[state],
      denominator, bhattacharyya
    );
    workspace.conserved_similarity[state] = compact_duration_similarity(
      workspace.first_conserved[state], workspace.second_conserved[state],
      denominator, bhattacharyya
    );
    workspace.independent_similarity[state] = compact_duration_similarity(
      workspace.first_independent[state], workspace.second_independent[state],
      denominator, bhattacharyya
    );
    const double post = compact_duration_similarity(
      workspace.first_conserved[state] + workspace.first_independent[state],
      workspace.second_conserved[state] + workspace.second_independent[state],
      denominator, bhattacharyya
    );
    workspace.homoplasy_similarity[state] = std::max(
      0.0,
      post - workspace.conserved_similarity[state] -
        workspace.independent_similarity[state]
    );
    baseline += workspace.shared_similarity[state] + post;
  }
  const double divergent = std::max(0.0, 1.0 - baseline);
  for (std::size_t depth_index = 0; depth_index < asynchronous.size();
       ++depth_index) {
    std::vector<double>& values = asynchronous[depth_index];
    std::fill(values.begin(), values.end(), 0.0);
    const int depth = data.maximum_transition_counts[depth_index];
    const int block = depth == 0 ? 4 : 4 + depth;
    std::fill(
      workspace.parallel_similarity.begin(),
      workspace.parallel_similarity.end(),
      0.0
    );
    for (std::size_t touched = 0; depth > 0 &&
         touched < workspace.touched_paths[depth_index].size();
         ++touched) {
      const int path_id = workspace.touched_paths[depth_index][touched];
      const int actual =
        workspace.path_transition_counts[depth_index][path_id];
      if (actual == 0) continue;
      const double similarity = compact_duration_similarity(
        workspace.first_path_durations[depth_index][path_id],
        workspace.second_path_durations[depth_index][path_id],
        denominator,
        bhattacharyya
      );
      if (similarity == 0.0) continue;
      const int state = workspace.path_terminal_states[depth_index][path_id];
      values[state * block + 3 + actual - 1] += similarity;
      workspace.parallel_similarity[state] += similarity;
    }
    for (int state = 0; state < static_cast<int>(data.state_levels.size()); ++state) {
      const int start = state * block;
      values[start] = workspace.shared_similarity[state];
      values[start + 1] = workspace.conserved_similarity[state];
      values[start + 2] = workspace.homoplasy_similarity[state];
      values[start + 3 + depth] = depth == 0 ?
        workspace.independent_similarity[state] :
        std::max(
          0.0,
          workspace.independent_similarity[state] -
            workspace.parallel_similarity[state]
        );
    }
    values.back() = divergent;
    guard_similarity_composition(values, data.time_tolerance);
  }
}

void append_state(std::vector<std::string>& path, const std::string& state) {
  if (path.empty() || path.back() != state) path.push_back(state);
}

std::string format_path(const std::vector<std::string>& states) {
  std::string result;
  for (std::size_t index = 0; index < states.size(); ++index) {
    result += "|" + states[index];
  }
  return result;
}

std::string project_pair_path(
    const std::string& path,
    int maximum_transition_count) {
  return project_transition_path(
    path_components(path),
    maximum_transition_count
  );
}

PairCalculationResult calculate_one_pair_v21(
    const PreparedScenarioData& data,
    const ComponentDefinitions& definitions,
    int first_tip,
    int second_tip,
    const std::string& metric,
    bool complete) {
  PairCalculationResult result;
  result.first_tip_index = first_tip;
  result.second_tip_index = second_tip;
  result.total_comparable_time = 0;
  result.first_total_history_time = 0;
  result.second_total_history_time = 0;

  int last_shared_interval = -1;
  bool left_shared_history = false;
  for (int interval = 0; interval < data.interval_count; ++interval) {
    const std::size_t first_cell =
      scenario_cell_index(first_tip, interval, data.interval_count);
    const std::size_t second_cell =
      scenario_cell_index(second_tip, interval, data.interval_count);
    const bool first_active = data.lineage_is_active[first_cell];
    const bool second_active = data.lineage_is_active[second_cell];
    if (first_active) {
      result.first_total_history_time += data.interval_durations[interval];
    }
    if (second_active) {
      result.second_total_history_time += data.interval_durations[interval];
    }
    if (!first_active || !second_active) continue;
    result.total_comparable_time += data.interval_durations[interval];
    const bool shared =
      data.phylogenetic_edge_identifiers[first_cell] ==
      data.phylogenetic_edge_identifiers[second_cell];
    if (shared) {
      if (left_shared_history) {
        Rcpp::stop("a tip pair cannot return to a shared phylogeny edge");
      }
      last_shared_interval = interval;
    } else if (last_shared_interval >= 0) {
      left_shared_history = true;
    }
  }
  if (result.total_comparable_time <= data.time_tolerance) {
    Rcpp::stop("pair has zero total comparable time");
  }

  const bool synthetic_root = root_is_synthetic(data);
  TextValue pair_mrca{true, std::string()};
  if (last_shared_interval >= 0) {
    pair_mrca = data.states[scenario_cell_index(
      first_tip, last_shared_interval, data.interval_count
    )];
  } else if (!synthetic_root) {
    pair_mrca = TextValue{false, data.root_anchor};
  }

  std::vector<std::string> first_post_states;
  std::vector<std::string> second_post_states;
  for (int interval = 0; interval < data.interval_count; ++interval) {
    const std::size_t first_cell =
      scenario_cell_index(first_tip, interval, data.interval_count);
    const std::size_t second_cell =
      scenario_cell_index(second_tip, interval, data.interval_count);
    const bool first_active = data.lineage_is_active[first_cell];
    const bool second_active = data.lineage_is_active[second_cell];
    if (!first_active && !second_active) continue;

    const bool both_active = first_active && second_active;
    const bool shared = both_active &&
      data.phylogenetic_edge_identifiers[first_cell] ==
      data.phylogenetic_edge_identifiers[second_cell];
    const TextValue first_state = first_active ?
      data.states[first_cell] : TextValue{true, std::string()};
    const TextValue second_state = second_active ?
      data.states[second_cell] : TextValue{true, std::string()};
    TextValue first_post_path{true, std::string()};
    TextValue second_post_path{true, std::string()};
    int first_count = NA_INTEGER;
    int second_count = NA_INTEGER;
    bool first_conserved = false;
    bool second_conserved = false;
    std::string first_async = "Inactive";
    std::string second_async = "Inactive";

    if (shared) {
      first_async = "Shared";
      second_async = "Shared";
    } else {
      if (first_active) {
        if (first_post_states.empty()) {
          append_state(
            first_post_states,
            synthetic_root && last_shared_interval < 0 ?
              first_state.value : pair_mrca.value
          );
        }
        append_state(first_post_states, first_state.value);
        first_count = static_cast<int>(first_post_states.size()) - 1;
        first_post_path = TextValue{false, format_path(first_post_states)};
        first_conserved = !synthetic_root && first_count == 0 &&
          first_state.value == pair_mrca.value;
        first_async = first_conserved ? "Conserved" : "Independent";
      }
      if (second_active) {
        if (second_post_states.empty()) {
          append_state(
            second_post_states,
            synthetic_root && last_shared_interval < 0 ?
              second_state.value : pair_mrca.value
          );
        }
        append_state(second_post_states, second_state.value);
        second_count = static_cast<int>(second_post_states.size()) - 1;
        second_post_path = TextValue{false, format_path(second_post_states)};
        second_conserved = !synthetic_root && second_count == 0 &&
          second_state.value == pair_mrca.value;
        second_async = second_conserved ? "Conserved" : "Independent";
      }
    }

    std::string phase = "Asynchronous_only";
    TextValue sync_class{true, std::string()};
    int requested_depth = NA_INTEGER;
    int actual_depth = NA_INTEGER;
    TextValue matching_path{true, std::string()};
    if (shared) {
      phase = "Shared";
      sync_class = TextValue{false, "Shared"};
    } else if (both_active) {
      if (first_conserved && second_conserved) {
        phase = "Conserved";
        sync_class = TextValue{false, "Conserved"};
      } else {
        phase = "Independent";
        if (first_state.value != second_state.value) {
          sync_class = TextValue{false, "Divergent"};
        } else if (first_conserved != second_conserved) {
          phase = "Conserved_Homoplasy";
          sync_class = TextValue{false, "Conserved_Homoplasy"};
        } else {
          int deepest_match = 0;
          for (std::size_t depth_index = 0;
               depth_index < data.maximum_transition_counts.size();
               ++depth_index) {
            const int depth = data.maximum_transition_counts[depth_index];
            if (depth <= 0) continue;
            const int actual_match = std::min(
              depth, std::min(first_count, second_count)
            );
            const std::string first_projection =
              project_pair_path(first_post_path.value, depth);
            const std::string second_projection =
              project_pair_path(second_post_path.value, depth);
            if (actual_match > deepest_match &&
                first_projection == second_projection) {
              deepest_match = actual_match;
              requested_depth = depth;
              matching_path = TextValue{false, first_projection};
            }
          }
          if (deepest_match > 0) {
            actual_depth = deepest_match;
            sync_class = TextValue{false, "Parallel"};
          } else {
            requested_depth = NA_INTEGER;
            matching_path = TextValue{true, std::string()};
            sync_class = TextValue{false, "Convergent"};
          }
        }
      }
    }

    result.interval_classifications.push_back(IntervalClassificationRecord{
      interval + 1,
      data.interval_start_times[interval],
      data.interval_end_times[interval],
      data.interval_durations[interval],
      first_state,
      second_state,
      first_active,
      second_active,
      first_active ? data.full_paths[first_cell] : TextValue{true, std::string()},
      second_active ? data.full_paths[second_cell] : TextValue{true, std::string()},
      first_post_path,
      second_post_path,
      first_active ? data.phylogenetic_edge_identifiers[first_cell] : NA_INTEGER,
      second_active ? data.phylogenetic_edge_identifiers[second_cell] : NA_INTEGER,
      first_active ? data.scenario_edge_identifiers[first_cell] : NA_INTEGER,
      second_active ? data.scenario_edge_identifiers[second_cell] : NA_INTEGER,
      first_active ? data.scenario_edge_step_identifiers[first_cell] : NA_INTEGER,
      second_active ? data.scenario_edge_step_identifiers[second_cell] : NA_INTEGER,
      pair_mrca,
      synthetic_root,
      phase,
      sync_class,
      first_async,
      second_async,
      first_count,
      second_count,
      requested_depth,
      actual_depth,
      matching_path
    });
  }

  for (std::size_t record_index = 0;
       record_index < result.interval_classifications.size();
       ++record_index) {
    const IntervalClassificationRecord& record =
      result.interval_classifications[record_index];
    if (record.synchronous_class.is_missing) continue;
    std::string component = "Divergent";
    if (record.synchronous_class.value != "Divergent") {
      component = record.first_state.value + "_" +
        record.synchronous_class.value;
      if (record.synchronous_class.value == "Parallel") {
        component += "_" + transition_count_label(
          record.actual_transition_count_used
        );
        if (complete) {
          result.synchronous_parallel_times_by_path.push_back(
            SynchronousParallelTimeRecord{
              record.requested_maximum_transition_count,
              record.actual_transition_count_used,
              record.matching_path.value,
              record.first_state.value,
              record.duration
            }
          );
        }
      }
    }
    add_named_amount(
      result.synchronous_similarity,
      component,
      record.duration / result.total_comparable_time
    );
  }
  guard_similarity_composition(
    result.synchronous_similarity,
    data.time_tolerance
  );

  std::map<std::string, double> first_shared;
  std::map<std::string, double> second_shared;
  std::map<std::string, double> first_conserved;
  std::map<std::string, double> second_conserved;
  std::map<std::string, double> first_independent;
  std::map<std::string, double> second_independent;
  for (std::size_t record_index = 0;
       record_index < result.interval_classifications.size();
       ++record_index) {
    const IntervalClassificationRecord& record =
      result.interval_classifications[record_index];
    if (record.first_asynchronous_class == "Shared") {
      add_named_amount(first_shared, record.first_state.value, record.duration);
    } else if (record.first_asynchronous_class == "Conserved") {
      add_named_amount(first_conserved, record.first_state.value, record.duration);
    } else if (record.first_asynchronous_class == "Independent") {
      add_named_amount(first_independent, record.first_state.value, record.duration);
    }
    if (record.second_asynchronous_class == "Shared") {
      add_named_amount(second_shared, record.second_state.value, record.duration);
    } else if (record.second_asynchronous_class == "Conserved") {
      add_named_amount(second_conserved, record.second_state.value, record.duration);
    } else if (record.second_asynchronous_class == "Independent") {
      add_named_amount(second_independent, record.second_state.value, record.duration);
    }
  }

  std::map<std::string, double> shared_similarity;
  std::map<std::string, double> conserved_similarity;
  std::map<std::string, double> independent_similarity;
  std::map<std::string, double> conserved_homoplasy_similarity;
  double asynchronous_baseline = 0;
  for (std::size_t state_index = 0;
       state_index < data.state_levels.size();
       ++state_index) {
    const std::string& state = data.state_levels[state_index];
    const double shared = duration_similarity(
      named_amount(first_shared, state),
      named_amount(second_shared, state),
      result.first_total_history_time,
      result.second_total_history_time,
      metric
    );
    const double conserved = duration_similarity(
      named_amount(first_conserved, state),
      named_amount(second_conserved, state),
      result.first_total_history_time,
      result.second_total_history_time,
      metric
    );
    const double independent = duration_similarity(
      named_amount(first_independent, state),
      named_amount(second_independent, state),
      result.first_total_history_time,
      result.second_total_history_time,
      metric
    );
    const double post_separation = duration_similarity(
      named_amount(first_conserved, state) +
        named_amount(first_independent, state),
      named_amount(second_conserved, state) +
        named_amount(second_independent, state),
      result.first_total_history_time,
      result.second_total_history_time,
      metric
    );
    double conserved_homoplasy = post_separation - conserved - independent;
    if (conserved_homoplasy < -data.time_tolerance) {
      Rcpp::stop("conserved-homoplasy similarity became negative");
    }
    shared_similarity[state] = shared;
    conserved_similarity[state] = conserved;
    independent_similarity[state] = independent;
    conserved_homoplasy_similarity[state] =
      std::max(0.0, conserved_homoplasy);
    asynchronous_baseline += shared + post_separation;
  }
  const double divergent = std::max(0.0, 1.0 - asynchronous_baseline);

  result.asynchronous_similarity.resize(data.maximum_transition_counts.size());
  for (std::size_t depth_index = 0;
       depth_index < data.maximum_transition_counts.size();
       ++depth_index) {
    const int requested = data.maximum_transition_counts[depth_index];
    std::map<std::string, double>& view =
      result.asynchronous_similarity[depth_index];
    for (std::size_t state_index = 0;
         state_index < data.state_levels.size();
         ++state_index) {
      const std::string& state = data.state_levels[state_index];
      add_named_amount(view, state + "_Shared", shared_similarity[state]);
      add_named_amount(view, state + "_Conserved", conserved_similarity[state]);
      add_named_amount(
        view,
        state + "_Conserved_Homoplasy",
        conserved_homoplasy_similarity[state]
      );
    }
    if (requested == 0) {
      for (std::size_t state_index = 0;
           state_index < data.state_levels.size();
           ++state_index) {
        const std::string& state = data.state_levels[state_index];
        add_named_amount(
          view, state + "_Convergent", independent_similarity[state]
        );
      }
    } else {
      std::map<std::string, double> first_paths;
      std::map<std::string, double> second_paths;
      std::map<std::string, std::string> path_states;
      std::map<std::string, int> path_counts;
      std::vector<std::string> first_path_order;
      for (std::size_t record_index = 0;
           record_index < result.interval_classifications.size();
           ++record_index) {
        const IntervalClassificationRecord& record =
          result.interval_classifications[record_index];
        if (record.first_asynchronous_class == "Independent") {
          const std::string path = project_pair_path(
            record.first_post_separation_path.value, requested
          );
          if (first_paths.find(path) == first_paths.end()) {
            first_path_order.push_back(path);
          }
          add_named_amount(first_paths, path, record.duration);
          path_states[path] = record.first_state.value;
          path_counts[path] =
            static_cast<int>(path_components(path).size()) - 1;
          if (complete) {
            result.asynchronous_independent_times_by_path.push_back(
              AsynchronousIndependentPathTimeRecord{
                "i",
                data.tip_labels[first_tip],
                requested,
                record.first_available_transition_count,
                path_counts[path],
                path,
                record.first_state.value,
                record.duration,
                record.duration / result.first_total_history_time
              }
            );
          }
        }
        if (record.second_asynchronous_class == "Independent") {
          const std::string path = project_pair_path(
            record.second_post_separation_path.value, requested
          );
          add_named_amount(second_paths, path, record.duration);
          path_states[path] = record.second_state.value;
          path_counts[path] =
            static_cast<int>(path_components(path).size()) - 1;
          if (complete) {
            result.asynchronous_independent_times_by_path.push_back(
              AsynchronousIndependentPathTimeRecord{
                "j",
                data.tip_labels[second_tip],
                requested,
                record.second_available_transition_count,
                path_counts[path],
                path,
                record.second_state.value,
                record.duration,
                record.duration / result.second_total_history_time
              }
            );
          }
        }
      }
      std::map<std::string, double> parallel_by_state;
      for (std::size_t path_index = 0;
           path_index < first_path_order.size();
           ++path_index) {
        const std::string& path = first_path_order[path_index];
        if (second_paths.find(path) == second_paths.end() ||
            path_counts[path] == 0) {
          continue;
        }
        const double similarity = duration_similarity(
          first_paths[path],
          second_paths[path],
          result.first_total_history_time,
          result.second_total_history_time,
          metric
        );
        if (similarity == 0) continue;
        const std::string& state = path_states[path];
        add_named_amount(
          view,
          state + "_Parallel_" + transition_depth_partition_label(
            path_counts[path], requested
          ),
          similarity
        );
        add_named_amount(parallel_by_state, state, similarity);
        if (complete) {
          result.asynchronous_parallel_similarities_by_path.push_back(
            AsynchronousParallelSimilarityRecord{
              requested,
              path_counts[path],
              path,
              state,
              similarity,
              metric
            }
          );
        }
      }
      for (std::size_t state_index = 0;
           state_index < data.state_levels.size();
           ++state_index) {
        const std::string& state = data.state_levels[state_index];
        double convergence =
          independent_similarity[state] - named_amount(parallel_by_state, state);
        if (convergence < -data.time_tolerance) {
          Rcpp::stop("parallel path similarity exceeds same-state similarity");
        }
        add_named_amount(
          view, state + "_Convergent", std::max(0.0, convergence)
        );
      }
    }
    add_named_amount(view, "Divergent", divergent);
    guard_similarity_composition(view, data.time_tolerance);
  }

  static_cast<void>(definitions);
  return result;
}

Rcpp::DataFrame component_definitions_to_r(
    const ComponentDefinitions& component_definitions) {
  const int row_count = component_definitions.metadata_rows.size();
  Rcpp::CharacterVector matrix_names(row_count);
  Rcpp::CharacterVector states(row_count);
  Rcpp::CharacterVector paths(row_count);
  Rcpp::CharacterVector phases(row_count);
  Rcpp::CharacterVector kinds(row_count);
  Rcpp::IntegerVector requested_maximum_transition_counts(row_count);
  Rcpp::IntegerVector actual_transition_counts_used(row_count);
  Rcpp::CharacterVector views(row_count);

  for (int row_index = 0; row_index < row_count; ++row_index) {
    const ComponentMetadataRow& row =
      component_definitions.metadata_rows[row_index];
    matrix_names[row_index] = utf8_character(row.matrix_name);
    states[row_index] = row.state.is_missing ?
      NA_STRING : utf8_character(row.state.value);
    paths[row_index] = row.path.is_missing ?
      NA_STRING : utf8_character(row.path.value);
    phases[row_index] = row.phase;
    kinds[row_index] = row.kind;
    requested_maximum_transition_counts[row_index] =
      row.requested_maximum_transition_count;
    actual_transition_counts_used[row_index] = row.actual_transition_count_used;
    views[row_index] = row.view;
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("matrix_name") = matrix_names,
    Rcpp::Named("state") = states,
    Rcpp::Named("path") = paths,
    Rcpp::Named("phase") = phases,
    Rcpp::Named("kind") = kinds,
    Rcpp::Named("requested_maximum_transition_count") =
      requested_maximum_transition_counts,
    Rcpp::Named("actual_transition_count_used") =
      actual_transition_counts_used,
    Rcpp::Named("view") = views
  );
}

Rcpp::DataFrame interval_classifications_v21_to_r(
    const std::vector<PairCalculationResult>& pairs,
    const PreparedScenarioData& data) {
  int row_count = 0;
  for (std::size_t pair = 0; pair < pairs.size(); ++pair) {
    row_count += pairs[pair].interval_classifications.size();
  }
  Rcpp::IntegerVector interval_id(row_count), phylo_i(row_count),
    phylo_j(row_count), scenario_i(row_count), scenario_j(row_count),
    step_i(row_count), step_j(row_count), available_i(row_count),
    available_j(row_count), requested(row_count), actual(row_count),
    row_i(row_count), row_j(row_count);
  Rcpp::NumericVector start(row_count), end(row_count), duration(row_count);
  Rcpp::LogicalVector active_i(row_count), active_j(row_count),
    synthetic(row_count);
  Rcpp::CharacterVector state_i(row_count), state_j(row_count),
    path_i(row_count), path_j(row_count), post_path_i(row_count),
    post_path_j(row_count), mrca(row_count), phase(row_count),
    sync_class(row_count), async_i(row_count), async_j(row_count),
    matching_path(row_count), tip_i(row_count), tip_j(row_count);
  Rcpp::List post_states_i(row_count), post_states_j(row_count);

  int output = 0;
  for (std::size_t pair_index = 0; pair_index < pairs.size(); ++pair_index) {
    const PairCalculationResult& pair = pairs[pair_index];
    std::vector<std::string> latest_i;
    std::vector<std::string> latest_j;
    for (std::size_t record_index = 0;
         record_index < pair.interval_classifications.size();
         ++record_index, ++output) {
      const IntervalClassificationRecord& record =
        pair.interval_classifications[record_index];
      if (!record.first_post_separation_path.is_missing) {
        latest_i = path_components(record.first_post_separation_path.value);
      }
      if (!record.second_post_separation_path.is_missing) {
        latest_j = path_components(record.second_post_separation_path.value);
      }
      interval_id[output] = record.interval_identifier;
      start[output] = record.start_time;
      end[output] = record.end_time;
      duration[output] = record.duration;
      state_i[output] = record.first_state.is_missing ?
        NA_STRING : utf8_character(record.first_state.value);
      state_j[output] = record.second_state.is_missing ?
        NA_STRING : utf8_character(record.second_state.value);
      active_i[output] = record.first_lineage_active;
      active_j[output] = record.second_lineage_active;
      path_i[output] = record.first_path.is_missing ?
        NA_STRING : utf8_character(record.first_path.value);
      path_j[output] = record.second_path.is_missing ?
        NA_STRING : utf8_character(record.second_path.value);
      post_path_i[output] = record.first_post_separation_path.is_missing ?
        NA_STRING : utf8_character(record.first_post_separation_path.value);
      post_path_j[output] = record.second_post_separation_path.is_missing ?
        NA_STRING : utf8_character(record.second_post_separation_path.value);
      phylo_i[output] = record.first_phylogenetic_edge_identifier;
      phylo_j[output] = record.second_phylogenetic_edge_identifier;
      scenario_i[output] = record.first_scenario_edge_identifier;
      scenario_j[output] = record.second_scenario_edge_identifier;
      step_i[output] = record.first_scenario_step_identifier;
      step_j[output] = record.second_scenario_step_identifier;
      mrca[output] = record.pair_mrca_state.is_missing ?
        NA_STRING : utf8_character(record.pair_mrca_state.value);
      synthetic[output] = record.root_is_synthetic;
      phase[output] = record.phase;
      sync_class[output] = record.synchronous_class.is_missing ?
        NA_STRING : utf8_character(record.synchronous_class.value);
      async_i[output] = record.first_asynchronous_class;
      async_j[output] = record.second_asynchronous_class;
      available_i[output] = record.first_available_transition_count;
      available_j[output] = record.second_available_transition_count;
      requested[output] = record.requested_maximum_transition_count;
      actual[output] = record.actual_transition_count_used;
      matching_path[output] = record.matching_path.is_missing ?
        NA_STRING : utf8_character(record.matching_path.value);
      post_states_i[output] = utf8_character_vector(latest_i);
      post_states_j[output] = utf8_character_vector(latest_j);
      row_i[output] = pair.first_tip_index + 1;
      row_j[output] = pair.second_tip_index + 1;
      tip_i[output] = source_tip_label(data, pair.first_tip_index);
      tip_j[output] = source_tip_label(data, pair.second_tip_index);
    }
  }
  post_states_i.attr("class") = "AsIs";
  post_states_j.attr("class") = "AsIs";

  return Rcpp::DataFrame::create(
    Rcpp::Named("interval_id") = interval_id,
    Rcpp::Named("start") = start,
    Rcpp::Named("end") = end,
    Rcpp::Named("duration") = duration,
    Rcpp::Named("state_i") = state_i,
    Rcpp::Named("state_j") = state_j,
    Rcpp::Named("first_lineage_active") = active_i,
    Rcpp::Named("second_lineage_active") = active_j,
    Rcpp::Named("path_i") = path_i,
    Rcpp::Named("path_j") = path_j,
    Rcpp::Named("post_separation_path_i") = post_path_i,
    Rcpp::Named("post_separation_path_j") = post_path_j,
    Rcpp::Named("phylo_edge_i") = phylo_i,
    Rcpp::Named("phylo_edge_j") = phylo_j,
    Rcpp::Named("scenario_edge_i") = scenario_i,
    Rcpp::Named("scenario_edge_j") = scenario_j,
    Rcpp::Named("scenario_step_i") = step_i,
    Rcpp::Named("scenario_step_j") = step_j,
    Rcpp::Named("pair_mrca_state") = mrca,
    Rcpp::Named("root_is_synthetic") = synthetic,
    Rcpp::Named("phase") = phase,
    Rcpp::Named("sync_class") = sync_class,
    Rcpp::Named("asynchronous_class_i") = async_i,
    Rcpp::Named("asynchronous_class_j") = async_j,
    Rcpp::Named("available_transition_count_i") = available_i,
    Rcpp::Named("available_transition_count_j") = available_j,
    Rcpp::Named("requested_maximum_transition_count") = requested,
    Rcpp::Named("actual_transition_count_used") = actual,
    Rcpp::Named("matching_path") = matching_path,
    Rcpp::Named("post_separation_states_i") = post_states_i,
    Rcpp::Named("post_separation_states_j") = post_states_j,
    Rcpp::Named("row_i") = row_i,
    Rcpp::Named("row_j") = row_j,
    Rcpp::Named("tip_i") = tip_i,
    Rcpp::Named("tip_j") = tip_j
  );
}

Rcpp::DataFrame synchronous_parallel_times_to_r(
    const std::vector<PairCalculationResult>& pair_results,
    const PreparedScenarioData& prepared_scenario_data) {
  int row_count = 0;
  for (std::size_t pair_index = 0;
       pair_index < pair_results.size();
       ++pair_index) {
    row_count += pair_results[pair_index]
      .synchronous_parallel_times_by_path.size();
  }
  if (row_count == 0) return Rcpp::DataFrame::create();

  Rcpp::IntegerVector requested_maximum_transition_counts(row_count);
  Rcpp::IntegerVector actual_transition_counts_used(row_count);
  Rcpp::CharacterVector paths(row_count);
  Rcpp::CharacterVector terminal_states(row_count);
  Rcpp::NumericVector durations(row_count);
  Rcpp::IntegerVector first_tip_indices(row_count);
  Rcpp::IntegerVector second_tip_indices(row_count);
  Rcpp::CharacterVector first_tip_labels(row_count);
  Rcpp::CharacterVector second_tip_labels(row_count);

  int output_row_index = 0;
  for (std::size_t pair_index = 0;
       pair_index < pair_results.size();
       ++pair_index) {
    const PairCalculationResult& pair_result = pair_results[pair_index];
    for (std::size_t record_index = 0;
         record_index < pair_result.synchronous_parallel_times_by_path.size();
         ++record_index, ++output_row_index) {
      const SynchronousParallelTimeRecord& record =
        pair_result.synchronous_parallel_times_by_path[record_index];
      requested_maximum_transition_counts[output_row_index] =
        record.requested_maximum_transition_count;
      actual_transition_counts_used[output_row_index] =
        record.actual_transition_count_used;
      paths[output_row_index] = utf8_character(record.path);
      terminal_states[output_row_index] = utf8_character(record.terminal_state);
      durations[output_row_index] = record.duration;
      first_tip_indices[output_row_index] = pair_result.first_tip_index + 1;
      second_tip_indices[output_row_index] = pair_result.second_tip_index + 1;
      first_tip_labels[output_row_index] =
        source_tip_label(prepared_scenario_data, pair_result.first_tip_index);
      second_tip_labels[output_row_index] =
        source_tip_label(prepared_scenario_data, pair_result.second_tip_index);
    }
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("requested_maximum_transition_count") =
      requested_maximum_transition_counts,
    Rcpp::Named("actual_transition_count_used") =
      actual_transition_counts_used,
    Rcpp::Named("path") = paths,
    Rcpp::Named("terminal_state") = terminal_states,
    Rcpp::Named("duration") = durations,
    Rcpp::Named("row_i") = first_tip_indices,
    Rcpp::Named("row_j") = second_tip_indices,
    Rcpp::Named("tip_i") = first_tip_labels,
    Rcpp::Named("tip_j") = second_tip_labels
  );
}

Rcpp::DataFrame asynchronous_independent_times_to_r(
    const std::vector<PairCalculationResult>& pair_results,
    const PreparedScenarioData& prepared_scenario_data) {
  int row_count = 0;
  for (std::size_t pair_index = 0;
       pair_index < pair_results.size();
       ++pair_index) {
    row_count += pair_results[pair_index]
      .asynchronous_independent_times_by_path.size();
  }
  if (row_count == 0) return Rcpp::DataFrame::create();

  Rcpp::CharacterVector lineages(row_count);
  Rcpp::CharacterVector tips(row_count);
  Rcpp::IntegerVector requested_maximum_transition_counts(row_count);
  Rcpp::IntegerVector available_transition_counts(row_count);
  Rcpp::IntegerVector actual_transition_counts_used(row_count);
  Rcpp::CharacterVector paths(row_count);
  Rcpp::CharacterVector terminal_states(row_count);
  Rcpp::NumericVector durations(row_count);
  Rcpp::NumericVector proportions(row_count);
  Rcpp::IntegerVector first_tip_indices(row_count);
  Rcpp::IntegerVector second_tip_indices(row_count);
  Rcpp::CharacterVector first_tip_labels(row_count);
  Rcpp::CharacterVector second_tip_labels(row_count);

  int output_row_index = 0;
  for (std::size_t pair_index = 0;
       pair_index < pair_results.size();
       ++pair_index) {
    const PairCalculationResult& pair_result = pair_results[pair_index];
    for (std::size_t record_index = 0;
         record_index <
           pair_result.asynchronous_independent_times_by_path.size();
         ++record_index, ++output_row_index) {
      const AsynchronousIndependentPathTimeRecord& record =
        pair_result.asynchronous_independent_times_by_path[record_index];
      lineages[output_row_index] = record.lineage;
      tips[output_row_index] = source_tip_label(
        prepared_scenario_data,
        record.lineage == "i" ?
          pair_result.first_tip_index : pair_result.second_tip_index
      );
      requested_maximum_transition_counts[output_row_index] =
        record.requested_maximum_transition_count;
      available_transition_counts[output_row_index] =
        record.available_transition_count;
      actual_transition_counts_used[output_row_index] =
        record.actual_transition_count_used;
      paths[output_row_index] = utf8_character(record.path);
      terminal_states[output_row_index] = utf8_character(record.terminal_state);
      durations[output_row_index] = record.duration;
      proportions[output_row_index] = record.proportion;
      first_tip_indices[output_row_index] = pair_result.first_tip_index + 1;
      second_tip_indices[output_row_index] = pair_result.second_tip_index + 1;
      first_tip_labels[output_row_index] =
        source_tip_label(prepared_scenario_data, pair_result.first_tip_index);
      second_tip_labels[output_row_index] =
        source_tip_label(prepared_scenario_data, pair_result.second_tip_index);
    }
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("lineage") = lineages,
    Rcpp::Named("tip") = tips,
    Rcpp::Named("requested_maximum_transition_count") =
      requested_maximum_transition_counts,
    Rcpp::Named("available_transition_count") = available_transition_counts,
    Rcpp::Named("actual_transition_count_used") =
      actual_transition_counts_used,
    Rcpp::Named("path") = paths,
    Rcpp::Named("terminal_state") = terminal_states,
    Rcpp::Named("duration") = durations,
    Rcpp::Named("proportion") = proportions,
    Rcpp::Named("row_i") = first_tip_indices,
    Rcpp::Named("row_j") = second_tip_indices,
    Rcpp::Named("tip_i") = first_tip_labels,
    Rcpp::Named("tip_j") = second_tip_labels
  );
}

Rcpp::DataFrame asynchronous_parallel_similarities_to_r(
    const std::vector<PairCalculationResult>& pair_results,
    const PreparedScenarioData& prepared_scenario_data) {
  int row_count = 0;
  for (std::size_t pair_index = 0;
       pair_index < pair_results.size();
       ++pair_index) {
    row_count += pair_results[pair_index]
      .asynchronous_parallel_similarities_by_path.size();
  }
  if (row_count == 0) return Rcpp::DataFrame::create();

  Rcpp::IntegerVector requested_maximum_transition_counts(row_count);
  Rcpp::IntegerVector actual_transition_counts_used(row_count);
  Rcpp::CharacterVector paths(row_count);
  Rcpp::CharacterVector terminal_states(row_count);
  Rcpp::NumericVector similarities(row_count);
  Rcpp::CharacterVector metrics(row_count);
  Rcpp::IntegerVector first_tip_indices(row_count);
  Rcpp::IntegerVector second_tip_indices(row_count);
  Rcpp::CharacterVector first_tip_labels(row_count);
  Rcpp::CharacterVector second_tip_labels(row_count);

  int output_row_index = 0;
  for (std::size_t pair_index = 0;
       pair_index < pair_results.size();
       ++pair_index) {
    const PairCalculationResult& pair_result = pair_results[pair_index];
    for (std::size_t record_index = 0;
         record_index <
           pair_result.asynchronous_parallel_similarities_by_path.size();
         ++record_index, ++output_row_index) {
      const AsynchronousParallelSimilarityRecord& record =
        pair_result.asynchronous_parallel_similarities_by_path[record_index];
      requested_maximum_transition_counts[output_row_index] =
        record.requested_maximum_transition_count;
      actual_transition_counts_used[output_row_index] =
        record.actual_transition_count_used;
      paths[output_row_index] = utf8_character(record.path);
      terminal_states[output_row_index] = utf8_character(record.terminal_state);
      similarities[output_row_index] = record.similarity;
      metrics[output_row_index] = record.metric;
      first_tip_indices[output_row_index] = pair_result.first_tip_index + 1;
      second_tip_indices[output_row_index] = pair_result.second_tip_index + 1;
      first_tip_labels[output_row_index] =
        source_tip_label(prepared_scenario_data, pair_result.first_tip_index);
      second_tip_labels[output_row_index] =
        source_tip_label(prepared_scenario_data, pair_result.second_tip_index);
    }
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("requested_maximum_transition_count") =
      requested_maximum_transition_counts,
    Rcpp::Named("actual_transition_count_used") =
      actual_transition_counts_used,
    Rcpp::Named("path") = paths,
    Rcpp::Named("terminal_state") = terminal_states,
    Rcpp::Named("similarity") = similarities,
    Rcpp::Named("metric") = metrics,
    Rcpp::Named("row_i") = first_tip_indices,
    Rcpp::Named("row_j") = second_tip_indices,
    Rcpp::Named("tip_i") = first_tip_labels,
    Rcpp::Named("tip_j") = second_tip_labels
  );
}

Rcpp::NumericMatrix total_comparable_time_matrix_to_r(
    const std::vector<PairCalculationResult>& pair_results,
    const PreparedScenarioData& prepared_scenario_data) {
  Rcpp::NumericMatrix total_comparable_time_by_pair(
    prepared_scenario_data.tip_count,
    prepared_scenario_data.tip_count
  );
  set_matrix_names(
    total_comparable_time_by_pair,
    prepared_scenario_data.tip_label_source,
    prepared_scenario_data.tip_label_source
  );
  for (std::size_t pair_index = 0;
       pair_index < pair_results.size();
       ++pair_index) {
    const PairCalculationResult& pair_result = pair_results[pair_index];
    total_comparable_time_by_pair(
      pair_result.first_tip_index,
      pair_result.second_tip_index
    ) = pair_result.total_comparable_time;
    total_comparable_time_by_pair(
      pair_result.second_tip_index,
      pair_result.first_tip_index
    ) = pair_result.total_comparable_time;
  }
  return total_comparable_time_by_pair;
}

Rcpp::NumericMatrix asynchronous_available_similarity_to_r(
    const std::vector<PairCalculationResult>& pairs,
    const PreparedScenarioData& data,
    const std::string& metric) {
  Rcpp::NumericMatrix result(data.tip_count, data.tip_count);
  set_matrix_names(result, data.tip_label_source, data.tip_label_source);
  for (std::size_t index = 0; index < pairs.size(); ++index) {
    const PairCalculationResult& pair = pairs[index];
    const double available = metric == "bhattacharyya" ?
      std::sqrt(
        pair.first_total_history_time * pair.second_total_history_time
      ) :
      std::min(
        pair.first_total_history_time, pair.second_total_history_time
      );
    result(pair.first_tip_index, pair.second_tip_index) = available;
    result(pair.second_tip_index, pair.first_tip_index) = available;
  }
  return result;
}

Rcpp::NumericVector tip_history_times_to_r(const PreparedScenarioData& data) {
  Rcpp::NumericVector result(data.tip_count);
  for (int tip = 0; tip < data.tip_count; ++tip) {
    double total = 0;
    for (int interval = 0; interval < data.interval_count; ++interval) {
      if (data.lineage_is_active[
            scenario_cell_index(tip, interval, data.interval_count)
          ]) {
        total += data.interval_durations[interval];
      }
    }
    result[tip] = total;
  }
  result.attr("names") = data.tip_label_source;
  return result;
}

Rcpp::NumericVector tip_history_times_to_r(
    const PreparedScenarioData& data,
    const CompactScenarioData& compact) {
  Rcpp::NumericVector result(data.tip_count);
  for (int tip = 0; tip < data.tip_count; ++tip) {
    result[tip] =
      compact.cumulative_duration[data.tip_end_columns[tip] - 1];
  }
  result.attr("names") = data.tip_label_source;
  return result;
}

CompactMatrixStorage initialize_compact_matrix_storage(
    const PreparedScenarioData& data,
    const ComponentDefinitions& definitions,
    bool record_matrix_output) {
  CompactMatrixStorage storage;
  storage.synchronous_totals.assign(
    definitions.synchronous_component_names.size(), 0.0
  );
  if (record_matrix_output) {
    for (std::size_t component = 0;
         component < definitions.synchronous_component_names.size();
         ++component) {
      Rcpp::NumericMatrix matrix(data.tip_count, data.tip_count);
      set_matrix_names(matrix, data.tip_label_source, data.tip_label_source);
      storage.synchronous.push_back(matrix);
    }
  }
  storage.asynchronous.resize(data.maximum_transition_counts.size());
  storage.asynchronous_totals.resize(data.maximum_transition_counts.size());
  for (std::size_t depth = 0;
       depth < data.maximum_transition_counts.size();
       ++depth) {
    storage.asynchronous_totals[depth].assign(
      definitions.asynchronous_component_names[depth].size(), 0.0
    );
    if (record_matrix_output) {
      for (std::size_t component = 0;
           component < definitions.asynchronous_component_names[depth].size();
           ++component) {
        Rcpp::NumericMatrix matrix(data.tip_count, data.tip_count);
        set_matrix_names(matrix, data.tip_label_source, data.tip_label_source);
        storage.asynchronous[depth].push_back(matrix);
      }
    }
  }
  if (record_matrix_output) {
    storage.synchronous_available =
      Rcpp::NumericMatrix(data.tip_count, data.tip_count);
    storage.asynchronous_available =
      Rcpp::NumericMatrix(data.tip_count, data.tip_count);
    set_matrix_names(
      storage.synchronous_available,
      data.tip_label_source,
      data.tip_label_source
    );
    set_matrix_names(
      storage.asynchronous_available,
      data.tip_label_source,
      data.tip_label_source
    );
  }
  storage.synchronous_available_total = 0.0;
  storage.asynchronous_available_total = 0.0;
  return storage;
}

void write_compact_pair(
    std::vector<Rcpp::NumericMatrix>& matrices,
    std::vector<double>& totals,
    const std::vector<double>& values,
    int first_tip,
    int second_tip,
    double weight) {
  for (std::size_t component = 0; component < values.size(); ++component) {
    if (values[component] == 0.0) continue;
    if (!matrices.empty()) {
      matrices[component](first_tip, second_tip) = values[component];
      matrices[component](second_tip, first_tip) = values[component];
    }
    totals[component] += values[component] * weight;
  }
}

CompactMatrixStorage calculate_compact_matrix_storage(
    const PreparedScenarioData& data,
    const ComponentDefinitions& definitions,
    const CompactScenarioData& compact,
    const std::string& metric,
    bool record_matrix_output) {
  CompactMatrixStorage storage =
    initialize_compact_matrix_storage(data, definitions, record_matrix_output);
  std::vector<double> synchronous(
    definitions.synchronous_component_names.size()
  );
  std::vector<std::vector<double> > asynchronous(
    data.maximum_transition_counts.size()
  );
  for (std::size_t depth = 0; depth < asynchronous.size(); ++depth) {
    asynchronous[depth].resize(
      definitions.asynchronous_component_names[depth].size()
    );
  }
  const int state_count = data.state_levels.size();
  const bool synthetic_root = root_is_synthetic(data);
  const std::vector<std::string>::const_iterator root_position = std::find(
    data.state_levels.begin(), data.state_levels.end(), data.root_anchor
  );
  const int root_anchor_state = root_position == data.state_levels.end() ?
    state_count : root_position - data.state_levels.begin();
  CompactWorkspace workspace;
  workspace.first_shared.resize(state_count);
  workspace.second_shared.resize(state_count);
  workspace.first_conserved.resize(state_count);
  workspace.second_conserved.resize(state_count);
  workspace.first_independent.resize(state_count);
  workspace.second_independent.resize(state_count);
  workspace.shared_similarity.resize(state_count);
  workspace.conserved_similarity.resize(state_count);
  workspace.independent_similarity.resize(state_count);
  workspace.homoplasy_similarity.resize(state_count);
  workspace.parallel_similarity.resize(state_count);
  workspace.first_path_durations.resize(asynchronous.size());
  workspace.second_path_durations.resize(asynchronous.size());
  workspace.path_stamps.resize(asynchronous.size());
  workspace.touched_paths.resize(asynchronous.size());
  workspace.path_transition_counts.resize(asynchronous.size());
  workspace.path_terminal_states.resize(asynchronous.size());
  for (std::size_t depth = 0; depth < asynchronous.size(); ++depth) {
    workspace.first_path_durations[depth].resize(compact.path_count);
    workspace.second_path_durations[depth].resize(compact.path_count);
    workspace.path_stamps[depth].resize(compact.path_count);
    workspace.path_transition_counts[depth].resize(compact.path_count);
    workspace.path_terminal_states[depth].resize(compact.path_count);
  }
  workspace.pair_stamp = 0;
  double shortest_history = std::numeric_limits<double>::infinity();
  double longest_history = 0.0;
  for (int tip = 0; tip < data.tip_count; ++tip) {
    const double history =
      compact.cumulative_duration[data.tip_end_columns[tip] - 1];
    shortest_history = std::min(shortest_history, history);
    longest_history = std::max(longest_history, history);
  }
  const bool weighted =
    longest_history - shortest_history > data.time_tolerance;

  for (int first_tip = 0; first_tip < data.tip_count - 1; ++first_tip) {
    const double first_history =
      compact.cumulative_duration[data.tip_end_columns[first_tip] - 1];
    for (int second_tip = first_tip + 1;
         second_tip < data.tip_count;
         ++second_tip) {
      const double second_history =
        compact.cumulative_duration[data.tip_end_columns[second_tip] - 1];
      calculate_compact_pair_similarities(
        data,
        compact,
        first_tip,
        second_tip,
        synthetic_root,
        root_anchor_state,
        metric,
        workspace,
        synchronous,
        asynchronous
      );
      const int comparable_end = std::min(
        data.tip_end_columns[first_tip], data.tip_end_columns[second_tip]
      ) - 1;
      const double comparable =
        compact.cumulative_duration[comparable_end];
      const double available = metric == "bhattacharyya" ?
        std::sqrt(first_history * second_history) :
        std::min(first_history, second_history);
      write_compact_pair(
        storage.synchronous,
        storage.synchronous_totals,
        synchronous,
        first_tip,
        second_tip,
        weighted ? comparable : 1.0
      );
      for (std::size_t depth = 0; depth < asynchronous.size(); ++depth) {
        write_compact_pair(
          storage.asynchronous[depth],
          storage.asynchronous_totals[depth],
          asynchronous[depth],
          first_tip,
          second_tip,
          weighted ? available : 1.0
        );
      }

      if (record_matrix_output) {
        storage.synchronous_available(first_tip, second_tip) = comparable;
        storage.synchronous_available(second_tip, first_tip) = comparable;
        storage.asynchronous_available(first_tip, second_tip) = available;
        storage.asynchronous_available(second_tip, first_tip) = available;
      }
      if (weighted) {
        storage.synchronous_available_total += comparable;
        storage.asynchronous_available_total += available;
      }
    }
  }
  return storage;
}

Rcpp::List compact_matrix_storage_to_r(
    const CompactMatrixStorage& storage,
    const PreparedScenarioData& data,
    const ComponentDefinitions& definitions) {
  Rcpp::List synchronous(storage.synchronous.size());
  Rcpp::CharacterVector synchronous_names(storage.synchronous.size());
  for (std::size_t component = 0;
       component < storage.synchronous.size();
       ++component) {
    synchronous[component] = storage.synchronous[component];
    synchronous_names[component] = utf8_character(
      definitions.synchronous_component_names[component]
    );
  }
  synchronous.attr("names") = synchronous_names;

  Rcpp::List asynchronous(storage.asynchronous.size());
  Rcpp::CharacterVector view_names(storage.asynchronous.size());
  for (std::size_t depth = 0; depth < storage.asynchronous.size(); ++depth) {
    Rcpp::List view(storage.asynchronous[depth].size());
    Rcpp::CharacterVector component_names(storage.asynchronous[depth].size());
    for (std::size_t component = 0;
         component < storage.asynchronous[depth].size();
         ++component) {
      view[component] = storage.asynchronous[depth][component];
      component_names[component] = utf8_character(
        definitions.asynchronous_component_names[depth][component]
      );
    }
    view.attr("names") = component_names;
    asynchronous[depth] = view;
    view_names[depth] = transition_view_name(
      data.maximum_transition_counts[depth]
    );
  }
  asynchronous.attr("names") = view_names;
  return Rcpp::List::create(
    Rcpp::Named("sync") = synchronous,
    Rcpp::Named("async") = asynchronous
  );
}

Rcpp::List compact_summary_family_to_r(
    const std::vector<double>& totals,
    const std::vector<std::string>& names,
    bool weighted,
    double available_total) {
  Rcpp::NumericVector total_values = Rcpp::wrap(totals);
  total_values.attr("names") = utf8_character_vector(names);
  const double total = Rcpp::sum(total_values);
  Rcpp::NumericVector proportions(total_values.size());
  if (total > 0.0) proportions = total_values / total;
  proportions.attr("names") = total_values.attr("names");
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("totals") = total_values,
    Rcpp::Named("tree_wide_proportions") = proportions,
    Rcpp::Named("means") = proportions,
    Rcpp::Named("total_across_matrices") = total
  );
  if (weighted) result["available_similarity_total"] = available_total;
  return result;
}

Rcpp::List compact_summaries_to_r(
    const CompactMatrixStorage& storage,
    const PreparedScenarioData& data,
    const ComponentDefinitions& definitions,
    const Rcpp::DataFrame& component_metadata,
    bool tree_is_ultrametric,
    const Rcpp::NumericVector& tip_history_times,
    bool include_availability_matrices) {
  const bool weighted = !tree_is_ultrametric;
  Rcpp::List asynchronous(storage.asynchronous_totals.size());
  Rcpp::CharacterVector view_names(storage.asynchronous_totals.size());
  for (std::size_t depth = 0;
       depth < storage.asynchronous_totals.size();
       ++depth) {
    asynchronous[depth] = compact_summary_family_to_r(
      storage.asynchronous_totals[depth],
      definitions.asynchronous_component_names[depth],
      weighted,
      storage.asynchronous_available_total
    );
    view_names[depth] = transition_view_name(
      data.maximum_transition_counts[depth]
    );
  }
  asynchronous.attr("names") = view_names;
  Rcpp::List weighting = Rcpp::List::create(
    Rcpp::Named("tree_is_ultrametric") = tree_is_ultrametric,
    Rcpp::Named("available_similarity_weighting_applied") = weighted,
    Rcpp::Named("tip_history_times") = tip_history_times
  );
  if (include_availability_matrices) {
    weighting["synchronous_available_similarity_by_pair"] =
      storage.synchronous_available;
    weighting["asynchronous_available_similarity_by_pair"] =
      storage.asynchronous_available;
  }
  return Rcpp::List::create(
    Rcpp::Named("sync") = compact_summary_family_to_r(
      storage.synchronous_totals,
      definitions.synchronous_component_names,
      weighted,
      storage.synchronous_available_total
    ),
    Rcpp::Named("async") = asynchronous,
    Rcpp::Named("components") = component_metadata,
    Rcpp::Named("weighting") = weighting
  );
}

Rcpp::List matrix_families_to_r(
    const std::vector<PairCalculationResult>& pair_results,
    const PreparedScenarioData& prepared_scenario_data,
    const ComponentDefinitions& component_definitions) {
  Rcpp::List synchronous_matrices(
    component_definitions.synchronous_component_names.size()
  );
  Rcpp::CharacterVector synchronous_matrix_names(
    component_definitions.synchronous_component_names.size()
  );
  std::map<std::string, int> synchronous_matrix_indices;
  for (std::size_t component_index = 0;
       component_index < component_definitions.synchronous_component_names.size();
       ++component_index) {
    const std::string& component_name =
      component_definitions.synchronous_component_names[component_index];
    Rcpp::NumericMatrix component_matrix(
      prepared_scenario_data.tip_count,
      prepared_scenario_data.tip_count
    );
    set_matrix_names(
      component_matrix,
      prepared_scenario_data.tip_label_source,
      prepared_scenario_data.tip_label_source
    );
    synchronous_matrices[component_index] = component_matrix;
    synchronous_matrix_names[component_index] = utf8_character(component_name);
    synchronous_matrix_indices[component_name] = component_index;
  }
  synchronous_matrices.attr("names") = synchronous_matrix_names;

  Rcpp::List asynchronous_views(
    prepared_scenario_data.maximum_transition_counts.size()
  );
  Rcpp::CharacterVector asynchronous_view_names(
    prepared_scenario_data.maximum_transition_counts.size()
  );
  std::vector<std::map<std::string, int> > asynchronous_matrix_indices(
    prepared_scenario_data.maximum_transition_counts.size()
  );
  for (std::size_t depth_index = 0;
       depth_index < prepared_scenario_data.maximum_transition_counts.size();
       ++depth_index) {
    const std::vector<std::string>& component_names =
      component_definitions.asynchronous_component_names[depth_index];
    Rcpp::List view_matrices(component_names.size());
    Rcpp::CharacterVector view_matrix_names(component_names.size());
    for (std::size_t component_index = 0;
         component_index < component_names.size();
         ++component_index) {
      Rcpp::NumericMatrix component_matrix(
        prepared_scenario_data.tip_count,
        prepared_scenario_data.tip_count
      );
      set_matrix_names(
        component_matrix,
        prepared_scenario_data.tip_label_source,
        prepared_scenario_data.tip_label_source
      );
      view_matrices[component_index] = component_matrix;
      view_matrix_names[component_index] = utf8_character(
        component_names[component_index]
      );
      asynchronous_matrix_indices[depth_index][component_names[component_index]] =
        component_index;
    }
    view_matrices.attr("names") = view_matrix_names;
    asynchronous_views[depth_index] = view_matrices;
    asynchronous_view_names[depth_index] = transition_view_name(
      prepared_scenario_data.maximum_transition_counts[depth_index]
    );
  }
  asynchronous_views.attr("names") = asynchronous_view_names;

  for (std::size_t pair_index = 0;
       pair_index < pair_results.size();
       ++pair_index) {
    const PairCalculationResult& pair_result = pair_results[pair_index];
    for (std::map<std::string, double>::const_iterator position =
           pair_result.synchronous_similarity.begin();
         position != pair_result.synchronous_similarity.end();
         ++position) {
      const std::map<std::string, int>::const_iterator matrix_position =
        synchronous_matrix_indices.find(position->first);
      if (matrix_position == synchronous_matrix_indices.end()) {
        std::ostringstream message;
        message << "synchronous component '" << position->first
                << "' is outside the declared domain:";
        for (std::map<std::string, int>::const_iterator declared =
               synchronous_matrix_indices.begin();
             declared != synchronous_matrix_indices.end();
             ++declared) {
          message << " '" << declared->first << "'";
        }
        Rcpp::stop(message.str());
      }
      Rcpp::NumericMatrix component_matrix =
        synchronous_matrices[matrix_position->second];
      component_matrix(
        pair_result.first_tip_index,
        pair_result.second_tip_index
      ) = position->second;
      component_matrix(
        pair_result.second_tip_index,
        pair_result.first_tip_index
      ) = position->second;
    }

    for (std::size_t depth_index = 0;
         depth_index < pair_result.asynchronous_similarity.size();
         ++depth_index) {
      Rcpp::List view_matrices = asynchronous_views[depth_index];
      for (std::map<std::string, double>::const_iterator position =
             pair_result.asynchronous_similarity[depth_index].begin();
           position != pair_result.asynchronous_similarity[depth_index].end();
           ++position) {
        const std::map<std::string, int>::const_iterator matrix_position =
          asynchronous_matrix_indices[depth_index].find(position->first);
        if (matrix_position == asynchronous_matrix_indices[depth_index].end()) {
          Rcpp::stop("asynchronous component is outside the declared domain");
        }
        Rcpp::NumericMatrix component_matrix =
          view_matrices[matrix_position->second];
        component_matrix(
          pair_result.first_tip_index,
          pair_result.second_tip_index
        ) = position->second;
        component_matrix(
          pair_result.second_tip_index,
          pair_result.first_tip_index
        ) = position->second;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("sync") = synchronous_matrices,
    Rcpp::Named("async") = asynchronous_views
  );
}

Rcpp::List summarize_matrix_family(
    const Rcpp::List& matrix_family,
    SEXP available_similarity) {
  const int component_count = matrix_family.size();
  Rcpp::NumericVector totals(component_count);
  Rcpp::NumericVector proportions(component_count);
  if (matrix_family.hasAttribute("names")) {
    totals.attr("names") = matrix_family.attr("names");
    proportions.attr("names") = matrix_family.attr("names");
  }
  if (component_count == 0) {
    return Rcpp::List::create(
      Rcpp::Named("totals") = totals,
      Rcpp::Named("tree_wide_proportions") = proportions,
      Rcpp::Named("means") = proportions,
      Rcpp::Named("total_across_matrices") = 0.0
    );
  }

  const Rcpp::NumericMatrix first_matrix = matrix_family[0];
  const bool weighted = available_similarity != R_NilValue;
  Rcpp::NumericMatrix weights;
  double available_total = 0;
  if (weighted) weights = Rcpp::as<Rcpp::NumericMatrix>(available_similarity);
  for (int component_index = 0;
       component_index < component_count;
       ++component_index) {
    const Rcpp::NumericMatrix component_matrix =
      matrix_family[component_index];
    double component_total = 0;
    for (int first_tip_index = 0;
         first_tip_index < component_matrix.nrow() - 1;
         ++first_tip_index) {
      for (int second_tip_index = first_tip_index + 1;
           second_tip_index < component_matrix.ncol();
           ++second_tip_index) {
        const double weight = weighted ?
          weights(first_tip_index, second_tip_index) : 1.0;
        component_total +=
          component_matrix(first_tip_index, second_tip_index) * weight;
        if (weighted && component_index == 0) available_total += weight;
      }
    }
    totals[component_index] = component_total;
  }
  const double total = Rcpp::sum(totals);
  if (total > 0) proportions = totals / total;
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("totals") = totals,
    Rcpp::Named("tree_wide_proportions") = proportions,
    Rcpp::Named("means") = proportions,
    Rcpp::Named("total_across_matrices") = total
  );
  if (weighted) result["available_similarity_total"] = available_total;
  return result;
}

Rcpp::List calculate_matrix_summaries(
    const Rcpp::List& matrix_families,
    const Rcpp::DataFrame& component_metadata,
    bool tree_is_ultrametric,
    const Rcpp::NumericVector& tip_history_times,
    const Rcpp::NumericMatrix& synchronous_available,
    const Rcpp::NumericMatrix& asynchronous_available) {
  SEXP synchronous_weights = R_NilValue;
  SEXP asynchronous_weights = R_NilValue;
  if (!tree_is_ultrametric) {
    synchronous_weights = synchronous_available;
    asynchronous_weights = asynchronous_available;
  }
  const Rcpp::List asynchronous_matrices = matrix_families["async"];
  Rcpp::List asynchronous_summaries(asynchronous_matrices.size());
  if (asynchronous_matrices.hasAttribute("names")) {
    asynchronous_summaries.attr("names") = asynchronous_matrices.attr("names");
  }
  for (int depth_index = 0;
       depth_index < asynchronous_matrices.size();
       ++depth_index) {
    asynchronous_summaries[depth_index] = summarize_matrix_family(
      Rcpp::as<Rcpp::List>(asynchronous_matrices[depth_index]),
      asynchronous_weights
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("sync") = summarize_matrix_family(
      Rcpp::as<Rcpp::List>(matrix_families["sync"]),
      synchronous_weights
    ),
    Rcpp::Named("async") = asynchronous_summaries,
    Rcpp::Named("components") = component_metadata,
    Rcpp::Named("weighting") = Rcpp::List::create(
      Rcpp::Named("tree_is_ultrametric") = tree_is_ultrametric,
      Rcpp::Named("available_similarity_weighting_applied") =
        !tree_is_ultrametric,
      Rcpp::Named("tip_history_times") = tip_history_times,
      Rcpp::Named("synchronous_available_similarity_by_pair") =
        synchronous_available,
      Rcpp::Named("asynchronous_available_similarity_by_pair") =
        asynchronous_available
    )
  );
}

Rcpp::CharacterMatrix prepared_text_matrix(
    const std::vector<TextValue>& values,
    const PreparedScenarioData& prepared_scenario_data,
    SEXP column_names) {
  Rcpp::CharacterMatrix result(
    prepared_scenario_data.tip_count,
    prepared_scenario_data.interval_count
  );
  for (int tip_index = 0;
       tip_index < prepared_scenario_data.tip_count;
       ++tip_index) {
    for (int interval_index = 0;
         interval_index < prepared_scenario_data.interval_count;
         ++interval_index) {
      const TextValue& value = values[scenario_cell_index(
        tip_index,
        interval_index,
        prepared_scenario_data.interval_count
      )];
      if (value.is_missing) {
        result(tip_index, interval_index) = NA_STRING;
      } else {
        result(tip_index, interval_index) = utf8_character(value.value);
      }
    }
  }
  set_matrix_names(result, prepared_scenario_data.tip_label_source, column_names);
  return result;
}

Rcpp::IntegerMatrix prepared_integer_matrix(
    const std::vector<int>& values,
    const PreparedScenarioData& prepared_scenario_data,
    SEXP column_names) {
  Rcpp::IntegerMatrix result(
    prepared_scenario_data.tip_count,
    prepared_scenario_data.interval_count
  );
  for (int tip_index = 0;
       tip_index < prepared_scenario_data.tip_count;
       ++tip_index) {
    for (int interval_index = 0;
         interval_index < prepared_scenario_data.interval_count;
         ++interval_index) {
      result(tip_index, interval_index) = values[scenario_cell_index(
        tip_index,
        interval_index,
        prepared_scenario_data.interval_count
      )];
    }
  }
  set_matrix_names(result, prepared_scenario_data.tip_label_source, column_names);
  return result;
}

Rcpp::CharacterMatrix subset_character_matrix(
    SEXP source_matrix,
    const PreparedScenarioData& prepared_scenario_data,
    SEXP column_names) {
  const Rcpp::CharacterMatrix source(source_matrix);
  Rcpp::CharacterMatrix result(
    prepared_scenario_data.tip_count,
    prepared_scenario_data.interval_count
  );
  for (int column_index = 0;
       column_index < prepared_scenario_data.interval_count;
       ++column_index) {
    for (int tip_index = 0;
         tip_index < prepared_scenario_data.tip_count;
         ++tip_index) {
      result(tip_index, column_index) = source(tip_index, column_index);
    }
  }
  set_matrix_names(result, prepared_scenario_data.tip_label_source, column_names);
  return result;
}

SEXP subset_identifier_matrix(
    SEXP source_matrix,
    const PreparedScenarioData& prepared_scenario_data,
    SEXP column_names) {
  if (TYPEOF(source_matrix) == INTSXP) {
    const Rcpp::IntegerMatrix source(source_matrix);
    Rcpp::IntegerMatrix result(
      prepared_scenario_data.tip_count,
      prepared_scenario_data.interval_count
    );
    for (int column_index = 0;
         column_index < prepared_scenario_data.interval_count;
         ++column_index) {
      for (int tip_index = 0;
           tip_index < prepared_scenario_data.tip_count;
           ++tip_index) {
        result(tip_index, column_index) = source(tip_index, column_index);
      }
    }
    set_matrix_names(result, prepared_scenario_data.tip_label_source, column_names);
    return result;
  }

  const Rcpp::NumericMatrix source(source_matrix);
  Rcpp::NumericMatrix result(
    prepared_scenario_data.tip_count,
    prepared_scenario_data.interval_count
  );
  for (int column_index = 0;
       column_index < prepared_scenario_data.interval_count;
       ++column_index) {
    for (int tip_index = 0;
         tip_index < prepared_scenario_data.tip_count;
         ++tip_index) {
      result(tip_index, column_index) = source(tip_index, column_index);
    }
  }
  set_matrix_names(result, prepared_scenario_data.tip_label_source, column_names);
  return result;
}

Rcpp::List prepared_scenario_data_to_r(
    const PreparedScenarioData& prepared_scenario_data) {
  const int interval_count = prepared_scenario_data.interval_count;
  Rcpp::IntegerVector interval_identifiers(interval_count);
  Rcpp::IntegerVector source_columns(interval_count);
  for (int interval_index = 0; interval_index < interval_count; ++interval_index) {
    interval_identifiers[interval_index] = interval_index + 1;
    source_columns[interval_index] = interval_index + 1;
  }

  Rcpp::List projected_path_matrices(
    prepared_scenario_data.maximum_transition_counts.size()
  );
  Rcpp::CharacterVector projected_path_matrix_names(
    prepared_scenario_data.maximum_transition_counts.size()
  );
  for (std::size_t depth_index = 0;
       depth_index < prepared_scenario_data.maximum_transition_counts.size();
       ++depth_index) {
    projected_path_matrices[depth_index] = prepared_text_matrix(
      prepared_scenario_data.paths_by_requested_depth[depth_index],
      prepared_scenario_data,
      prepared_scenario_data.path_column_names
    );
    projected_path_matrix_names[depth_index] = transition_view_name(
      prepared_scenario_data.maximum_transition_counts[depth_index]
    );
  }
  projected_path_matrices.attr("names") = projected_path_matrix_names;

  Rcpp::LogicalMatrix active_matrix(
    prepared_scenario_data.tip_count,
    prepared_scenario_data.interval_count
  );
  for (int tip_index = 0;
       tip_index < prepared_scenario_data.tip_count;
       ++tip_index) {
    for (int interval_index = 0;
         interval_index < prepared_scenario_data.interval_count;
         ++interval_index) {
      active_matrix(tip_index, interval_index) =
        prepared_scenario_data.lineage_is_active[scenario_cell_index(
          tip_index,
          interval_index,
          prepared_scenario_data.interval_count
        )];
    }
  }
  set_matrix_names(
    active_matrix,
    prepared_scenario_data.tip_label_source,
    R_NilValue
  );

  Rcpp::NumericVector named_tip_end_times = Rcpp::wrap(
    prepared_scenario_data.tip_end_times
  );
  named_tip_end_times.attr("names") = prepared_scenario_data.tip_label_source;
  Rcpp::IntegerVector named_tip_end_columns = Rcpp::wrap(
    prepared_scenario_data.tip_end_columns
  );
  named_tip_end_columns.attr("names") = prepared_scenario_data.tip_label_source;

  return Rcpp::List::create(
    Rcpp::Named("tip_labels") = prepared_scenario_data.tip_label_source,
    Rcpp::Named("time") = Rcpp::DataFrame::create(
      Rcpp::Named("interval_id") = interval_identifiers,
      Rcpp::Named("source_column") = source_columns,
      Rcpp::Named("start") = Rcpp::wrap(prepared_scenario_data.interval_start_times),
      Rcpp::Named("end") = Rcpp::wrap(prepared_scenario_data.interval_end_times),
      Rcpp::Named("duration") = Rcpp::wrap(prepared_scenario_data.interval_durations)
    ),
    Rcpp::Named("states") = subset_character_matrix(
      prepared_scenario_data.state_source,
      prepared_scenario_data,
      prepared_scenario_data.scenario_column_names
    ),
    Rcpp::Named("paths") = subset_character_matrix(
      prepared_scenario_data.path_source,
      prepared_scenario_data,
      prepared_scenario_data.path_column_names
    ),
    Rcpp::Named("paths_using_up_to_transition_count") = projected_path_matrices,
    Rcpp::Named("available_transition_count") = prepared_integer_matrix(
      prepared_scenario_data.available_transition_counts,
      prepared_scenario_data,
      prepared_scenario_data.path_column_names
    ),
    Rcpp::Named("phylo_edge_ids") = subset_identifier_matrix(
      prepared_scenario_data.phylogenetic_edge_identifier_source,
      prepared_scenario_data,
      prepared_scenario_data.phylogenetic_edge_identifier_column_names
    ),
    Rcpp::Named("scenario_edge_ids") = subset_identifier_matrix(
      prepared_scenario_data.scenario_edge_identifier_source,
      prepared_scenario_data,
      prepared_scenario_data.scenario_edge_identifier_column_names
    ),
    Rcpp::Named("scenario_edge_step_ids") = subset_identifier_matrix(
      prepared_scenario_data.scenario_edge_step_identifier_source,
      prepared_scenario_data,
      prepared_scenario_data.scenario_edge_step_identifier_column_names
    ),
    Rcpp::Named("active") = active_matrix,
    Rcpp::Named("duration") = Rcpp::wrap(prepared_scenario_data.interval_durations),
    Rcpp::Named("tip_end_times") = named_tip_end_times,
    Rcpp::Named("tip_end_columns") = named_tip_end_columns,
    Rcpp::Named("root_anchor") =
      Rcpp::as<Rcpp::List>(prepared_scenario_data.root_policy)["anchor"],
    Rcpp::Named("root_policy") = prepared_scenario_data.root_policy,
    Rcpp::Named("maximum_transition_counts") = Rcpp::wrap(
      prepared_scenario_data.maximum_transition_counts
    ),
    Rcpp::Named("state_levels") = prepared_scenario_data.state_level_source,
    Rcpp::Named("time_tolerance") = prepared_scenario_data.time_tolerance
  );
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List scenario_mats_v3_calculate_cpp(
    Rcpp::List trajectory_obj,
    Rcpp::IntegerVector maximum_transition_counts,
    std::string async_metric,
    bool record_complete_output,
    bool record_matrix_output,
    double time_tolerance) {
  const std::vector<int> requested_maximum_transition_counts =
    Rcpp::as<std::vector<int> >(maximum_transition_counts);
  const PreparedScenarioData prepared_scenario_data = prepare_scenario_data(
    trajectory_obj,
    requested_maximum_transition_counts,
    record_complete_output,
    time_tolerance
  );
  const ComponentDefinitions component_definitions =
    build_component_definitions(prepared_scenario_data);
  const Rcpp::DataFrame component_metadata =
    component_definitions_to_r(component_definitions);

  if (!record_complete_output) {
    const CompactScenarioData compact =
      build_compact_scenario_data(prepared_scenario_data);
    const CompactMatrixStorage storage = calculate_compact_matrix_storage(
      prepared_scenario_data,
      component_definitions,
      compact,
      async_metric,
      record_matrix_output
    );
    const Rcpp::NumericVector tip_history_times =
      tip_history_times_to_r(prepared_scenario_data, compact);
    const bool tree_is_ultrametric =
      Rcpp::max(tip_history_times) - Rcpp::min(tip_history_times) <=
        prepared_scenario_data.time_tolerance;
    const Rcpp::List summaries = compact_summaries_to_r(
      storage,
      prepared_scenario_data,
      component_definitions,
      component_metadata,
      tree_is_ultrametric,
      tip_history_times,
      record_matrix_output
    );
    if (!record_matrix_output) {
      return Rcpp::List::create(
        Rcpp::Named("summaries") = summaries
      );
    }
    const Rcpp::List matrices = compact_matrix_storage_to_r(
      storage,
      prepared_scenario_data,
      component_definitions
    );
    return Rcpp::List::create(
      Rcpp::Named("matrices") = matrices,
      Rcpp::Named("summaries") = summaries
    );
  }

  std::vector<PairCalculationResult> pair_results;
  pair_results.reserve(
    static_cast<std::size_t>(prepared_scenario_data.tip_count) *
      static_cast<std::size_t>(prepared_scenario_data.tip_count - 1) / 2
  );
  for (int first_tip_index = 0;
       first_tip_index < prepared_scenario_data.tip_count - 1;
       ++first_tip_index) {
    for (int second_tip_index = first_tip_index + 1;
         second_tip_index < prepared_scenario_data.tip_count;
         ++second_tip_index) {
      pair_results.push_back(calculate_one_pair_v21(
        prepared_scenario_data,
        component_definitions,
        first_tip_index,
        second_tip_index,
        async_metric,
        record_complete_output
      ));
    }
  }
  const Rcpp::List matrix_families = matrix_families_to_r(
    pair_results,
    prepared_scenario_data,
    component_definitions
  );
  const Rcpp::NumericMatrix synchronous_available =
    total_comparable_time_matrix_to_r(
      pair_results,
      prepared_scenario_data
    );
  const Rcpp::NumericMatrix asynchronous_available =
    asynchronous_available_similarity_to_r(
      pair_results,
      prepared_scenario_data,
      async_metric
    );
  const Rcpp::NumericVector tip_history_times =
    tip_history_times_to_r(prepared_scenario_data);
  const bool tree_is_ultrametric =
    Rcpp::max(tip_history_times) - Rcpp::min(tip_history_times) <=
      prepared_scenario_data.time_tolerance;
  const Rcpp::List matrix_summaries = calculate_matrix_summaries(
    matrix_families,
    component_metadata,
    tree_is_ultrametric,
    tip_history_times,
    synchronous_available,
    asynchronous_available
  );
  if (!record_complete_output) {
    return Rcpp::List::create(
      Rcpp::Named("matrices") = matrix_families,
      Rcpp::Named("summaries") = matrix_summaries
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("intervals") = prepared_scenario_data_to_r(prepared_scenario_data),
    Rcpp::Named("matrices") = matrix_families,
    Rcpp::Named("total_comparable_time_by_pair") = synchronous_available,
    Rcpp::Named("available_similarity_by_pair") = Rcpp::List::create(
      Rcpp::Named("synchronous") = synchronous_available,
      Rcpp::Named("asynchronous") = asynchronous_available
    ),
    Rcpp::Named("pair_details") = Rcpp::List::create(
      Rcpp::Named("interval_classification") = interval_classifications_v21_to_r(
        pair_results,
        prepared_scenario_data
      ),
      Rcpp::Named("sync_parallel_time_by_path") = synchronous_parallel_times_to_r(
        pair_results,
        prepared_scenario_data
      ),
      Rcpp::Named("async_independent_time_by_path") =
        asynchronous_independent_times_to_r(
          pair_results,
          prepared_scenario_data
        ),
      Rcpp::Named("async_parallel_similarity_by_path") =
        asynchronous_parallel_similarities_to_r(
          pair_results,
          prepared_scenario_data
        )
    ),
    Rcpp::Named("summaries") = matrix_summaries,
    Rcpp::Named("components") = component_metadata,
    Rcpp::Named("metadata") = Rcpp::List::create(
      Rcpp::Named("maximum_transition_counts") = maximum_transition_counts,
      Rcpp::Named("async_metric") = async_metric,
      Rcpp::Named("tree_is_ultrametric") = tree_is_ultrametric,
      Rcpp::Named("available_similarity_weighting_applied") =
        !tree_is_ultrametric,
      Rcpp::Named("normalization") = Rcpp::List::create(
        Rcpp::Named("synchronous") = "shared_active_time",
        Rcpp::Named("asynchronous_bhattacharyya") =
          "sqrt(first_history_time * second_history_time)",
        Rcpp::Named("asynchronous_minimum") =
          "min(first_history_time, second_history_time)"
      ),
      Rcpp::Named("time_tolerance") = prepared_scenario_data.time_tolerance
    )
  );
}
