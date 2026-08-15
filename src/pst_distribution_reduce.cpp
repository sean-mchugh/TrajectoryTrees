#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <utility>
#include <vector>

using namespace Rcpp;

static bool pst_distribution_integer_vectors_equal(
    const Rcpp::IntegerVector& left,
    const Rcpp::IntegerVector& right) {
  if (left.size() != right.size()) {
    return false;
  }
  for (R_xlen_t index = 0; index < left.size(); ++index) {
    if (left[index] != right[index]) {
      return false;
    }
  }
  return true;
}

static double pst_distribution_type7(
    std::vector<double>& values,
    const double probability) {
  if (values.empty()) {
    return NA_REAL;
  }
  const double h = (static_cast<double>(values.size()) - 1.0) * probability;
  const std::size_t lower_rank = static_cast<std::size_t>(std::floor(h));
  const std::size_t upper_rank = static_cast<std::size_t>(std::ceil(h));
  std::nth_element(values.begin(), values.begin() + lower_rank, values.end());
  const double lower = values[lower_rank];
  if (lower_rank == upper_rank) {
    return lower;
  }
  std::nth_element(values.begin(), values.begin() + upper_rank, values.end());
  const double upper = values[upper_rank];
  return lower + (h - static_cast<double>(lower_rank)) * (upper - lower);
}

static double pst_distribution_implicit_zero_order_stat(
    std::vector<double>& present,
    const std::size_t rank,
    const std::size_t implicit_zeros,
    const std::size_t negative_count,
    const std::size_t explicit_zero_count) {
  const std::size_t total_zero_count = implicit_zeros + explicit_zero_count;
  if (rank < negative_count) {
    std::nth_element(present.begin(), present.begin() + rank, present.end());
    return present[rank];
  }
  if (rank < negative_count + total_zero_count) {
    return 0.0;
  }
  const std::size_t present_rank = rank - implicit_zeros;
  std::nth_element(present.begin(), present.begin() + present_rank, present.end());
  return present[present_rank];
}

static double pst_distribution_type7_implicit_zeros(
    std::vector<double>& present,
    const std::size_t implicit_zeros,
    const std::size_t negative_count,
    const std::size_t explicit_zero_count,
    const double probability) {
  const std::size_t count = present.size() + implicit_zeros;
  if (count == 0) {
    return NA_REAL;
  }
  const double h = (static_cast<double>(count) - 1.0) * probability;
  const std::size_t lower_rank = static_cast<std::size_t>(std::floor(h));
  const std::size_t upper_rank = static_cast<std::size_t>(std::ceil(h));
  const double lower = pst_distribution_implicit_zero_order_stat(
    present, lower_rank, implicit_zeros, negative_count, explicit_zero_count
  );
  if (lower_rank == upper_rank) {
    return lower;
  }
  const double upper = pst_distribution_implicit_zero_order_stat(
    present, upper_rank, implicit_zeros, negative_count, explicit_zero_count
  );
  return lower + (h - static_cast<double>(lower_rank)) * (upper - lower);
}

static double pst_distribution_weighted_order_stat(
    const std::vector<std::pair<double, int> >& sorted,
    const std::size_t rank) {
  std::size_t cumulative = 0;
  for (const auto& item : sorted) {
    cumulative += static_cast<std::size_t>(item.second);
    if (rank < cumulative) {
      return item.first;
    }
  }
  return sorted.back().first;
}

static double pst_distribution_type7_weighted(
    const std::vector<std::pair<double, int> >& sorted,
    const std::size_t count,
    const double probability) {
  if (count == 0) {
    return NA_REAL;
  }
  const double h = (static_cast<double>(count) - 1.0) * probability;
  const std::size_t lower_rank = static_cast<std::size_t>(std::floor(h));
  const std::size_t upper_rank = static_cast<std::size_t>(std::ceil(h));
  const double lower = pst_distribution_weighted_order_stat(sorted, lower_rank);
  if (lower_rank == upper_rank) {
    return lower;
  }
  const double upper = pst_distribution_weighted_order_stat(sorted, upper_rank);
  return lower + (h - static_cast<double>(lower_rank)) * (upper - lower);
}

// [[Rcpp::export]]
Rcpp::List pst_distribution_cpp_reduce_exact(
    const Rcpp::NumericMatrix& values,
    const Rcpp::NumericVector& probs) {
  if (probs.size() != 2 || !R_finite(probs[0]) || !R_finite(probs[1]) ||
      probs[0] < 0.0 || probs[1] > 1.0 || probs[0] > probs[1]) {
    Rcpp::stop("probs must contain ordered lower and upper probabilities");
  }
  const int sample_count = values.nrow();
  const int element_count = values.ncol();
  if (sample_count < 1) {
    Rcpp::stop("C++ distribution reduction requires at least one sample");
  }
  NumericVector lower(element_count, NA_REAL);
  NumericVector mean(element_count, NA_REAL);
  NumericVector median(element_count, NA_REAL);
  NumericVector upper(element_count, NA_REAL);
  IntegerVector observed(element_count, 0);
  std::vector<double> finite;
  finite.reserve(sample_count);
  // Reduce each output element independently while reusing one scratch vector.
  for (int element = 0; element < element_count; ++element) {
    finite.clear();
    double total = 0.0;
    // Scan posterior samples once for mean, count, infinity rejection, and selection input.
    for (int sample = 0; sample < sample_count; ++sample) {
      const double value = values(sample, element);
      if (R_IsNA(value) || R_IsNaN(value)) {
        continue;
      }
      if (!R_finite(value)) {
        Rcpp::stop("C++ distribution block contains an infinite value at element %d", element + 1);
      }
      finite.push_back(value);
      total += value;
    }
    observed[element] = static_cast<int>(finite.size());
    if (finite.empty()) {
      continue;
    }
    mean[element] = total / static_cast<double>(finite.size());
    lower[element] = pst_distribution_type7(finite, probs[0]);
    median[element] = pst_distribution_type7(finite, 0.5);
    upper[element] = pst_distribution_type7(finite, probs[1]);
  }
  return List::create(
    Named("lower") = lower,
    Named("mean") = mean,
    Named("median") = median,
    Named("upper") = upper,
    Named("n_observed") = observed
  );
}

// [[Rcpp::export]]
Rcpp::List pst_distribution_cpp_align_reduce_tt(
    const Rcpp::List& matrices,
    const Rcpp::List& row_maps,
    const Rcpp::List& column_maps,
    const Rcpp::LogicalVector& structural_time,
    const Rcpp::NumericVector& canonical_time,
    const Rcpp::NumericVector& probs,
    const int elements_per_block) {
  const int sample_count = matrices.size();
  const int canonical_rows = canonical_time.size();
  const int canonical_columns = structural_time.size();
  if (sample_count < 1 || row_maps.size() != sample_count ||
      column_maps.size() != sample_count) {
    Rcpp::stop("TT C++ alignment maps are not sample-aligned");
  }
  if (elements_per_block < 1) {
    Rcpp::stop("TT C++ block must contain at least one element");
  }
  const int element_count = canonical_rows * canonical_columns;
  NumericVector lower(element_count, NA_REAL);
  NumericVector mean(element_count, NA_REAL);
  NumericVector median(element_count, NA_REAL);
  NumericVector upper(element_count, NA_REAL);
  IntegerVector observed(element_count, 0);
  std::vector<NumericMatrix> sample_matrices;
  std::vector<IntegerVector> sample_row_maps;
  std::vector<IntegerVector> sample_column_maps;
  std::vector<int> sample_weights;
  sample_matrices.reserve(sample_count);
  sample_row_maps.reserve(sample_count);
  sample_column_maps.reserve(sample_count);
  std::map<SEXP, std::vector<int> > representatives_by_matrix;
  // Collapse repeated matrix objects to weighted representatives before hot loops.
  for (int sample = 0; sample < sample_count; ++sample) {
    SEXP matrix_key = matrices[sample];
    const IntegerVector current_rows = as<IntegerVector>(row_maps[sample]);
    const IntegerVector current_columns = as<IntegerVector>(column_maps[sample]);
    int matching_representative = -1;
    // Group one shared matrix only when its row and column interpretations also match.
    for (const int representative : representatives_by_matrix[matrix_key]) {
      if (pst_distribution_integer_vectors_equal(
            current_rows, sample_row_maps[representative]) &&
          pst_distribution_integer_vectors_equal(
            current_columns, sample_column_maps[representative])) {
        matching_representative = representative;
        break;
      }
    }
    if (matching_representative >= 0) {
      sample_weights[matching_representative] += 1;
      continue;
    }
    const int representative = static_cast<int>(sample_matrices.size());
    representatives_by_matrix[matrix_key].push_back(representative);
    sample_matrices.push_back(as<NumericMatrix>(matrices[sample]));
    sample_row_maps.push_back(current_rows);
    sample_column_maps.push_back(current_columns);
    sample_weights.push_back(1);
    if (sample_row_maps.back().size() != canonical_rows ||
        sample_column_maps.back().size() != canonical_columns) {
      Rcpp::stop("TT C++ alignment index dimensions do not match");
    }
  }
  const int representative_count = static_cast<int>(sample_matrices.size());
  const bool weighted_mode = representative_count < sample_count;
  std::vector<std::vector<int> > present_samples(canonical_columns);
  std::vector<int> present_weights(canonical_columns, 0);
  // Index only samples that physically contain each union column; all others are implicit zeros.
  for (int canonical_column = 0; canonical_column < canonical_columns;
       ++canonical_column) {
    if (structural_time[canonical_column] == TRUE) {
      continue;
    }
    for (int sample = 0; sample < representative_count; ++sample) {
      if (sample_column_maps[sample][canonical_column] >= 0) {
        present_samples[canonical_column].push_back(sample);
        present_weights[canonical_column] += sample_weights[sample];
      }
    }
  }
  std::vector<double> finite;
  finite.reserve(sample_count);
  std::vector<std::pair<double, int> > weighted_values;
  weighted_values.reserve(representative_count + 1);
  // Process canonical TT cells in bounded blocks without materializing aligned matrices.
  for (int block_start = 0; block_start < element_count;
       block_start += elements_per_block) {
    const int block_end = std::min(element_count, block_start + elements_per_block);
    // Reduce every column-major canonical TT cell directly from source matrices.
    for (int element = block_start; element < block_end; ++element) {
      const int canonical_row = element % canonical_rows;
      const int canonical_column = element / canonical_rows;
      if (structural_time[canonical_column] == TRUE) {
        lower[element] = canonical_time[canonical_row];
        mean[element] = canonical_time[canonical_row];
        median[element] = canonical_time[canonical_row];
        upper[element] = canonical_time[canonical_row];
        observed[element] = sample_count;
        continue;
      }
      finite.clear();
      double total = 0.0;
      std::size_t negative_count = 0;
      std::size_t explicit_zero_count = 0;
      const std::size_t implicit_zeros = static_cast<std::size_t>(
        sample_count - present_weights[canonical_column]
      );
      weighted_values.clear();
      std::size_t finite_weight = 0;
      // Resolve only samples containing this union column; absent samples are exact zeros.
      for (const int sample : present_samples[canonical_column]) {
        const int source_column = sample_column_maps[sample][canonical_column];
        const int source_row = sample_row_maps[sample][canonical_row];
        if (source_row < 0 || source_row >= sample_matrices[sample].nrow() ||
            source_column >= sample_matrices[sample].ncol()) {
          Rcpp::stop("TT C++ alignment index is outside its source matrix");
        }
        const double value = sample_matrices[sample](source_row, source_column);
        if (R_IsNA(value) || R_IsNaN(value)) {
          continue;
        }
        if (!R_finite(value)) {
          Rcpp::stop("TT C++ alignment contains an infinite value at element %d", element + 1);
        }
        finite.push_back(value);
        total += value * static_cast<double>(sample_weights[sample]);
        finite_weight += static_cast<std::size_t>(sample_weights[sample]);
        if (weighted_mode) {
          weighted_values.push_back(std::make_pair(value, sample_weights[sample]));
        }
        negative_count += value < 0.0;
        explicit_zero_count += value == 0.0;
      }
      const std::size_t observed_count = finite_weight + implicit_zeros;
      observed[element] = static_cast<int>(observed_count);
      if (observed_count == 0) {
        continue;
      }
      mean[element] = total / static_cast<double>(observed_count);
      if (weighted_mode) {
        if (implicit_zeros > 0) {
          weighted_values.push_back(std::make_pair(0.0, static_cast<int>(implicit_zeros)));
        }
        std::sort(weighted_values.begin(), weighted_values.end(),
                  [](const auto& left, const auto& right) {
                    return left.first < right.first;
                  });
        lower[element] = pst_distribution_type7_weighted(
          weighted_values, observed_count, probs[0]
        );
        median[element] = pst_distribution_type7_weighted(
          weighted_values, observed_count, 0.5
        );
        upper[element] = pst_distribution_type7_weighted(
          weighted_values, observed_count, probs[1]
        );
      } else {
        lower[element] = pst_distribution_type7_implicit_zeros(
          finite, implicit_zeros, negative_count, explicit_zero_count, probs[0]
        );
        median[element] = pst_distribution_type7_implicit_zeros(
          finite, implicit_zeros, negative_count, explicit_zero_count, 0.5
        );
        upper[element] = pst_distribution_type7_implicit_zeros(
          finite, implicit_zeros, negative_count, explicit_zero_count, probs[1]
        );
      }
    }
    Rcpp::checkUserInterrupt();
  }
  return List::create(
    Named("lower") = lower,
    Named("mean") = mean,
    Named("median") = median,
    Named("upper") = upper,
    Named("n_observed") = observed
  );
}
