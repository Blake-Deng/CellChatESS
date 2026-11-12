
#include <cstdlib>
#include <fstream>
#include <cstdio>
#include <utility>
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

static inline double quantile_type7_sorted(const std::vector<double>& x, double p) {
  const int n = x.size();
  if (n == 0) return NA_REAL;
  if (n == 1) return x[0];
  const double h = 1.0 + (n - 1.0) * p;
  const int hf = static_cast<int>(std::floor(h));
  const double frac = h - hf;
  const int i0 = hf - 1;
  if (i0 < 0) return x.front();
  if (i0 >= n - 1) return x.back();
  return x[i0] + frac * (x[i0 + 1] - x[i0]);
}

static inline double tri_mean_sorted_copy(const std::vector<double>& x) {
  if (x.empty()) return NA_REAL;
  std::vector<double> sorted(x);
  std::sort(sorted.begin(), sorted.end());
  const double q1 = quantile_type7_sorted(sorted, 0.25);
  const double q2 = quantile_type7_sorted(sorted, 0.50);
  const double q3 = quantile_type7_sorted(sorted, 0.75);
  return (q1 + 2.0 * q2 + q3) / 4.0;
}

static inline double tri_mean_multiselect(const std::vector<double>& x) {
  const int n = x.size();
  if (n == 0) return NA_REAL;
  if (n == 1) return x[0];

  const double probs[3] = {0.25, 0.50, 0.75};
  int lower[3];
  double frac[3];
  std::vector<int> needed;
  needed.reserve(6);
  for (int q = 0; q < 3; ++q) {
    const double h = 1.0 + (n - 1.0) * probs[q];
    const int hf = static_cast<int>(std::floor(h));
    lower[q] = hf - 1;
    frac[q] = h - hf;
    needed.push_back(lower[q]);
    needed.push_back(lower[q] + 1);
  }
  std::sort(needed.begin(), needed.end());
  needed.erase(std::unique(needed.begin(), needed.end()), needed.end());

  std::vector<double> buffer(x);
  std::vector<double> selected(needed.size());
  for (int i = 0; i < static_cast<int>(needed.size()); ++i) {
    const int k = needed[i];
    std::nth_element(buffer.begin(), buffer.begin() + k, buffer.end());
    selected[i] = buffer[k];
  }
  const auto selected_value = [&](int k) {
    const auto it = std::lower_bound(needed.begin(), needed.end(), k);
    return selected[static_cast<int>(it - needed.begin())];
  };

  double quantiles[3];
  for (int q = 0; q < 3; ++q) {
    const double lo = selected_value(lower[q]);
    const double hi = selected_value(lower[q] + 1);
    quantiles[q] = lo + frac[q] * (hi - lo);
  }
  return (quantiles[0] + 2.0 * quantiles[1] + quantiles[2]) / 4.0;
}

static inline double tri_mean_unsorted(const std::vector<double>& x) {
  return x.size() < 4096 ? tri_mean_sorted_copy(x) : tri_mean_multiselect(x);
}

// [[Rcpp::export]]
NumericMatrix group_tri_mean_cpp(NumericMatrix data, IntegerVector group, int K) {
  const int G = data.nrow();
  const int C = data.ncol();
  NumericMatrix out(G, K);

  std::vector< std::vector<int> > cells(K);
  for (int c = 0; c < C; ++c) {
    int g = group[c] - 1;
    if (g >= 0 && g < K) cells[g].push_back(c);
  }

  std::vector<double> vals;
  for (int gene = 0; gene < G; ++gene) {
    for (int k = 0; k < K; ++k) {
      vals.clear();
      vals.reserve(cells[k].size());
      for (int idx : cells[k]) {
        double v = data(gene, idx);
        if (!NumericVector::is_na(v)) vals.push_back(v);
      }
      out(gene, k) = tri_mean_unsorted(vals);
    }
  }
  return out;
}

// [[Rcpp::export]]
List group_tri_mean_boot_cpp(NumericMatrix data, IntegerMatrix groupBoot, int K) {
  const int G = data.nrow();
  const int C = data.ncol();
  const int B = groupBoot.ncol();
  NumericVector out(Dimension(G, K, B));

  std::vector< std::vector<int> > cells(K);
  std::vector<double> vals;
  for (int b = 0; b < B; ++b) {
    for (int k = 0; k < K; ++k) cells[k].clear();
    for (int c = 0; c < C; ++c) {
      int g = groupBoot(c, b) - 1;
      if (g >= 0 && g < K) cells[g].push_back(c);
    }
    const int offset = G * K * b;
    for (int gene = 0; gene < G; ++gene) {
      for (int k = 0; k < K; ++k) {
        vals.clear();
        vals.reserve(cells[k].size());
        for (int idx : cells[k]) {
          double v = data(gene, idx);
          if (!NumericVector::is_na(v)) vals.push_back(v);
        }
        out[gene + G * k + offset] = tri_mean_unsorted(vals);
      }
    }
  }
  return List::create(Named("avg_boot") = out);
}

// Resolve the gene/complex distinction once, before any OpenMP region.
static inline std::vector<unsigned char> direct_gene_flags(const IntegerMatrix& index) {
  std::vector<unsigned char> flags(index.nrow(), 1);
  if (index.hasAttribute("direct_gene")) {
    LogicalVector direct = index.attr("direct_gene");
    if (direct.size() != index.nrow()) stop("Invalid direct-gene metadata");
    for (int i = 0; i < index.nrow(); ++i) flags[i] = direct[i] == TRUE;
  }
  return flags;
}

static inline double avg_expr_matrix(const NumericMatrix& avg, const IntegerMatrix& idx, int row, int k, bool direct_gene = true) {
  double sum_log = 0.0;
  double single_value = NA_REAL;
  int n = 0;
  for (int j = 0; j < idx.ncol(); ++j) {
    const int g = idx(row, j) - 1;
    if (g < 0) continue;
    const double x = avg(g, k);
    if (NumericVector::is_na(x)) continue;
    // A single gene uses the group summary directly, as in CellChat.
    // For multiple components, retain the original log-sum order.
    if (n == 0) {
      single_value = x;
    } else {
      if (n == 1) sum_log += std::log(single_value);
      sum_log += std::log(x);
    }
    ++n;
  }
  if (n == 0) return NA_REAL;
  if (n == 1) return direct_gene ? single_value : std::exp(std::log(single_value));
  return std::exp(sum_log / static_cast<double>(n));
}

static inline double avg_expr_boot(const NumericVector& avgBoot, int G, int K,
                                   const IntegerMatrix& idx, int row, int k, int b, bool direct_gene = true) {
  double sum_log = 0.0;
  double single_value = NA_REAL;
  int n = 0;
  const int boot_offset = G * K * b;
  for (int j = 0; j < idx.ncol(); ++j) {
    const int g = idx(row, j) - 1;
    if (g < 0) continue;
    const double x = avgBoot[g + G * k + boot_offset];
    if (NumericVector::is_na(x)) continue;
    // A single gene uses the group summary directly, as in CellChat.
    // For multiple components, retain the original log-sum order.
    if (n == 0) {
      single_value = x;
    } else {
      if (n == 1) sum_log += std::log(single_value);
      sum_log += std::log(x);
    }
    ++n;
  }
  if (n == 0) return NA_REAL;
  if (n == 1) return direct_gene ? single_value : std::exp(std::log(single_value));
  return std::exp(sum_log / static_cast<double>(n));
}

static inline double coreceptor_factor_matrix(const NumericMatrix& avg, const IntegerMatrix& idx, int row, int k) {
  double out = 1.0;
  for (int j = 0; j < idx.ncol(); ++j) {
    const int g = idx(row, j) - 1;
    if (g < 0) continue;
    const double x = avg(g, k);
    if (NumericVector::is_na(x)) continue;
    out *= 1.0 + x;
  }
  return out;
}

static inline double coreceptor_factor_boot(const NumericVector& avgBoot, int G, int K,
                                            const IntegerMatrix& idx, int row, int k, int b) {
  double out = 1.0;
  const int boot_offset = G * K * b;
  for (int j = 0; j < idx.ncol(); ++j) {
    const int g = idx(row, j) - 1;
    if (g < 0) continue;
    const double x = avgBoot[g + G * k + boot_offset];
    if (NumericVector::is_na(x)) continue;
    out *= 1.0 + x;
  }
  return out;
}

static inline double agonist_factor_value(double x, double Kh_n, double n_power) {
  const double xn = std::pow(x, n_power);
  return 1.0 + xn / (Kh_n + xn);
}

static inline double antagonist_factor_value(double x, double Kh_n, double n_power) {
  const double xn = std::pow(x, n_power);
  return Kh_n / (Kh_n + xn);
}

static inline double cofactor_factor_matrix(const NumericMatrix& avg, const IntegerMatrix& idx, int row, int k,
                                            bool agonist_mode, double Kh_n, double n_power) {
  double out = 1.0;
  for (int j = 0; j < idx.ncol(); ++j) {
    const int g = idx(row, j) - 1;
    if (g < 0) continue;
    const double x = avg(g, k);
    if (NumericVector::is_na(x)) continue;
    out *= agonist_mode ? agonist_factor_value(x, Kh_n, n_power) : antagonist_factor_value(x, Kh_n, n_power);
  }
  return out;
}

static inline double cofactor_factor_boot(const NumericVector& avgBoot, int G, int K,
                                          const IntegerMatrix& idx, int row, int k, int b,
                                          bool agonist_mode, double Kh_n, double n_power) {
  double out = 1.0;
  const int boot_offset = G * K * b;
  for (int j = 0; j < idx.ncol(); ++j) {
    const int g = idx(row, j) - 1;
    if (g < 0) continue;
    const double x = avgBoot[g + G * k + boot_offset];
    if (NumericVector::is_na(x)) continue;
    out *= agonist_mode ? agonist_factor_value(x, Kh_n, n_power) : antagonist_factor_value(x, Kh_n, n_power);
  }
  return out;
}

static inline bool index_row_empty(const IntegerMatrix& idx, int row) {
  for (int j = 0; j < idx.ncol(); ++j) {
    if (idx(row, j) > 0) return false;
  }
  return true;
}

// [[Rcpp::export]]
List cellchat_prob_from_avg_cpp(
    NumericMatrix avg,
    NumericVector avgBoot,
    IntegerMatrix ligandIdx,
    IntegerMatrix receptorIdx,
    IntegerMatrix coAIdx,
    IntegerMatrix coIIdx,
    IntegerMatrix agonistIdx,
    IntegerMatrix antagonistIdx,
    LogicalVector hasAgonist,
    LogicalVector hasAntagonist,
    double Kh,
    double n_power) {

  const IntegerVector dims = avgBoot.attr("dim");
  const int G = dims[0];
  const int K = dims[1];
  const int nboot = dims[2];
  const int nLR = ligandIdx.nrow();
  const auto ligand_direct = direct_gene_flags(ligandIdx);
  const auto receptor_direct = direct_gene_flags(receptorIdx);
  const int KK = K * K;
  const R_xlen_t KK_x = static_cast<R_xlen_t>(KK);
  const double Kh_n = std::pow(Kh, n_power);

  NumericVector prob(Dimension(K, K, nLR));
  NumericVector pval(Dimension(K, K, nLR));

  std::vector<double> L(K), R(K), ago(K), ant(K);
  std::vector<double> Lb(K), Rb(K), agoB(K), antB(K);

  for (int lr = 0; lr < nLR; ++lr) {
    for (int k = 0; k < K; ++k) {
      L[k] = avg_expr_matrix(avg, ligandIdx, lr, k, ligand_direct[lr]);
      R[k] = avg_expr_matrix(avg, receptorIdx, lr, k, receptor_direct[lr]);
      const double coA = coreceptor_factor_matrix(avg, coAIdx, lr, k);
      const double coI = coreceptor_factor_matrix(avg, coIIdx, lr, k);
      R[k] = R[k] * coA / coI;
      ago[k] = hasAgonist[lr] ? cofactor_factor_matrix(avg, agonistIdx, lr, k, true, Kh_n, n_power) : 1.0;
      ant[k] = hasAntagonist[lr] ? cofactor_factor_matrix(avg, antagonistIdx, lr, k, false, Kh_n, n_power) : 1.0;
    }

    std::vector<double> pnull(KK);
    bool all_zero = true;
    for (int s = 0; s < K; ++s) {
      for (int r = 0; r < K; ++r) {
        const double dataLR = L[s] * R[r];
        double p1 = 0.0;
        if (dataLR > 0.0) {
          const double dataLR_n = std::pow(dataLR, n_power);
          p1 = dataLR_n / (Kh_n + dataLR_n);
        }
        double val = p1 * (ago[s] * ago[r]) * (ant[s] * ant[r]);
        const int idx = s + K * r;
        pnull[idx] = val;
        prob[static_cast<R_xlen_t>(idx) + KK_x * lr] = val;
        if (val != 0.0) all_zero = false;
      }
    }

    if (all_zero) {
      for (int idx = 0; idx < KK; ++idx) pval[static_cast<R_xlen_t>(idx) + KK_x * lr] = 1.0;
      continue;
    }

    std::vector<int> reject(KK, 0);
    for (int b = 0; b < nboot; ++b) {
      for (int k = 0; k < K; ++k) {
        Lb[k] = avg_expr_boot(avgBoot, G, K, ligandIdx, lr, k, b, ligand_direct[lr]);
        Rb[k] = avg_expr_boot(avgBoot, G, K, receptorIdx, lr, k, b, receptor_direct[lr]);
        const double coA = coreceptor_factor_boot(avgBoot, G, K, coAIdx, lr, k, b);
        const double coI = coreceptor_factor_boot(avgBoot, G, K, coIIdx, lr, k, b);
        Rb[k] = Rb[k] * coA / coI;
        agoB[k] = hasAgonist[lr] ? cofactor_factor_boot(avgBoot, G, K, agonistIdx, lr, k, b, true, Kh_n, n_power) : 1.0;
        antB[k] = hasAntagonist[lr] ? cofactor_factor_boot(avgBoot, G, K, antagonistIdx, lr, k, b, false, Kh_n, n_power) : 1.0;
      }
      for (int s = 0; s < K; ++s) {
        for (int r = 0; r < K; ++r) {
          const double dataLR = Lb[s] * Rb[r];
          double p1 = 0.0;
          if (dataLR > 0.0) {
            const double dataLR_n = std::pow(dataLR, n_power);
            p1 = dataLR_n / (Kh_n + dataLR_n);
          }
          const double val = p1 * (agoB[s] * agoB[r]) * (antB[s] * antB[r]);
          const int idx = s + K * r;
          if (val - pnull[idx] > 0.0) reject[idx] += 1;
        }
      }
    }

    for (int idx = 0; idx < KK; ++idx) {
      double pv = static_cast<double>(reject[idx]) / static_cast<double>(nboot);
      if (pnull[idx] == 0.0) pv = 1.0;
      pval[static_cast<R_xlen_t>(idx) + KK_x * lr] = pv;
    }
  }

  return List::create(Named("prob") = prob, Named("pval") = pval);
}

static inline double tri_mean_gene_cells(const NumericMatrix& data, int gene0,
                                         const std::vector<int>& cell_idx) {
  std::vector<double> vals;
  vals.reserve(cell_idx.size());
  for (int c : cell_idx) {
    const double v = data(gene0, c);
    if (!NumericVector::is_na(v)) vals.push_back(v);
  }
  return tri_mean_unsorted(vals);
}

static inline void add_index_genes(std::vector<int>& genes, const IntegerMatrix& idx, int row) {
  for (int j = 0; j < idx.ncol(); ++j) {
    const int g0 = idx(row, j) - 1;
    if (g0 >= 0) genes.push_back(g0);
  }
}

static inline double avg_expr_local(const std::vector<double>& bootAvg,
                                    const std::vector<int>& gene_to_local,
                                    const IntegerMatrix& idx,
                                    int row, int k, int K, bool direct_gene = true) {
  double sum_log = 0.0;
  double single_value = NA_REAL;
  int n = 0;
  for (int j = 0; j < idx.ncol(); ++j) {
    const int g0 = idx(row, j) - 1;
    if (g0 < 0) continue;
    const int local = gene_to_local[g0];
    if (local < 0) continue;
    const double val = bootAvg[local * K + k];
    if (NumericVector::is_na(val)) continue;
    // A single gene uses the group summary directly, as in CellChat.
    // For multiple components, retain the original log-sum order.
    if (n == 0) {
      single_value = val;
    } else {
      if (n == 1) sum_log += std::log(single_value);
      sum_log += std::log(val);
    }
    ++n;
  }
  if (n == 0) return NA_REAL;
  if (n == 1) return direct_gene ? single_value : std::exp(std::log(single_value));
  return std::exp(sum_log / static_cast<double>(n));
}

static inline double coreceptor_factor_local(const std::vector<double>& bootAvg,
                                             const std::vector<int>& gene_to_local,
                                             const IntegerMatrix& idx,
                                             int row, int k, int K) {
  double out = 1.0;
  for (int j = 0; j < idx.ncol(); ++j) {
    const int g0 = idx(row, j) - 1;
    if (g0 < 0) continue;
    const int local = gene_to_local[g0];
    if (local < 0) continue;
    const double x = bootAvg[local * K + k];
    if (NumericVector::is_na(x)) continue;
    out *= 1.0 + x;
  }
  return out;
}

static inline double cofactor_factor_local(const std::vector<double>& bootAvg,
                                           const std::vector<int>& gene_to_local,
                                           const IntegerMatrix& idx,
                                           int row, int k, int K,
                                           bool agonist_mode, double Kh_n, double n_power) {
  double out = 1.0;
  for (int j = 0; j < idx.ncol(); ++j) {
    const int g0 = idx(row, j) - 1;
    if (g0 < 0) continue;
    const int local = gene_to_local[g0];
    if (local < 0) continue;
    const double x = bootAvg[local * K + k];
    if (NumericVector::is_na(x)) continue;
    out *= agonist_mode ? agonist_factor_value(x, Kh_n, n_power) : antagonist_factor_value(x, Kh_n, n_power);
  }
  return out;
}

// [[Rcpp::export]]
List cellchat_prob_ess_cpp(
    NumericMatrix data,
    IntegerVector group,
    IntegerMatrix groupBoot,
    NumericMatrix avg,
    IntegerMatrix ligandIdx,
    IntegerMatrix receptorIdx,
    IntegerMatrix coAIdx,
    IntegerMatrix coIIdx,
    IntegerMatrix agonistIdx,
    IntegerMatrix antagonistIdx,
    LogicalVector hasAgonist,
    LogicalVector hasAntagonist,
    double Kh,
    double n_power,
    int threads) {

  const int G = data.nrow();
  const int C = data.ncol();
  const int K = Rcpp::max(group);
  const int nboot = groupBoot.ncol();
  const int nLR = ligandIdx.nrow();
  const auto ligand_direct = direct_gene_flags(ligandIdx);
  const auto receptor_direct = direct_gene_flags(receptorIdx);
  const int KK = K * K;
  const R_xlen_t KK_x = static_cast<R_xlen_t>(KK);
  const double Kh_n = std::pow(Kh, n_power);
  const int requested_threads = std::max(1, threads);
  int threads_used = 1;
#ifdef _OPENMP
  threads_used = std::min(requested_threads, std::max(1, nboot));
#endif


  const char* target_path=std::getenv("CELLCHATESS_TRACE_TARGETS");
  const char* output_path=std::getenv("CELLCHATESS_TRACE_OUTPUT");
  std::vector<std::pair<int,int> > trace_targets;
  if(target_path && target_path[0]) {
    std::ifstream input(target_path); if(!input)stop("Cannot open trace targets");
    int lr0,pair0; while(input>>lr0>>pair0)trace_targets.emplace_back(lr0,pair0);
  }
  std::ofstream trace_output;
  if(!trace_targets.empty()) {
    if(threads_used!=1)stop("Diagnostic tracing requires one thread");
    if(!output_path || !output_path[0])stop("Trace output path missing");
    trace_output.open(output_path,std::ios::app);if(!trace_output)stop("Cannot open trace output");
  }
  auto trace_decision=[&](int lr,int idx,int b,double observed,double boot) {
    if(std::find(trace_targets.begin(),trace_targets.end(),std::make_pair(lr,idx))==trace_targets.end())return;
    char line[1024];
    std::snprintf(line,sizeof(line),"%d\t%d\t%d\t%.17g\t%.17g\t%.17g\t%d\t%a\t%a\n",lr+1,idx,b+1,observed,boot,boot-observed,(boot-observed>0.0)?1:0,observed,boot);
    trace_output<<line;
  };

  std::vector< std::vector<int> > boot_cells(nboot * K);
  for (int b = 0; b < nboot; ++b) {
    for (int c = 0; c < C; ++c) {
      const int k = groupBoot(c, b) - 1;
      if (k >= 0 && k < K) boot_cells[b * K + k].push_back(c);
    }
  }

  std::vector<int> cache_genes;
  cache_genes.reserve(G);
  for (int lr = 0; lr < nLR; ++lr) {
    add_index_genes(cache_genes, ligandIdx, lr);
    add_index_genes(cache_genes, receptorIdx, lr);
    add_index_genes(cache_genes, coAIdx, lr);
    add_index_genes(cache_genes, coIIdx, lr);
    if (hasAgonist[lr]) add_index_genes(cache_genes, agonistIdx, lr);
    if (hasAntagonist[lr]) add_index_genes(cache_genes, antagonistIdx, lr);
  }
  std::sort(cache_genes.begin(), cache_genes.end());
  cache_genes.erase(std::unique(cache_genes.begin(), cache_genes.end()), cache_genes.end());
  const int U = cache_genes.size();
  std::vector<int> gene_to_cache(G, -1);
  for (int u = 0; u < U; ++u) gene_to_cache[cache_genes[u]] = u;

  // Layout: [bootstrap][gene][group]. The innermost group dimension is
  // contiguous, matching the ESS access pattern for one LR gene.
  std::vector<double> tri_cache(static_cast<size_t>(U) * K * nboot, NA_REAL);
  std::vector<unsigned char> tri_seen(static_cast<size_t>(U) * K * nboot, 0);

  NumericVector prob(Dimension(K, K, nLR));
  NumericVector pval(Dimension(K, K, nLR));
  std::vector<double> L(K), R(K), ago(K), ant(K);
  // Reuse serial scratch buffers across LR interactions.
  std::vector<double> Lb_serial(K), Rb_serial(K), agoB_serial(K), antB_serial(K);
  std::vector<int> gene_to_local(G, -1);

  long long active_pairs = 0;
  long long skipped_pairs = 0;
  long long boot_tri_mean_evals = 0;
  long long cache_hits = 0;
  long long streamed_lr = 0;
  long long max_lr_genes = 0;

  for (int lr = 0; lr < nLR; ++lr) {
    const bool simple_lr = !hasAgonist[lr] && !hasAntagonist[lr] &&
      index_row_empty(coAIdx, lr) && index_row_empty(coIIdx, lr);
    const R_xlen_t lr_offset = KK_x * lr;

    for (int idx = 0; idx < KK; ++idx) {
      pval[static_cast<R_xlen_t>(idx) + lr_offset] = 1.0;
    }

    for (int k = 0; k < K; ++k) {
      L[k] = avg_expr_matrix(avg, ligandIdx, lr, k, ligand_direct[lr]);
      R[k] = avg_expr_matrix(avg, receptorIdx, lr, k, receptor_direct[lr]);
      if (simple_lr) {
        ago[k] = 1.0;
        ant[k] = 1.0;
      } else {
        const double coA = coreceptor_factor_matrix(avg, coAIdx, lr, k);
        const double coI = coreceptor_factor_matrix(avg, coIIdx, lr, k);
        R[k] = R[k] * coA / coI;
        ago[k] = hasAgonist[lr] ? cofactor_factor_matrix(avg, agonistIdx, lr, k, true, Kh_n, n_power) : 1.0;
        ant[k] = hasAntagonist[lr] ? cofactor_factor_matrix(avg, antagonistIdx, lr, k, false, Kh_n, n_power) : 1.0;
      }
    }

    std::vector<double> pnull(KK, 0.0);
    std::vector<int> active_idx;
    std::vector<int> active_groups;
    std::vector<int> active_senders;
    std::vector<int> active_receivers;
    std::vector<char> group_seen(K, 0);
    active_senders.reserve(K);
    active_receivers.reserve(K);
    active_groups.reserve(K);

    for (int k = 0; k < K; ++k) {
      if (L[k] > 0.0) active_senders.push_back(k);
      if (R[k] > 0.0) active_receivers.push_back(k);
    }
    active_idx.reserve(active_senders.size() * active_receivers.size());

    for (int s : active_senders) {
      for (int r : active_receivers) {
        const double dataLR = L[s] * R[r];
        const double dataLR_n = std::pow(dataLR, n_power);
        const double p1 = dataLR_n / (Kh_n + dataLR_n);
        const double val = p1 * (ago[s] * ago[r]) * (ant[s] * ant[r]);
        const int idx = s + K * r;
        pnull[idx] = val;
        prob[static_cast<R_xlen_t>(idx) + lr_offset] = val;
        if (val != 0.0) {
          active_idx.push_back(idx);
          if (!group_seen[s]) {
            group_seen[s] = 1;
            active_groups.push_back(s);
          }
          if (!group_seen[r]) {
            group_seen[r] = 1;
            active_groups.push_back(r);
          }
        }
      }
    }

    active_pairs += static_cast<long long>(active_idx.size());
    skipped_pairs += static_cast<long long>(KK - active_idx.size());
    if (active_idx.empty()) continue;
    std::sort(active_groups.begin(), active_groups.end());
    streamed_lr += 1;

    std::vector<int> lr_genes;
    add_index_genes(lr_genes, ligandIdx, lr);
    add_index_genes(lr_genes, receptorIdx, lr);
    if (!simple_lr) {
      add_index_genes(lr_genes, coAIdx, lr);
      add_index_genes(lr_genes, coIIdx, lr);
      if (hasAgonist[lr]) add_index_genes(lr_genes, agonistIdx, lr);
      if (hasAntagonist[lr]) add_index_genes(lr_genes, antagonistIdx, lr);
    }
    std::sort(lr_genes.begin(), lr_genes.end());
    lr_genes.erase(std::unique(lr_genes.begin(), lr_genes.end()), lr_genes.end());
    max_lr_genes = std::max(max_lr_genes, static_cast<long long>(lr_genes.size()));
    for (int i = 0; i < static_cast<int>(lr_genes.size()); ++i) {
      gene_to_local[lr_genes[i]] = i;
    }

    std::vector<int> reject(active_idx.size(), 0);
    long long lr_boot_tri_mean_evals = 0;
    long long lr_cache_hits = 0;

    auto process_bootstrap = [&](
      const int b,
      std::vector<double>& bootAvg,
      std::vector<double>& Lb,
      std::vector<double>& Rb,
      std::vector<double>& agoB,
      std::vector<double>& antB,
      std::vector<int>& reject_local,
      long long& tri_mean_evals_local,
      long long& cache_hits_local
    ) {
      std::fill(bootAvg.begin(), bootAvg.end(), NA_REAL);
      for (int local = 0; local < static_cast<int>(lr_genes.size()); ++local) {
        const int gene0 = lr_genes[local];
        const int cache_gene = gene_to_cache[gene0];
        if (cache_gene < 0) continue;
        const size_t cache_base = static_cast<size_t>(K) *
          (static_cast<size_t>(cache_gene) + static_cast<size_t>(U) * static_cast<size_t>(b));
        for (int k : active_groups) {
          const size_t cache_idx = cache_base + static_cast<size_t>(k);
          if (!tri_seen[cache_idx]) {
            tri_cache[cache_idx] = tri_mean_gene_cells(data, gene0, boot_cells[b * K + k]);
            tri_seen[cache_idx] = 1;
            tri_mean_evals_local += 1;
          } else {
            cache_hits_local += 1;
          }
          bootAvg[local * K + k] = tri_cache[cache_idx];
        }
      }

      for (int k : active_groups) {
        Lb[k] = avg_expr_local(bootAvg, gene_to_local, ligandIdx, lr, k, K, ligand_direct[lr]);
        Rb[k] = avg_expr_local(bootAvg, gene_to_local, receptorIdx, lr, k, K, receptor_direct[lr]);
        if (simple_lr) {
          agoB[k] = 1.0;
          antB[k] = 1.0;
        } else {
          const double coA = coreceptor_factor_local(bootAvg, gene_to_local, coAIdx, lr, k, K);
          const double coI = coreceptor_factor_local(bootAvg, gene_to_local, coIIdx, lr, k, K);
          Rb[k] = Rb[k] * coA / coI;
          agoB[k] = hasAgonist[lr] ? cofactor_factor_local(bootAvg, gene_to_local, agonistIdx, lr, k, K, true, Kh_n, n_power) : 1.0;
          antB[k] = hasAntagonist[lr] ? cofactor_factor_local(bootAvg, gene_to_local, antagonistIdx, lr, k, K, false, Kh_n, n_power) : 1.0;
        }
      }

      for (int j = 0; j < static_cast<int>(active_idx.size()); ++j) {
        const int idx = active_idx[j];
        const int s = idx % K;
        const int r = idx / K;
        const double dataLR = Lb[s] * Rb[r];
        double p1 = 0.0;
        if (dataLR > 0.0) {
          const double dataLR_n = std::pow(dataLR, n_power);
          p1 = dataLR_n / (Kh_n + dataLR_n);
        }
        const double val = p1 * (agoB[s] * agoB[r]) * (antB[s] * antB[r]);
        trace_decision(lr,idx,b,pnull[idx],val);
        if (val - pnull[idx] > 0.0) reject_local[j] += 1;
      }
    };

    if (threads_used == 1) {
      // Avoid OpenMP setup overhead for the default single-thread path.
      std::vector<double> bootAvg(lr_genes.size() * K, NA_REAL);
      // bootAvg is LR-sized; the K-sized buffers are reused across LRs.
      for (int b = 0; b < nboot; ++b) {
        process_bootstrap(
          b, bootAvg, Lb_serial, Rb_serial, agoB_serial, antB_serial, reject,
          lr_boot_tri_mean_evals, lr_cache_hits
        );
      }
    } else {
      std::vector< std::vector<int> > reject_by_thread(
        threads_used, std::vector<int>(active_idx.size(), 0)
      );

#ifdef _OPENMP
#pragma omp parallel num_threads(threads_used) reduction(+:lr_boot_tri_mean_evals,lr_cache_hits)
#endif
      {
        int thread_id = 0;
#ifdef _OPENMP
        thread_id = omp_get_thread_num();
#endif
        std::vector<int>& reject_local = reject_by_thread[thread_id];
        std::vector<double> bootAvg(lr_genes.size() * K, NA_REAL);
        std::vector<double> Lb(K), Rb(K), agoB(K), antB(K);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (int b = 0; b < nboot; ++b) {
          process_bootstrap(
            b, bootAvg, Lb, Rb, agoB, antB, reject_local,
            lr_boot_tri_mean_evals, lr_cache_hits
          );
        }
      }

      for (int thread_id = 0; thread_id < threads_used; ++thread_id) {
        for (int j = 0; j < static_cast<int>(active_idx.size()); ++j) {
          reject[j] += reject_by_thread[thread_id][j];
        }
      }
    }

    boot_tri_mean_evals += lr_boot_tri_mean_evals;
    cache_hits += lr_cache_hits;

    for (int j = 0; j < static_cast<int>(active_idx.size()); ++j) {
      const int idx = active_idx[j];
      pval[static_cast<R_xlen_t>(idx) + lr_offset] = static_cast<double>(reject[j]) / static_cast<double>(nboot);
    }

    for (int g0 : lr_genes) {
      gene_to_local[g0] = -1;
    }
  }

  const long long total_pairs = static_cast<long long>(nLR) * static_cast<long long>(KK);
  return List::create(
    Named("prob") = prob,
    Named("pval") = pval,
    Named("active_pairs") = static_cast<double>(active_pairs),
    Named("skipped_pairs") = static_cast<double>(skipped_pairs),
    Named("total_pairs") = static_cast<double>(total_pairs),
    Named("active_fraction") = total_pairs > 0 ? static_cast<double>(active_pairs) / static_cast<double>(total_pairs) : NA_REAL,
    Named("boot_tri_mean_evals") = static_cast<double>(boot_tri_mean_evals),
    Named("cache_hits") = static_cast<double>(cache_hits),
    Named("cache_genes") = static_cast<double>(U),
    Named("cache_slots") = static_cast<double>(static_cast<long long>(U) * static_cast<long long>(K) * static_cast<long long>(nboot)),
    Named("streamed_lr") = static_cast<double>(streamed_lr),
    Named("max_lr_genes") = static_cast<double>(max_lr_genes),
    Named("threads_used") = threads_used
  );
}

// [[Rcpp::export]]
List pathway_sum_cpp(NumericVector prob, NumericVector pval, IntegerVector pathway, int n_pathways, double thresh) {
  const IntegerVector dims = prob.attr("dim");
  const int K = dims[0];
  const int nLR = dims[2];
  const int KK = K * K;
  const R_xlen_t KK_x = static_cast<R_xlen_t>(KK);
  NumericVector prob_pathway(Dimension(K, K, n_pathways));
  NumericVector lr_sum(nLR);
  NumericVector pathway_sum(n_pathways);

  for (int lr = 0; lr < nLR; ++lr) {
    const int p = pathway[lr] - 1;
    if (p < 0 || p >= n_pathways) continue;
    double lr_total = 0.0;
    for (int idx = 0; idx < KK; ++idx) {
      const int src_tgt = idx;
      const R_xlen_t pos = static_cast<R_xlen_t>(src_tgt) + KK_x * lr;
      double val = prob[pos];
      if (pval[pos] > thresh) val = 0.0;
      prob_pathway[static_cast<R_xlen_t>(src_tgt) + KK_x * p] += val;
      lr_total += val;
      pathway_sum[p] += val;
    }
    lr_sum[lr] = lr_total;
  }

  return List::create(
    Named("prob_pathway") = prob_pathway,
    Named("lr_sum") = lr_sum,
    Named("pathway_sum") = pathway_sum
  );
}

// [[Rcpp::export]]
List aggregate_net_cpp(NumericVector prob, NumericVector pval, double thresh) {
  const IntegerVector dims = prob.attr("dim");
  const int K = dims[0];
  const int nLR = dims[2];
  const int KK = K * K;
  const R_xlen_t KK_x = static_cast<R_xlen_t>(KK);
  NumericMatrix count(K, K);
  NumericMatrix weight(K, K);

  for (int lr = 0; lr < nLR; ++lr) {
    for (int idx = 0; idx < KK; ++idx) {
      const R_xlen_t pos = static_cast<R_xlen_t>(idx) + KK_x * lr;
      const double val = prob[pos];
      double pv = pval[pos];
      if (val == 0.0) pv = 1.0;
      if (pv < thresh && val > 0.0) {
        const int s = idx % K;
        const int r = idx / K;
        count(s, r) += 1.0;
        weight(s, r) += val;
      }
    }
  }

  return List::create(Named("count") = count, Named("weight") = weight);
}
