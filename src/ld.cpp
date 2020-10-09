
#include "ld.h"

namespace ld_estimator {

Linkage pairwise(
    std::vector<std::string> & var_a1, std::vector<std::string> & var_a2,
    std::vector<std::string> & var_b1, std::vector<std::string> & var_b2,
    std::vector<bool> & ploidy) {
  // compute linkage disequilibrium between two variants
  //
  // Args:
  //     var1: list of genotypes for first variant
  //     var2: list of genotypes for second variant, ordered as per var1
  //     ploidy: list of ploidy states, ordered as per var1
  if (is_monomorphic(var_a1, var_a2) or is_monomorphic(var_b1, var_b2)) {
    throw std::invalid_argument("monomorphic variant");
  }
  
  std::array<std::string, 2> alleles1;
  std::array<std::string, 2> alleles2;
  get_alleles(var_a1, var_a2, alleles1);
  get_alleles(var_b1, var_b2, alleles2);
  Phase phase = {alleles1[0], alleles2[0]};
  std::pair<Haps<int>, int> counts = tally_haplotypes(var_a1, var_a2, var_b1, var_b2, ploidy, alleles1[0], alleles1[1], alleles2[0], alleles1[1]);
  Haps<int> known = counts.first;
  int unknown = counts.second;
  if (lacks_haplotypes(known, unknown)) {
    Haps<double> freqs = {0.0, 0.0, 0.0, 0.0};
    return Linkage {1.0, 0.0, 0.0, freqs, phase};
  }

  Freqs f = get_allele_freqs(known, unknown);

  Haps<double> freqs = get_frequencies(known, unknown, 0.00000001);
  double loglike1 = likelihood_freqs(known, freqs, unknown);
  double loglike0 = likelihood_null(known, f.A1, f.A2, f.B1, f.B2, unknown);

  double d = get_d(freqs);
  if (d < 0) {
    flip_freqs(freqs, known, f.A2, f.B2, d);
    phase = {alleles1[0], alleles2[1]};
  }
  double denom = get_denominator(freqs);

  double rsq = (d * d) / (f.A1 * (1.0 - f.A1) * f.A2 * (1.0 - f.A2));
  return Linkage {d / denom, loglike1 - loglike0, rsq, freqs, phase};
}

}
