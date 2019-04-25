
#include "ld.h"

namespace ld_estimator {

Linkage pairwise(std::vector<std::vector<std::string> > var1,
    std::vector<std::vector<std::string> > var2,
    std::vector<bool> ploidy) {
  // compute linkage disequilibrium between two variants
  //
  // Args:
  //     var1: list of genotypes for first variant
  //     var2: list of genotypes for second variant, ordered as per var1
  //     ploidy: list of ploidy states, ordered as per var1
  if (is_monomorphic(var1) or is_monomorphic(var2)) {
    throw std::invalid_argument("monomorphic variant");
  }
  
  std::pair<Haps<int>, int> counts = tally_haplotypes(var1, var2, ploidy);
  Haps<int> known = counts.first;
  int unknown = counts.second;
  if (lacks_haplotypes(known, unknown)) {
    Haps<double> freqs = {0.0, 0.0, 0.0, 0.0};
    return Linkage {1.0, 0.0, 0.0, 0.0, 0.0, freqs};
  }

  std::vector<double> f = get_allele_freqs(known, unknown);
  double pA1 = f[0];
  double pB1 = f[1];
  double pA2 = f[2];
  double pB2 = f[3];

  Haps<double> freqs = get_frequencies(known, unknown, 0.00000001);
  double loglike1 = likelihood_freqs(known, freqs, unknown);
  double loglike0 = likelihood_null(known, pA1, pA2, pB1, pB2, unknown);

  double d = get_d(freqs);
  std::vector<std::string> alleles1 = get_alleles(var1);
  std::vector<std::string> alleles2 = get_alleles(var2);
  std::vector<std::string> phase = {alleles1[0], alleles2[0]};
  if (d < 0) {
    flip_freqs(freqs, known, pA2, pB2, d);
    phase = {alleles1[0], alleles2[1]};
  }
  double denom = get_denominator(freqs);

  double rsq = (d * d) / (pA1 * (1.0 - pA1) * pA2 * (1.0 - pA2));
  std::vector<double> surface = get_lsurface(pA1, pA2, pB1, denom, known, unknown);
  std::vector<double> ci = ld_confidence_interval(surface, loglike1, 0.05);
  auto low_ci = ci[0];
  auto high_ci = ci[1];
  return Linkage {d / denom, loglike1 - loglike0, rsq, low_ci, high_ci, freqs, phase};
}

}
