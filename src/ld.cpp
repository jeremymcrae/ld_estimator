
#include "ld.h"

namespace ld_estimator {

Linkage pairwise(char ** a1s, char ** a2s, char ** b1s, char ** b2s, int size, std::vector<bool> & ploidy) {
  // compute linkage disequilibrium between two variants
  //
  // Args:
  //     var1: list of genotypes for first variant
  //     var2: list of genotypes for second variant, ordered as per var1
  //     ploidy: list of ploidy states, ordered as per var1
  if (is_monomorphic(a1s, a2s, size) or is_monomorphic(b1s, b2s, size)) {
    throw std::invalid_argument("monomorphic variant");
  }
  
  std::array<const char *, 2> alleles1;
  std::array<const char *, 2> alleles2;
  get_alleles(a1s, a2s, size, ploidy, alleles1);
  get_alleles(b1s, b2s, size, ploidy, alleles2);
  
  Phase phase = {std::string(alleles1[0]), std::string(alleles2[0])};
  std::pair<Haps<int>, int> counts = tally_haplotypes(a1s, a2s, b1s, b2s, size, ploidy, alleles1[0], alleles1[1], alleles2[0], alleles1[1]);
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
