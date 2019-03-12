
#include "frequencies.h"

namespace ld_estimator {

Haps<double> count_haplotypes(int iteration, Haps<int> known, Haps<double> freqs, int unknown) {
  // count haplotypes
  Haps<double> counts = Haps<double> {(double)known.aa, (double)known.ab, (double)known.ba, (double)known.bb};
  if (iteration > 0) {
    double obligates = freqs.aa * freqs.bb;
    double hets = freqs.ab * freqs.ba;
    counts.aa += unknown * obligates / (obligates + hets);
    counts.bb += unknown * obligates / (obligates + hets);
    counts.ab += unknown * hets / (obligates + hets);
    counts.ba += unknown * hets / (obligates + hets);
  }

  return counts;
}

Haps<double> estimate_frequencies(Haps<double> counts, double prob) {
  // estimate haplotype frequencies from haplotype counts
  double total = (counts.aa + counts.ab + counts.ba + counts.bb) + (4.0 * prob);

  // estimate frequency of each haplotype, but keep them above 1e-10
  double aa = std::max(1e-10, ((counts.aa + prob) / total));
  double ab = std::max(1e-10, ((counts.ab + prob) / total));
  double ba = std::max(1e-10, ((counts.ba + prob) / total));
  double bb = std::max(1e-10, ((counts.bb + prob) / total));

  return Haps<double> {aa, ab, ba, bb};
}

}
