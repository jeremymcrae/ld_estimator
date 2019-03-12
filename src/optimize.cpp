
#include "optimize.h"

namespace ld_estimator {

Haps<double> get_frequencies(Haps<int> known, int unknown, double epsilon) {
  // run expectation maximisation to get optimal haplotype frequencies
  Haps<double> initial = Haps<double> {0.01, 0.01, 0.01, 0.01};
  Haps<double> counts = count_haplotypes(0, known, initial, unknown);
  Haps<double> frequencies = estimate_frequencies(counts, 0.1);

  int iteration = 1;
  double current = -999999999.0;
  double previous = 5;
  while (std::abs(current - previous) > epsilon and iteration < 1000) {
    previous = current;
    counts = count_haplotypes(iteration, known, frequencies, unknown);
    current = likelihood_freqs(known, frequencies, unknown);
    frequencies = estimate_frequencies(counts, 0.0);
    iteration += 1;
  }

  return frequencies;
}

}
