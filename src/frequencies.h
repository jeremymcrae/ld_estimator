#ifndef LD_ESTIMATOR_FREQUENCIES_H
#define LD_ESTIMATOR_FREQUENCIES_H

#include <algorithm>

#include "haps.h"

namespace ld_estimator {
  Haps<double> count_haplotypes(int iteration, Haps<int> known, Haps<double> freqs, int unknown);
  Haps<double> estimate_frequencies(Haps<double> counts, double prob);
}

#endif
