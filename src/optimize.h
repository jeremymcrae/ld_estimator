#ifndef LD_ESTIMATOR_OPTIMIZE_H
#define LD_ESTIMATOR_OPTIMIZE_H

#include <cmath>

#include "likelihoods.h"
#include "frequencies.h"
#include "haps.h"

namespace ld_estimator {
  Haps<double> get_frequencies(Haps<int> known, int unknown, double epsilon=0.00000001);
}

#endif
