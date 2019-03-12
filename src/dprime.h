#ifndef LD_ESTIMATOR_DPRIME_H
#define LD_ESTIMATOR_DPRIME_H

#include <algorithm>

#include "haps.h"

namespace ld_estimator {
  double get_d(Haps<double> freqs);
  double get_denominator(Haps<double> freqs);
}

#endif
