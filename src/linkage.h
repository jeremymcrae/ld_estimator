#ifndef LD_ESTIMATOR_LINKAGE_H
#define LD_ESTIMATOR_LINKAGE_H

#include <vector>

#include "haps.h"

namespace ld_estimator {
  template <typename T>
  struct LD {
    double dprime;
    double loglikelihood;
    double r_squared;
    double ci_low;
    double ci_high;
    Haps<T> freqs;
  };
}

#endif
