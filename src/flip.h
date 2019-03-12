#ifndef LD_ESTIMATOR_FLIP_H
#define LD_ESTIMATOR_FLIP_H

#include <utility>

#include "dprime.h"
#include "haps.h"

namespace ld_estimator {
  template <typename T>
  void fliphap(Haps<T>& hap);
  void flip_freqs(Haps<double>& freqs, Haps<int>& known, double& pA2, double& pB2, double& d);
}

#endif
