#ifndef LD_ESTIMATOR_LIKELIHOODS_H
#define LD_ESTIMATOR_LIKELIHOODS_H

#include <vector>
#include <cmath>

#include "haps.h"

namespace ld_estimator {
  static const double LN10 = std::log(10.0);
  
  double likelihood_freqs(Haps<int> known, Haps<double> freqs, int unknown);
  double likelihood_null(Haps<int> known, double a1, double a2, double b1, double b2,
    int unknown);
  std::vector<double> get_lsurface(double pA1, double pA2, double pB1,
    double denom, Haps<int> known, int unknown);
}

#endif
