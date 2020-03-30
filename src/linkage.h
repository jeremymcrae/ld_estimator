#ifndef LD_ESTIMATOR_LINKAGE_H
#define LD_ESTIMATOR_LINKAGE_H

#include <vector>
#include <string>

#include "haps.h"

namespace ld_estimator {
  struct Linkage {
    double dprime;
    double loglikelihood;
    double r_squared;
    Haps<double> freqs;
    std::vector<std::string> phase;
  };
}

#endif
