#ifndef LD_ESTIMATOR_LINKAGE_H
#define LD_ESTIMATOR_LINKAGE_H

#include <string>

#include "haps.h"

namespace ld_estimator {

struct Phase {
  std::string var1_allele;
  std::string var2_allele;
};

struct Linkage {
  double dprime;
  double loglikelihood;
  double r_squared;
  Haps<double> freqs;
  Phase phase;
};
}

#endif
