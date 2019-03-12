#ifndef LD_ESTIMATOR_AF_H
#define LD_ESTIMATOR_AF_H

#include <vector>

#include "haps.h"

namespace ld_estimator {
  std::vector<double> get_allele_freqs(Haps<int> known, double unknown);
}

#endif
