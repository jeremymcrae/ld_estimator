#ifndef LD_ESTIMATOR_AF_H
#define LD_ESTIMATOR_AF_H

#include <vector>

#include "haps.h"

namespace ld_estimator {

struct Freqs {
  double A1;
  double B1;
  double A2;
  double B2;
};

Freqs get_allele_freqs(Haps<int> known, double unknown);

}

#endif
