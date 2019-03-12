#ifndef LD_ESTIMATOR_LD_H
#define LD_ESTIMATOR_LD_H

#include <vector>
#include <string>
#include <stdexcept>

#include "af.h"
#include "confidence.h"
#include "dprime.h"
#include "flip.h"
#include "haps.h"
#include "likelihoods.h"
#include "linkage.h"
#include "optimize.h"
#include "tallies.h"
#include "utils.h"
#include "linkage.h"

namespace ld_estimator {
  // template <typename T>
  // double estimate_ld(std::vector<std::vector<T> > var1, std::vector<std::vector<T> > var2,
  //   std::vector<bool> ploidy);
  double estimate_ld(std::vector<std::vector<std::string> > var1,
    std::vector<std::vector<std::string> > var2,
    std::vector<bool> ploidy);
}

#endif
