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
  Linkage pairwise(std::vector<std::vector<std::string> > var1,
    std::vector<std::vector<std::string> > var2,
    std::vector<bool> ploidy);
}

#endif
