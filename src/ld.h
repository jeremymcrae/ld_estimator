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

Linkage pairwise(char ** a1s, char ** a2s, char ** b1s, char ** b2s, int size, std::vector<bool> & ploidy);

}

#endif
