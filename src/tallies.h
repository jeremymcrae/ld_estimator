#ifndef LD_ESTIMATOR_TALLIES_H
#define LD_ESTIMATOR_TALLIES_H

#include <vector>
#include <utility>
#include <string>
#include <map>

#include "utils.h"
#include "haps.h"

namespace ld_estimator {
  // template <typename T>
std::pair<Haps<int>, int> tally_haplotypes(char ** a1s, char ** a2s,
    char ** b1s, char ** b2s, int size, std::vector<bool> & ploidy,
    const char * major_1, const char * minor_1, const char * major_2, const char * minor_2);
}

#endif
