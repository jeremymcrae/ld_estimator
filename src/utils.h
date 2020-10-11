#ifndef LD_ESTIMATOR_UTILS_H
#define LD_ESTIMATOR_UTILS_H

#include <array>
#include <vector>
#include <string>
#include <set>
#include <utility>
#include <algorithm>
#include <map>

#include "haps.h"

namespace ld_estimator {
  bool is_monomorphic(char ** a1, char ** a2, int size);
  bool lacks_haplotypes(Haps<int> & counts, int unknown);
  void get_alleles(char ** a1, char ** a2, int size, std::vector<bool> & haploid, std::array<const char *, 2> & alleles);
}

#endif
