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
  bool is_monomorphic(std::vector<std::vector<std::string> > & var);
  bool lacks_haplotypes(Haps<int> & counts, int unknown);
  void get_alleles(std::vector<std::vector<std::string> > & var, std::array<std::string, 2> & alleles);
}

#endif
