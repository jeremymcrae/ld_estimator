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
  std::pair<Haps<int>, int> tally_haplotypes(std::vector<std::vector<std::string>> &var1,
      std::vector<std::vector<std::string>> &var2, std::vector<bool> &ploidy, 
      std::string major_1, std::string minor_1, std::string major_2, std::string minor_2);
}

#endif
