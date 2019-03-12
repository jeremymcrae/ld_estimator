
#include "af.h"

namespace ld_estimator {
  std::vector<double> get_allele_freqs(Haps<int> known, double unknown) {
    double chroms = (known.aa + known.ab + known.ba + known.bb + (2 * unknown));
    double pA1 = (known.aa + known.ab + unknown) / chroms;
    double pA2 = (known.aa + known.ba + unknown) / chroms;
    return std::vector<double> {pA1, 1.0 - pA1, pA2, 1.0 - pA2};
  }
}
