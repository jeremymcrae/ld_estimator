#ifndef LD_ESTIMATOR_HAPS_H
#define LD_ESTIMATOR_HAPS_H

namespace ld_estimator {
  template <typename T>
  struct Haps {
    // store haplotype counts or haplotype frequencies
    T aa;
    T ab;
    T ba;
    T bb;
  };
}

#endif
