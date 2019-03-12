
#include "flip.h"

namespace ld_estimator {

template <typename T>
void fliphap(Haps<T>& hap) {
  // flip AA with AB and BA with BB
  T aa = hap.aa;
  T ab = hap.ab;
  T ba = hap.ba;
  T bb = hap.bb;

  // flip the values
  hap.ab = aa;
  hap.aa = ab;
  hap.ba = bb;
  hap.bb = ba;
}

void flip_freqs(Haps<double>& freqs, Haps<int>& known, double& pA2, double& pB2, double& d) {
  //flip frequencies to get positive D'
  fliphap(freqs);
  fliphap(known);
  d = get_d(freqs);

  // flip frequency of second allele
  std::swap(pA2, pB2);
}

}
