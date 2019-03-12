
#include "dprime.h"

namespace ld_estimator {

double get_d(Haps<double> freqs){
  return freqs.aa * freqs.bb - freqs.ab * freqs.ba;
}

double get_denominator(Haps<double> freqs) {
  double a = (freqs.aa + freqs.ba) * (freqs.ba + freqs.bb);
  double b = (freqs.aa + freqs.ab) * (freqs.ab + freqs.bb);
  return std::min(a, b);
}

}
