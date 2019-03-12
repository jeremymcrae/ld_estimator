
#include "likelihoods.h"

namespace ld_estimator {

double likelihood_freqs(Haps<int> known, Haps<double> freqs, int unknown) {
    // get loglikelihood for model from estimated haplotype frequencies
    double aa = freqs.aa;
    double ab = freqs.ab;
    double ba = freqs.ba;
    double bb = freqs.bb;

    return (known.aa * std::log(aa) +
      known.ab * std::log(ab) +
      known.ba * std::log(ba) +
      known.bb * std::log(bb) +
      unknown * std::log(aa * bb + ab * ba)) / LN10;
}

double likelihood_null(Haps<int> known, double a1, double a2, double b1, double b2, int unknown) {
  // get loglikelihood for null model from allele frequencies
  return (known.aa * std::log(a1 * a2) +
    known.ab * std::log(a1 * b2) +
    known.ba * std::log(b1 * a2) +
    known.bb * std::log(b1 * b2) +
    unknown * std::log(2 * a1 * a2 * b1 * b2)) / LN10;
}

std::vector<double> get_lsurface(double pA1, double pA2, double pB1,
    double denom, Haps<int> known, int unknown) {
  // get log likelihood surface
  std::vector<double> lsurface;
  for (int i=0; i < 100; i++) {
    double dpr = i * 0.01;
    double AA = dpr * denom + pA1 * pA2;
    double AB = pA1 - AA;
    double BA = pA2 - AA;
    double BB = pB1 - BA;
    Haps<double> freqs = Haps<double> {AA, AB, BA, BB};
    if (i == 100) {
      // one value will be 0
      freqs.aa = std::max(1e-10, freqs.aa);
      freqs.ab = std::max(1e-10, freqs.ab);
      freqs.ba = std::max(1e-10, freqs.ba);
      freqs.bb = std::max(1e-10, freqs.bb);
    }
    double likelihood = likelihood_freqs(known, freqs, unknown);
    lsurface.push_back(likelihood);
  }
  return lsurface;
}

}
