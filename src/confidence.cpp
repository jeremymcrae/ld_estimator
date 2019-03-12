
#include "confidence.h"

namespace ld_estimator {

std::vector<double> ld_confidence_interval(std::vector<double> lsurface,
    double loglike1, double alpha=0.05) {
  // estimate confidence interval for D'
  double total = 0.0;
  int size = lsurface.size();
  for (auto i=0; i < size; i++) {
    lsurface[i] -= loglike1;
    lsurface[i] = std::pow(10, lsurface[i]);
    total += lsurface[i];
  }

  int low_i = 0;
  double summed = 0.0;
  for (auto i=0; i < size; i++) {
    summed += lsurface[i];
    if (summed > alpha * total and (summed - lsurface[i]) < alpha * total) {
      low_i = i - 1;
      break;
    }
  }

  int high_i = 0;
  summed = 0.0;
  for (auto i=99; i >= 0; i--) {
    summed += lsurface[i];
    if (summed > alpha * total and (summed - lsurface[i]) < alpha * total) {
      high_i = i + 1;
      break;
    }
  }

  if (high_i > 100) {
    high_i = 100;
  }

  return std::vector<double> {(double)low_i / 100, (double)high_i / 100};
}

}
