#ifndef LD_ESTIMATOR_CONFIDENCE_H
#define LD_ESTIMATOR_CONFIDENCE_H

#include <cmath>
#include <utility>
#include <vector>

namespace ld_estimator {
  std::vector<double> ld_confidence_interval(std::vector<double> lsurface,
    double loglike1, double alpha);
}

#endif
