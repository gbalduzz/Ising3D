#include "time_series.hpp"

// TODO: unroll summation loops for efficiency, if necessary.

Real mean(const std::vector<Real> &v) {
  double sum = 0;
  for (Real x : v)
    sum += x;

  return sum / v.size();
}

Real var2(const std::vector<Real> &v, Real mean) {
  double sum = 0;
  for (Real x : v)
    sum += (x - mean) * (x - mean);

  return sum / v.size();
}

Real var2(const std::vector<Real> &v) { return var2(v, mean(v)); }

std::vector<Real> autocorrelation(const std::vector<Real> &v, int max_tau) {
  const Real mu = mean(v);
  std::vector<Real> result(max_tau, 0);

  for (int i = 0; i < v.size() - max_tau; ++i) {
    const Real diff_t0 = (v[i] - mu);
    for (int j = 0; j < max_tau; ++j)
      result[j] += diff_t0 * (v[i + j] - mu);
  }

  const Real norm = 1. / result[0]; // results[0] = n * sigma2.
  for (Real &x : result)
    x *= norm;

  return result;
}

Real integratedAutocorrelation(const std::vector<Real> &autocorr_func) {
  Real tau = autocorr_func[0];

  for (int i = 1; i < autocorr_func.size(); ++i) {
    if (autocorr_func[i] < 0 || autocorr_func[i] > autocorr_func[i - 1])
      break;
    tau += autocorr_func[i];
  }

  return tau;
}
