#pragma once

#include <vector>

#include "typedefs.hpp"

Real mean(const std::vector<Real> &v);

Real var2(const std::vector<Real> &v, Real mean);
Real var2(const std::vector<Real> &v);

std::vector<Real> autocorrelation(const std::vector<Real> &v, int max_tau);
Real integratedAutocorrelation(const std::vector<Real> &autocorr_func);
