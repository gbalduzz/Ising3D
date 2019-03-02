#pragma once

#include <random>
#include <vector>

#include "typedefs.hpp"

class IsingLattice {
public:
  IsingLattice(int L, Real beta);

  void doSweep();

  Real getE() const { return E_; }
  Real getM() const { return M_; }
  auto size() const { return n_; }

  void setBeta(Real beta) { beta_ = beta; }

private:
  void doStep();
  int haloMagnetization(int idx) const;

  const int n_;
  const int L_;
  Real beta_;
  std::mt19937_64 rng_;

  std::vector<bool> spins_;

  long int E_; // in units of J.
  long int M_; // not normalized.
};
