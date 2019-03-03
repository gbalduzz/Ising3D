#pragma once

#include <random>
#include <vector>

#include "typedefs.hpp"

class IsingLattice {
public:
  IsingLattice(int L, Real beta);

  void doSweep();

  Real getE() const { return Real(E_) / n_; } // Return energy density.
  Real getM() const { return Real(M_) / n_; }
  auto size() const { return n_; }

  void setBeta(Real beta) { beta_ = beta; }

private:
  void doStep();
  int haloMagnetization(int idx) const;
  long int computeE() const;

  const int L_;
  const int n_;
  Real beta_;
  std::vector<short int> spins_;

  std::mt19937_64 rng_;

  long int E_; // in units of J.
  long int M_; // not normalized.
};
