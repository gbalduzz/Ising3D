#pragma once

#include <array>
#include <random>
#include <vector>

#include "typedefs.hpp"

class IsingLattice {
public:
  IsingLattice(int L);

  void doSweep();

  Real getE() const { return Real(E_) / n_; } // Return energy density.
  Real getM() const { return Real(M_) / n_; }
  auto size() const { return n_; }

  void setBeta(Real beta);

private:
  void doStep();
  int haloMagnetization(int idx) const;
  long int computeE() const;
  Real computeProb(int halo) const;

  const int L_;
  const int n_;
  Real beta_;
  std::vector<bool> spins_;

  std::mt19937_64 rng_;

  long int E_; // in units of J.
  long int M_; // not normalized.

  std::array<Real, 7> exp_table_;
};
