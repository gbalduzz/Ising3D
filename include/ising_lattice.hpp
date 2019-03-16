// Base lattice class.
#pragma once

#include <array>
#include <random>
#include <vector>

#include "typedefs.hpp"

class IsingLattice {
public:
  IsingLattice(int L);

  virtual void doSweep() = 0;
  virtual void setBeta(Real beta) = 0;
  virtual void markThermalized(){}

  Real getE() const { return Real(E_) / n_; } // Return energy density.
  Real getM() const { return Real(M_) / n_; } // Return magnetization density.
  auto size() const { return n_; }

protected:
  template <bool right_only = false> int haloMagnetization(int idx) const;
  void computeEandM();

  const int L_;
  const int n_;
  Real beta_;
  std::vector<std::int8_t> spins_;

  std::mt19937_64 rng_;

  long int E_; // in units of J.
  long int M_; // not normalized.
};
