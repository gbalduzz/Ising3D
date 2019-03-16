#pragma once

#include <array>

#include "ising_lattice.hpp"
#include "typedefs.hpp"

class Metropolis final : public IsingLattice {
public:
  Metropolis(int L);

  void doSweep() override;

  void setBeta(Real beta) override;

private:
  void doStep();
  Real computeProb(int halo) const;

  std::array<Real, 7> exp_table_;
};
