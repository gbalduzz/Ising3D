#include "metropolis.hpp"

#include <cassert>

Metropolis::Metropolis(int L)
    : IsingLattice(L) {}

void Metropolis::doSweep() {
  for (int i = 0; i < n_; ++i)
    doStep();
}

void Metropolis::setBeta(Real beta) {
  beta_ = beta;

  // Set table.
  for (int halo = 0; halo < exp_table_.size(); ++halo) {
    exp_table_[halo] = std::exp(-2 * beta_ * halo);
  }
}

Real Metropolis::computeProb(const int halo) const {
  if (halo <= 0)
    return 1;
  assert(halo <= 6);
  return exp_table_[halo];
}

void Metropolis::doStep() {
  std::uniform_int_distribution<int> distro_int(0, n_ - 1);
  std::uniform_real_distribution<Real> distro_real(0, 1);
  const int candidate = distro_int(rng_);

  const int halo = haloMagnetization(candidate);
  const int s_old = spins_[candidate];

  const Real prob = computeProb(s_old * halo);
  const bool accept = distro_real(rng_) < prob;

  if (accept) {
    spins_[candidate] *= -1; // flip spin.
    const int delta_E = 2 * s_old * halo;
    M_ += -2 * s_old;
    E_ += delta_E;

    // Note: these are built in tests, they are disabled if you build in Release mode, which defines the flag NDEBUG.
    assert(E_ == computeE());
    assert(prob == std::min(std::exp(-beta_ * delta_E), Real(1)));
  }
}
