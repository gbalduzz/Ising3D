#include "ising_lattice.hpp"

#include <cassert>

IsingLattice::IsingLattice(int L)
    : L_(L), n_(L * L * L), spins_(n_), rng_(0), E_(0), M_(0) {
  std::uniform_int_distribution<unsigned short> distro(0, 1);

  for (int i = 0; i < n_; ++i) {
    spins_[i] = distro(rng_);
    M_ += 2 * spins_[i] - 1;
  }

  E_ = computeE();
}

long int IsingLattice::computeE() const {
  long int E = 0;
  for (int i = 0; i < n_; ++i) {
    const int spin = spins_[i] ? 1 : -1;
    E += -spin * haloMagnetization(i);
  }

  assert(E >= -6 * n_ && E <= 6 * n_);
  assert((E % 2) == 0);
  return E / 2;
}

void IsingLattice::doSweep() {
  for (int i = 0; i < n_; ++i)
    doStep();
}

int IsingLattice::haloMagnetization(int idx) const {
  // Take care of periodic boundary conditions
  // Note: don't exclude that implementing this in terms of '%' operator is faster on certain architectures.
  auto pbs = [&](const int i) {
    if (i >= 0 && i < L_)
      return i;
    else if (i >= L_)
      return i - L_;
    else
      return i + L_;
  };
  auto index = [&](int i, int j, int k) {
    return pbs(i) + L_ * pbs(j) + L_ * L_ * pbs(k);
  };

  // Compute linear indices.
  const int k = idx / (L_ * L_);
  idx -= k * (L_ * L_);
  const int j = idx / L_;
  const int i = idx - L_ * j;

  const int n_up = spins_[index(i + 1, j, k)] + spins_[index(i - 1, j, k)] +
                   spins_[index(i, j + 1, k)] + spins_[index(i, j - 1, k)] +
                   spins_[index(i, j, k + 1)] + spins_[index(i, j, k - 1)];

  assert(n_up <= 6 && n_up >= 0);
  return 2 * n_up - 6;
}

void IsingLattice::setBeta(Real beta) {
  beta_ = beta;

  // Set table.
  for (int halo = 0; halo < exp_table_.size(); ++halo) {
    exp_table_[halo] = std::exp(-2 * beta_ * halo);
  }
}

Real IsingLattice::computeProb(const int halo) const {
  if (halo <= 0)
    return 1;
  assert(halo <= 6);
  return exp_table_[halo];
}

void IsingLattice::doStep() {
  std::uniform_int_distribution<int> distro_int(0, n_ - 1);
  std::uniform_real_distribution<Real> distro_real(0, 1);
  const int candidate = distro_int(rng_);

  const int halo = haloMagnetization(candidate);
  const int s_old = spins_[candidate] ? 1 : -1;

  const Real prob = computeProb(s_old * halo);
  const bool accept = distro_real(rng_) < prob;

  if (accept) {
    spins_[candidate] = !spins_[candidate]; // flip spin.
    const int delta_E = 2 * s_old * halo;
    M_ += -2 * s_old;
    E_ += delta_E;

    // Note: these are built in tests, they are disabled if you build in Release mode, which defines the flag NDEBUG.
    assert(E_ == computeE());
    assert(prob == std::min(std::exp(-beta_ * delta_E), Real(1)));
  }
}
