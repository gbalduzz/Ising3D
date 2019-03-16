#include "ising_lattice.hpp"

#include <cassert>

IsingLattice::IsingLattice(int L)
    : L_(L), n_(L * L * L), spins_(n_), rng_(0), E_(0), M_(0) {
  std::uniform_int_distribution<std::int8_t> distro(0, 1);
  for (int i = 0; i < n_; ++i) {
    spins_[i] = 2 * distro(rng_) - 1;
  }
  computeEandM();
}

void IsingLattice::computeEandM() {
  E_ = M_ = 0;
  unsigned int linindex = 0;
  for (int k = 0; k < L_; ++k)
    for (int j = 0; j < L_; ++j)
      for (int i = 0; i < L_; ++i) {
        M_ += spins_[linindex];
        E_ -= spins_[linindex] * rightMagnetization(i, j, k);
        ++linindex;
      }
}

int IsingLattice::haloMagnetization(int idx) const {
  // Take care of periodic boundary conditions
  // Note: don't exclude that implementing this in terms of '%' operator is
  // faster on certain architectures.
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

  return spins_[index(i + 1, j, k)] + spins_[index(i - 1, j, k)] +
         spins_[index(i, j + 1, k)] + spins_[index(i, j - 1, k)] +
         spins_[index(i, j, k + 1)] + spins_[index(i, j, k - 1)];
}

int IsingLattice::rightMagnetization(int i, int j, int k) const {
  // TODO: Remove copy.
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

  return spins_[index(i + 1, j, k)] + spins_[index(i, j + 1, k)] +
         spins_[index(i, j, k + 1)];
}
