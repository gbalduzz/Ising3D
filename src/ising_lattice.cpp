#include "ising_lattice.hpp"

#include <cassert>

IsingLattice::IsingLattice(int L)
    : L_(L), n_(L * L * L), spins_(n_), rng_(0), E_(0), M_(0) {
  std::uniform_int_distribution<std::int8_t> distro(0, 1);

  for (int i = 0; i < n_; ++i) {
    spins_[i] = 2 * distro(rng_) - 1;
    M_ += spins_[i];
  }

  E_ = computeE();
}

long int IsingLattice::computeE() const {
  long int E = 0;
  for (int i = 0; i < n_; ++i) {
    E += -spins_[i] * haloMagnetization<true>(i);
  }
  return E;
}

template <bool right_only> int IsingLattice::haloMagnetization(int idx) const {
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

  if
    constexpr(right_only) {
      return spins_[index(i + 1, j, k)] + spins_[index(i, j + 1, k)] +
             spins_[index(i, j, k + 1)];
    }
  else {
    return spins_[index(i + 1, j, k)] + spins_[index(i - 1, j, k)] +
           spins_[index(i, j + 1, k)] + spins_[index(i, j - 1, k)] +
           spins_[index(i, j, k + 1)] + spins_[index(i, j, k - 1)];
  }
}

template int IsingLattice::haloMagnetization<true>(int idx) const;
template int IsingLattice::haloMagnetization<false>(int idx) const;
