#include "ising_lattice.hpp"

IsingLattice::IsingLattice(int L, Real beta)
    : L_(L), n_(L * L * L), beta_(beta), spins_(n_), rng_(0), E_(0), M_(0) {
  std::uniform_int_distribution<unsigned short> distro(0, 1);

  for (int i = 0; i < n_; ++i)
    spins_[i] = distro(rng_);

  // Initialize energy and magnetization.
  for (int i = 0; i < n_; ++i) {
    const int spin = spins_[i] ? 1 : -1;
    M_ += spin;
    E_ += -spin * haloMagnetization(i);
  }
}

void IsingLattice::doSweep() {
  for (int i = 0; i < n_; ++i)
    doStep();
}

int IsingLattice::haloMagnetization(int idx) const {
  // Take care of periodic boundary conditions
  auto pbs = [&](const int i) {
    if (i > 0 && i < L_)
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

  return 2 * n_up - 6;
}

void IsingLattice::doStep() {
  std::uniform_int_distribution<int> distro_int(0, n_ - 1);
  std::uniform_real_distribution<Real> distro_real(0, 1);
  const int candidate = distro_int(rng_);

  const int halo = haloMagnetization(candidate);
  const int delta_s = spins_[candidate] ? -2 : 2;

  // prob = Exp[ -beta delta_spin].
  // TODO: store in table;
  const Real prob = std::exp(-beta_ * (delta_s * halo));
  const bool accept = distro_real(rng_) < prob;

  if (accept) {
    spins_[candidate] = !spins_[candidate]; // flip spin.
    M_ += delta_s;
    E_ += delta_s * halo;
  }
}
