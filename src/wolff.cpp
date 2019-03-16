#include "wolff.hpp"

#include <cassert>
#include <iostream>

Wolff::Wolff(int L) : IsingLattice(L) {}

void Wolff::doSweep() {
  int current_steps;

  if (steps_per_sweep_) // system is thermalized.
    current_steps = steps_per_sweep_;
  else // Use current average cluster size to compute n_ / avg_cluster_size
    current_steps =
        n_steps_ ? std::ceil(Real(n_) / cumulative_cluster_size_ * n_steps_)
                 : 1;

  assert(current_steps > 0);
  for (int i = 0; i < current_steps; ++i) {
    cumulative_cluster_size_ += doStep();
  }

  n_steps_ += current_steps;
  if (steps_per_sweep_) { // Get ready for a measurement.
    computeEandM();
  }
}

void Wolff::setBeta(Real beta) {
  beta_ = beta;

  // Mark as not thermalized.
  steps_per_sweep_ = 0;
  cumulative_cluster_size_ = 0;
  n_steps_ = 0;

  prob_ = 1 - std::exp(-2 * beta_);
}

void Wolff::markThermalized() {
  steps_per_sweep_ = std::ceil(Real(n_) / cumulative_cluster_size_ * n_steps_);

  std::cout << "\tAverage cluster size: "
            << Real(cumulative_cluster_size_) / n_steps_ << "\n";
}

std::array<unsigned int, 3> Wolff::unrollIndex(unsigned int idx) {
  std::array<unsigned int, 3> unrolled;
  // Compute linear indices.
  unrolled[2] = idx / (L_ * L_);
  idx -= unrolled[2] * (L_ * L_);
  unrolled[1] = idx / L_;
  unrolled[0] = idx - L_ * unrolled[1];

  return unrolled;
}

auto Wolff::neighbours(const std::array<unsigned int, 3> &site) {
  // Take care of periodic boundary conditions
  // Note: don't exclude that implementing this in terms of '%' operator is
  // faster on certain architectures.
  auto pbs = [&](const int i) -> unsigned int {
    if (i >= 0 && i < L_)
      return i;
    else if (i >= L_)
      return i - L_;
    else
      return i + L_;
  };
  auto index = [&](int i, int j, int k) -> std::array<unsigned int, 3> {
    return std::array<unsigned int, 3>{pbs(i), pbs(j), pbs(k)};
  };

  const int i = site[0];
  const int j = site[1];
  const int k = site[2];

  return std::array<std::array<unsigned int, 3>, 6>{
      index(i + 1, j, k), index(i - 1, j, k), index(i, j + 1, k),
      index(i, j - 1, k), index(i, j, k + 1), index(i, j, k - 1)};
}

std::size_t Wolff::doStep() {
  std::uniform_int_distribution<unsigned int> distro_int(0, n_ - 1);
  std::uniform_real_distribution<Real> distro_real(0, 1);

  const auto root_idx = distro_int(rng_);
  const std::int8_t s_old = spins_[root_idx];
  const std::int8_t new_spin = -s_old;

  std::size_t cluster_size = 0;
  assert(cluster_queue_.empty());

  auto add_to_cluster = [&](const std::array<unsigned int, 3> &site,
                            unsigned int linindex) {
    assert(linindex == site[0] + site[1] * L_ + site[2] * L_ * L_);
    cluster_queue_.push(site);
    spins_[linindex] = new_spin;
    ++cluster_size;
  };

  add_to_cluster(unrollIndex(root_idx), root_idx);

  // Perform a breath first search of the lattice and flip cluster members on the flight.
  while (!cluster_queue_.empty()) {
    const auto site = cluster_queue_.front();
    cluster_queue_.pop();

    for (auto next : neighbours(site)) {
      const unsigned int linindex = next[0] + next[1] * L_ + next[2] * L_ * L_;
      if (spins_[linindex] == s_old && distro_real(rng_) < prob_)
        add_to_cluster(next, linindex);
    }
  }

  return cluster_size;
}
