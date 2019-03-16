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

void Wolff::addToCluster(unsigned int id) {
  cluster_queue_.push(id);
  cluster_members_.push_back(id);
  spins_[id] = 0; // Mark as part of the cluster.
}

std::array<unsigned int, 6> Wolff::neighbours(unsigned int idx) {
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
  auto index = [&](int i, int j, int k) -> unsigned int {
    return pbs(i) + L_ * pbs(j) + L_ * L_ * pbs(k);
  };

  // Compute linear indices.
  const  int k = idx / (L_ * L_);
  idx -= k * (L_ * L_);
  const int j = idx / L_;
  const int i = idx - L_ * j;

  const std::array<unsigned int, 6> nn{index(i + 1, j, k), index(i - 1, j, k),
                                     index(i, j + 1, k), index(i, j - 1, k),
                                     index(i, j, k + 1), index(i, j, k - 1)};
  for(auto el : nn)
      assert(el < n_);
  return nn;
}

std::size_t Wolff::doStep() {
  std::uniform_int_distribution<unsigned int> distro_int(0, n_ - 1);
  std::uniform_real_distribution<Real> distro_real(0, 1);
  const unsigned int candidate = distro_int(rng_);
  const auto s_old = spins_[candidate];

  cluster_members_.clear();
  assert(cluster_queue_.empty());

  addToCluster(candidate);

  while (!cluster_queue_.empty()) {
    const unsigned int site = cluster_queue_.front();
    cluster_queue_.pop();

    for (auto next : neighbours(site)) {
      if (spins_[next] == s_old && distro_real(rng_) < prob_)
        addToCluster(next);
    }
  }

  for (unsigned int idx : cluster_members_) {
    spins_[idx] = -s_old; // flip cluster spins.
  }

  return cluster_members_.size();
}
