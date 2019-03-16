#pragma once

#include <array>
#include <queue>

#include "ising_lattice.hpp"
#include "typedefs.hpp"

class Wolff final : public IsingLattice {
public:
  Wolff(int L);

  void doSweep() override;

  void setBeta(Real beta) override;

  void markThermalized() override;

private:
  // Returns the size of the modified cluster.
  std::size_t doStep();

  auto neighbours(const std::array<unsigned int, 3>& site);
  std::array<unsigned int, 3> unrollIndex(unsigned int idx);

  Real prob_;
  std::size_t cumulative_cluster_size_ = 0;
  std::size_t n_steps_ = 0;
  int steps_per_sweep_ = 0;

  std::queue<std::array<unsigned int, 3>> cluster_queue_;
};
