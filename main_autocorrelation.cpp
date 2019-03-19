#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "clock.hpp"
#include "metropolis.hpp"
#include "time_series.hpp"
#include "wolff.hpp"

int main(int argc, char **argv) {

  const int n_meas = 10000;
  const int thermalization_sweeps = 100;

  std::array<std::ofstream, 2> files{std::ofstream("timing_metropolis.txt"),
                                     std::ofstream("timing_wolff.txt")};
  for (auto &file : files)
    file << "# L\tbeta\ttime[s/sweep]\tautocorr\t MC speed" << std::endl;

  auto make_table = [&](const std::vector<int> &Ls,
                        const std::vector<Real> &Ts) {
    for (int L : Ls) {
      std::cout << "Size: " << L << std::endl;

      // Build the two MC solvers.
      std::array<std::unique_ptr<IsingLattice>, 2> solvers{
          std::make_unique<Metropolis>(L), std::make_unique<Wolff>(L)};

      for (int alg_id = 0; alg_id < solvers.size(); ++alg_id) {
        auto &solver = solvers[alg_id];
        auto &file = files[alg_id];

        for (Real T : Ts) {
          const Real beta = 1. / T;
          std::cout << "Inverse temperature: " << beta << std::endl;

          std::vector<Real> energies(n_meas);

          solver->setBeta(beta);

          // Thermalize.
          for (int i = 0; i < thermalization_sweeps; ++i) {
            solver->doSweep();
          }

          solver->markThermalized();

          // Timed measure.
          Clock::start();
          for (int i = 0; i < n_meas; ++i) {
            solver->doSweep();
            energies[i] = solver->getE();
          }

          const Real time_sweep = Clock::stop() / Real(n_meas);

          // Compute autocorrelation function.
          auto rho = autocorrelation(energies, 100);
          // Integrated autocorrelation time.
          Real tau = integratedAutocorrelation(rho);

          // Store the results.
          file << L << "\t" << beta << "\t" << time_sweep << "\t" << tau << "\t"
               << 1. / (time_sweep * tau) << std::endl;
        }
      }
    }
  };

  constexpr Real Tc = 4.51;
  make_table(std::vector<int>{12}, std::vector<Real>{Tc - 1., Tc, Tc + 1.});
  make_table(std::vector<int>{16, 20}, std::vector<Real>{Tc});

  return 0;
}
