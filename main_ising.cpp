#include <fstream>
#include <iostream>
#include <string>

#include "ising_lattice.hpp"

int main() {
  const std::vector<Real> betas{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
  const int thermalization_sweep = 500;
  const int measurements = 10000;
  const int L = 20;
  const std::string outname = "output.txt";

  IsingLattice lattice(L, betas[0]);
  std::ofstream out(outname);
  out << "# beta\tE\tM\tE2\tM2\n";

  for (auto beta : betas) {
    std::cout << "Inverse temperature: " << beta << std::endl;
    lattice.setBeta(beta);

    // Thermalize.
    for (int i = 0; i < thermalization_sweep; ++i)
      lattice.doSweep();

    // Measure.
    Real E(0), E2(0), M(0), M2(0);
    for (int i = 0; i < measurements; ++i) {
      lattice.doSweep();
      E += lattice.getE();
      E2 += lattice.getE() * lattice.getE();
      M += lattice.getM();
      M2 += lattice.getM() * lattice.getM();
    }

    E /= lattice.size() * measurements; // Store energy density.
    E2 /= lattice.size() * lattice.size() * measurements;
    M /= lattice.size() * measurements;
    M2 /= lattice.size() * lattice.size() * measurements;
    out << beta << "\t" << E << "\t" << M << "\t" << E2 << "\t" << M2 << std::endl;
  }

  return 0;
}
