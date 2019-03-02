#include <fstream>
#include <iostream>
#include <string>

#include "ising_lattice.hpp"
#include "json11.hpp"

int main(int argc, char **argv) {
  const std::string inputname = argc > 1 ? argv[1] : "input.json";
  std::cout << "Reading " << inputname << "\n";
  json11::Json reader = json11::Json::parseFile(inputname);

  const Real start_beta = reader["initial-beta"].number_value();
  const Real end_beta = reader["final-beta"].number_value();
  const int n_beta = reader["n-betas"].int_value();
  const int thermalization_sweep = reader["thermalization-sweep"].int_value();
  const int measurements = reader["measurements"].int_value();
  const int L = reader["L"].int_value();

  const std::string outname = "output.txt";
  const Real delta_beta = (end_beta - start_beta) / n_beta;

  IsingLattice lattice(L, start_beta);
  std::ofstream out(outname);
  out << "# beta\tE\tM\tE2\tM2\n";

  for (Real beta = start_beta; beta < end_beta; beta += delta_beta) {
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
    out << beta << "\t" << E << "\t" << M << "\t" << E2 << "\t" << M2
        << std::endl;
  }

  return 0;
}
