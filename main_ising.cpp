#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

#include "json11.hpp"
#include "metropolis.hpp"
#include "wolff.hpp"

int main(int argc, char **argv) {
  const std::string inputname = argc > 1 ? argv[1] : "input.json";
  std::cout << "Reading " << inputname << "\n";
  json11::Json reader = json11::Json::parseFile(inputname);

  const int thermalization_sweep = reader["thermalization-sweep"].int_value();
  const int measurements = reader["measurements"].int_value();
  const int L = reader["L"].int_value();

  std::vector<int> Ls;
  for (const auto &l : reader["Ls"].array_items()) {
    Ls.push_back(l.int_value());
  }

  std::vector<Real> betas;
  for (const auto &T : reader["Ts"].array_items()) {
    betas.push_back(1. / T.number_value());
  }
  std::sort(betas.begin(), betas.end());

  for (int L : Ls) {
    Wolff lattice(L);
    std::cout << "Size: " << L << std::endl;

    std::ofstream e_file("outputs/energies_L" + std::to_string(L) + ".txt");
    e_file << "# beta\tE\n";
    e_file.precision(10);

    std::ofstream m_file("outputs/magnetization_L" + std::to_string(L) +
                         ".txt");
    m_file << "# beta\tM\n";
    m_file.precision(10);

    for (Real beta : betas) {
      std::cout << "Inverse temperature: " << beta << std::endl;
      lattice.setBeta(beta);

      e_file << beta << "\t";
      m_file << beta << "\t";

      // Thermalize.
      for (int i = 0; i < thermalization_sweep; ++i) {
        lattice.doSweep();
      }

      lattice.markThermalized();

      // Measure.
      for (int i = 0; i < measurements; ++i) {
        lattice.doSweep();
        e_file << lattice.getE() << "\t";
        m_file << lattice.getM() << "\t";
      }
      e_file << "\n";
      m_file << "\n";
    }
  }

  return 0;
}
