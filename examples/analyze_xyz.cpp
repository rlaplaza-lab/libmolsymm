#include "molsymm/molsymm.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static int atomic_number_from_symbol(const std::string& sym) {
  static const char* symbols[] = {"",   "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",
                                  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",
                                  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
                                  "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr"};
  for (int z = 1; z <= 36; ++z) {
    if (sym == symbols[z]) {
      return z;
    }
  }
  return 0;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "usage: analyze_xyz file.xyz\n";
    return 1;
  }
  std::ifstream in(argv[1]);
  if (!in) {
    std::cerr << "cannot open " << argv[1] << "\n";
    return 1;
  }
  int nat = 0;
  std::string line;
  std::getline(in, line);
  nat = std::stoi(line);
  std::getline(in, line); // comment
  std::vector<int> z;
  std::vector<double> xyz;
  for (int i = 0; i < nat; ++i) {
    std::getline(in, line);
    std::istringstream ss(line);
    std::string sym;
    double x;
    double y;
    double zz;
    ss >> sym >> x >> y >> zz;
    z.push_back(atomic_number_from_symbol(sym));
    xyz.push_back(x);
    xyz.push_back(y);
    xyz.push_back(zz);
  }

  molsymm::Molecule mol(z, xyz);
  molsymm::Analysis a = molsymm::analyze(mol);
  std::cout << "point group: " << a.point_group << "\n";
  std::cout << "symmetry number: " << a.symmetry_number << "\n";
  std::cout << "order: " << a.order << "\n";
  std::cout << "rmsd: " << a.rmsd << "\n";
  auto disp = molsymm::unique_displacements(a);
  std::cout << "unique displacements: " << disp.displacements.size() << "\n";
  return 0;
}
