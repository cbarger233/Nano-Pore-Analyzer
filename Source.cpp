#include "Contcar.h"
#include "Atom.h"
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>

int main() {
	Contcar C;
	if (C.is_crystal()) {
		C.read_contcar_crystal();
		C.write_crystal_SiC();
		C.get_crystal();
		C.write_header();
		C.get_name();
		C.get_atomic_percents();
		C.get_density();
		C.get_parameters();
		C.get_band_gap();
		C.get_free_energy();
		C.get_cluster_sizes();
		C.write_data();
	}
	else if (C.check_files()) {
	        C.get_name();
		C.get_crystal();
		C.write_header();
		C.read_contcar();
		C.calc_volume();
		C.get_bond_densities();
		C.get_atomic_percents();
		C.get_density();
		C.get_cluster_sizes();
		C.get_parameters();
		C.get_band_gap();
		C.get_free_energy();
		C.get_bond_lengths();
		C.write_data();
		C.write_carbon_cluster();
		C.write_silicon_cluster();
	}
	
	else {
	  std::cout << "Failed to open file. Files necessary to run are CONTCAR, EIGENVAL, OSZICAR, NAME.txt, and CRYSTAL.txt." << std::endl;
	}

	return 0;
}
