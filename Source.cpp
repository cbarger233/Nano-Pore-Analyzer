#include "Contcar.h"
#include "Atom.h"
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>

int main() {
	Contcar C;
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

	return 0;
}