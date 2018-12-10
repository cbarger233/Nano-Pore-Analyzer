#include"Atom.h"
#include<iostream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#ifndef Contcar_H
#define Contcar_H

//constant for calculating bond length
const double SiC_bond_length = 2.24;	//maximum bond length between silicon and carbon atoms
const double CH_bond_length = 1.2;	//maximum bond length between carbon and hydrogen atoms
const double SiH_bond_length = 1.6;	//etc
const double CC_bond_length = 1.9;
const double SiSi_bond_length = 2.5;
const double HH_bond_length = 1.0;
const double mass_H = 1.00797;		//mass of a hydrogen atom in amu, used in density calculation
const double mass_C = 12.011;
const double mass_Si = 28.0855;



class Contcar {
public:

//protected:
	double SiC_crystalline = 0;		//SiC bond density for the crystalline structure to be read in by get_crystal();
	std::string system_name;		
	double lattice_constant = 0;		//self-explanatory
	double a1[3], a2[3], a3[3];		//basis vectors for the system
	std::string coordinate_type;		//direct or cartesian
	std::string types[3];			//array to store atom types
	int amounts[3];				//parallel array to store amounts of each atom
	void read_contcar();			//function that allows you to read the VASP CONTCAR file
	
	//atoms that are in the CONTCAR file
	std::vector<Atom> carbons;		//the above read_contcar() function reads in the atoms and stores them in these vectors
	std::vector<Atom> silicons;
	std::vector<Atom> hydrogens;
	
	double volume = 0;	//volume of the cell in cm^3
	void calc_volume();	//function to calculate the volume
	
	double SiC = 0, CH = 0, SiH = 0, CC = 0, SiSi = 0, HH = 0;	//number of bonds
	double percentSi = 0, percentC = 0, percentH = 0;		//atomic percentages
	double density = 0;						//density in g/cm^3

	//bond densities.. for example, SiC is the bond density of Si-C bonds
	//in bonds per cm^3
	double SiC_p = 0, CH_p = 0, SiH_p = 0, CC_p = 0, SiSi_p = 0, HH_p = 0;

	//"n" variables used to describe the data
	double n_SiC = 0, n_CH = 0, n_SiH = 0, n_CC = 0, n_SiSi = 0;
	double band_gap;

	//free energy of the system
	double free_e = 0;
	
	//boolean variable to see if there are coordination defects
	bool coordination_defect = false;

	//vectors of cluster sizes
	std::vector<int> silicon_cluster_sizes;
	std::vector<int> carbon_cluster_sizes;
	int silicon_clusters[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	int carbon_clusters[12] = {0,0,0,0,0,0,0,0,0,0,0,0};

	//statistical info
	double mean_SiSi = 0;	//mean Si-Si bond length
	double mean_CC = 0;
	double mean_SiC = 0;
	double MAD_SiSi = 0;	//mean average deviation of the SiSi bonds
	double MAD_CC = 0;
	double MAD_SiC = 0;
	double STD_SiSi = 0;	//standard deviation of the SiSi bonds
	double STD_CC = 0;
	double STD_SiC = 0;

	//stuff for bond length distributions. bond lengths are calculated and stored in these vectors
	std::vector<double> sisi_bondlengths;
	std::vector<double> cc_bondlengths;
	std::vector<double> sic_bondlengths;
	
	bool is_crystal();				//checks to see if the system is the crystal
	void read_contcar_crystal();			//reads the crystal's contcar
	void write_crystal_SiC();			//writes out the crystal SiC density to a file
	void get_crystal();				//get the SiC bond density from the CRYSTAL.txt file
	void get_name();                                //gets the name of the system
	bool check_files();				//checks to see if all of the necessary files can be opened and looked at
	void write_header();				//writes the header file if one isnt already written
	double find_distance(Atom one, Atom two);	//calculates the distance between two atoms
	void get_coordination();			//function to find coordination defects
	void get_bond_densities();			//calculates the bond densities in the system
	void get_atomic_percents();			//calculates atomic percents
	void get_density();				//calculates the desnity of the system
	void get_cluster_sizes();			//calculates the cluster sizes of CC and SiSi bonds
	void get_parameters();				//calculates parameters that describe the system
	void get_band_gap();				//calculates the band gap, requires the VASP EIGENVAL file
	void write_silicon_cluster();			//writes a POSCAR-like file that shows all of the SiSi bonds
	void write_carbon_cluster();			//same as above but for CC
	void get_free_energy();				//function that calculates the free energy of the system
	void write_data();				//writes the data out to the header file
	void get_bond_lengths();			//calculates the bond lengths
};

void Contcar::get_name() {
  std::ifstream fin;
  fin.open("NAME.txt");
  getline(fin, system_name);
}

bool Contcar::is_crystal() {
	bool ok = true;
	
	std::ifstream fin;
	fin.open("IS_CRYSTAL.txt");
	if (!fin.is_open())
		ok = false;
	return ok;
}

bool Contcar::check_files() {
	bool ok = true;
	std::ifstream fin;
	fin.open("CONTCAR");
	std::ifstream in;
	in.open("OSZICAR");
	std::ifstream name;
	name.open("NAME.txt");
	std::ifstream eig;
	eig.open("EIGENVAL");
	std::ifstream cry;
	cry.open("CRYSTAL.txt");
	if (!fin.is_open() || !in.is_open() || !name.is_open() || !eig.is_open() || !cry.is_open()) {
		ok = false;
	}
	return ok;
}

void Contcar::write_header() {
	std::ifstream fin;
	fin.open("DATA.csv");
	if (!fin.is_open()) {
		std::ofstream fout;
		fout.open("DATA.csv");
		fout << "Name, Band Gap (eV), n*_SiC, n_CH, n_SiH, n_SiSi, n_CC, SiC density, CH density, SiH density, SiSi density, CC density, HH density, Mean SiSi Length, MAD SiSi, STD SiSi, Mean CC Length, MAD CC, STD CC, Mean SiC Length, MAD SiC, STD SiC, %C, %Si, %H, Density(g/cm^3), Volume(cm^3), Free Energy, Coordination_Defects << ,";
		for (unsigned i = 1; i < 12; i++) { fout << i << "C, "; };
		for (unsigned j = 1; j < 12; j++) { fout << j << "Si, "; };
	}
}


void Contcar::read_contcar_crystal() {
	std::ifstream fin;
	fin.open("CONTCAR");
	std::string d = "";
	getline(fin, d);
	fin >> lattice_constant;
	fin >> a1[0] >> a1[1] >> a1[2];
	fin >> a2[0] >> a2[1] >> a2[2];
	fin >> a3[0] >> a3[1] >> a3[2];
	fin >> types[0] >> types[1];
	fin >> amounts[0] >> amounts[1];
	fin >> coordinate_type;

	Atom temp;
	double x, y, z;

	//reading in the first type of atom
	char test = types[0][0];
	for (int i = 0; i < amounts[0]; i++) {
		fin >> x >> y >> z;
		temp.set_xdirect(x);
		temp.set_ydirect(y);
		temp.set_zdirect(z);
		temp.set_color(1);
		switch (test)
		{
		case 'C':
			carbons.push_back(temp);
			break;
		case 'S':
			silicons.push_back(temp);
			break;
		case 'H':
			hydrogens.push_back(temp);
		default:
			break;
		}
	}

	//reading in the second type of atom
	test = types[1][0];
	for (int i = 0; i < amounts[1]; i++) {
		fin >> x >> y >> z;
		temp.set_xdirect(x);
		temp.set_ydirect(y);
		temp.set_zdirect(z);
		temp.set_color(1);
		switch (test)
		{
		case 'C':
			carbons.push_back(temp);
			break;
		case 'S':
			silicons.push_back(temp);
			break;
		case 'H':
			hydrogens.push_back(temp);
		default:
			break;
		}
	}
	
	//set the cartesian coordinates for all of the atoms
	for (unsigned i = 0; i < carbons.size(); i++) { carbons[i].set_cartesian(a1, a2, a3); }
	for (unsigned j = 0; j < silicons.size(); j++) { silicons[j].set_cartesian(a1, a2, a3); }
	
	volume = a1[0] * (a2[1] * a3[2] - a2[2] * a3[1]) - a1[1] * (a2[0] * a3[2] - a2[2] * a3[0]) + a1[2] * (a2[0] * a3[1] - a2[1] * a3[0]);
	volume = volume*pow(10.0, -24.0);
	
	double distance = 0;
	
	//setting SiC bond density
	for (unsigned i = 0; i < silicons.size(); i++) {
		for (unsigned j = 0; j < carbons.size(); j++) {
			distance = find_distance(silicons[i], carbons[j]);
			if (distance < SiC_bond_length && distance > 0.1) {
				//silicons[i].partners.push_back(&carbons[j]);
				//carbons[j].partners.push_back(&silicons[i]);
				SiC = SiC + 1.0;
			}
		}
	}
	SiC_p = SiC * pow(10, -21) / volume;
	
}


void Contcar::write_crystal_SiC() {
	std::ofstream fout;
	fout.open("CRYSTAL.txt");
	fout << SiC_p;
}


void Contcar::get_crystal() {
	std::ifstream fin;
	fin.open("CRYSTAL.txt");
	fin >> SiC_crystalline;
}


void Contcar::read_contcar() {
	std::ifstream fin;
	fin.open("CONTCAR");
	std::string d = "";
	getline(fin, d);
	fin >> lattice_constant;
	fin >> a1[0] >> a1[1] >> a1[2];
	fin >> a2[0] >> a2[1] >> a2[2];
	fin >> a3[0] >> a3[1] >> a3[2];
	fin >> types[0] >> types[1] >> types[2];
	fin >> amounts[0] >> amounts[1] >> amounts[2];
	fin >> coordinate_type;

	Atom temp;
	double x, y, z;

	//reading in the first type of atom
	char test = types[0][0];
	for (int i = 0; i < amounts[0]; i++) {
		fin >> x >> y >> z;
		temp.set_xdirect(x);
		temp.set_ydirect(y);
		temp.set_zdirect(z);
		temp.set_color(1);
		switch (test)
		{
		case 'C':
			carbons.push_back(temp);
			break;
		case 'S':
			silicons.push_back(temp);
			break;
		case 'H':
			hydrogens.push_back(temp);
		default:
			break;
		}
	}

	//reading in the second type of atom
	test = types[1][0];
	for (int i = 0; i < amounts[1]; i++) {
		fin >> x >> y >> z;
		temp.set_xdirect(x);
		temp.set_ydirect(y);
		temp.set_zdirect(z);
		temp.set_color(1);
		switch (test)
		{
		case 'C':
			carbons.push_back(temp);
			break;
		case 'S':
			silicons.push_back(temp);
			break;
		case 'H':
			hydrogens.push_back(temp);
		default:
			break;
		}
	}

	//reading in the third atom type
	test = types[2][0];
	for (int i = 0; i < amounts[2]; i++) {
		fin >> x >> y >> z;
		temp.set_xdirect(x);
		temp.set_ydirect(y);
		temp.set_zdirect(z);
		temp.set_color(1);
		switch (test)
		{
		case 'C':
			carbons.push_back(temp);
			break;
		case 'S':
			silicons.push_back(temp);
			break;
		case 'H':
			hydrogens.push_back(temp);
		default:
			break;
		}
	}

	//set the cartesian coordinates for all of the atoms
	for (unsigned i = 0; i < carbons.size(); i++) { carbons[i].set_cartesian(a1, a2, a3); }
	for (unsigned j = 0; j < silicons.size(); j++) { silicons[j].set_cartesian(a1, a2, a3); }
	for (unsigned k = 0; k < hydrogens.size(); k++) { hydrogens[k].set_cartesian(a1, a2, a3); }
}

//function that calculates the volume. Uses the lattice vectors
//and calculates the volume of the parallelepiped formed by them
void Contcar::calc_volume() {
	volume = a1[0] * (a2[1] * a3[2] - a2[2] * a3[1]) - a1[1] * (a2[0] * a3[2] - a2[2] * a3[0]) + a1[2] * (a2[0] * a3[1] - a2[1] * a3[0]);
	volume = volume*pow(10.0, -24.0);
}

//function to find the minimum distance between two atoms
//takes into account the periodic boundary conditions
double Contcar::find_distance(Atom one, Atom two) {
	long double distance, d;
	long double x_periodic, y_periodic, z_periodic;

	long double tempx = one.get_xcartesian() - two.get_xcartesian();
	long double tempy = one.get_ycartesian() - two.get_ycartesian();
	long double tempz = one.get_zcartesian() - two.get_zcartesian();

	distance = sqrt(pow(tempx, 2.0) + pow(tempy, 2.0) + pow(tempz, 2.0));
	
	//calculating the distance thus far has been pretty straightforward
	//but now we calculate the length between all of the atoms
	//considering the periodic boundary conditions and takes the minimum
	//of the distance calculated above and the one calculated below

	for (int f = -1; f < 2; f++) {
		for (int g = -1; g < 2; g++) {
			for (int h = -1; h < 2; h++) {
				x_periodic = tempx + a1[0] * f + a2[0] * g + a3[0] * h;
				y_periodic = tempy + a1[1] * f + a2[1] * g + a3[1] * h;
				z_periodic = tempz + a1[2] * f + a2[2] * g + a3[2] * h;
				d = sqrt(pow(x_periodic, 2.0) + pow(y_periodic, 2.0) + pow(z_periodic, 2.0));
				if (d < distance)
					distance = d;
			}
		}
	}
	return distance;
}

void Contcar::get_coordination() {
	int bonds = 0;
	double distance = 0;
	
	for (unsigned i = 0; i < carbons.size(); i++) {
		for (unsigned j = 0; j < silicons.size(); j++) {
			for (unsigned k = 0; k < hydrogens.size(); k++) {
				distance = find_distance(silicons[j], carbons[i]);
				if (distance < SiC_bond_length && distance > 0.1)
					bonds++;
				distance = find_distance(carbons[i], hydrogens[k]);
				if (distance < CH_bond_length && distance > 0.1)
					bonds++;
			}
		}
		if (bonds < 4)
			coordination_defect = true;
		bonds = 0;
	}
	
	
	for (unsigned i = 0; i < silicons.size(); i++) {
		for (unsigned j = 0; j < carbons.size(); j++) {
			for (unsigned k = 0; k < hydrogens.size(); k++) {
				distance = find_distance(silicons[i], carbons[j]);
				if (distance < SiC_bond_length && distance > 0.1)
					bonds++;
				distance = find_distance(silicons[i], hydrogens[k]);
				if (distance < SiH_bond_length && distance > 0.1)
					bonds++;
			}
		}
		if (bonds < 4)
			coordination_defect = true;
		bonds = 0;
	}
}
//function to get the bond densities
//we use nested for-loops that loop over each type of atom
//and calculates the distance between them in the process
//the distance is compared to the maximum bond lengths
//and then bonds are assigned to the ones that meet the criteria
void Contcar::get_bond_densities() {

	double distance = 0;

	//setting SiC bond density
	for (unsigned i = 0; i < silicons.size(); i++) {
		for (unsigned j = 0; j < carbons.size(); j++) {
			distance = find_distance(silicons[i], carbons[j]);
			if (distance < SiC_bond_length && distance > 0.1) {
				//silicons[i].partners.push_back(&carbons[j]);
				//carbons[j].partners.push_back(&silicons[i]);
				sic_bondlengths.push_back(distance);
				SiC = SiC + 1.0;
			}
		}
	}
	SiC_p = SiC * pow(10, -21) / volume;

	//setting CH bond density
	for (unsigned i = 0; i < carbons.size(); i++) {
		for (unsigned j = 0; j < hydrogens.size(); j++) {
			distance = find_distance(carbons[i], hydrogens[j]);
			if (distance < CH_bond_length && distance > 0.1) {
				//carbons[i].add_partner(hydrogens[j]);
				//hydrogens[j].add_partner(carbons[i]);
				CH += 1.0;
			}
		}
	}
	CH_p = CH * pow(10, -21) / volume;

	//setting SiH bond density
	for (unsigned i = 0; i < silicons.size(); i++) {
		for (unsigned j = 0; j < hydrogens.size(); j++) {
			distance = find_distance(silicons[i], hydrogens[j]);
			if (distance < SiH_bond_length && distance > 0.1) {
				//silicons[i].add_partner(hydrogens[j]);
				//hydrogens[j].add_partner(silicons[i]);
				SiH += 1.0;
			}
		}
	}
	SiH_p = SiH * pow(10, -21) / volume;

	//setting CC bond density
	for (unsigned i = 0; i < carbons.size(); i++) {
		for (unsigned j = 0; j < carbons.size(); j++) {
			distance = find_distance(carbons[i], carbons[j]);
			if (distance < CC_bond_length && distance > 0.1) {
				carbons[i].partners.push_back(&carbons[j]);
				carbons[i].bond_lengths.push_back(distance);
				carbons[j].partners.push_back(&carbons[i]);
				carbons[j].bond_lengths.push_back(0);
				cc_bondlengths.push_back(distance);
				CC += 1.0;
			}
		}
	}
	CC_p = CC * pow(10, -21) / volume;

	//setting SiSi bond density
	for (unsigned i = 0; i < silicons.size(); i++) {
		for (unsigned j = 0; j < silicons.size(); j++) {
			distance = find_distance(silicons[i], silicons[j]);
			if (distance < SiSi_bond_length && distance > 0.1) {
				silicons[i].partners.push_back(&silicons[j]);
				silicons[i].bond_lengths.push_back(distance);
				silicons[j].partners.push_back(&silicons[i]);
				silicons[j].bond_lengths.push_back(0);			//push back 0 so that when we loop over all bond lengths for atoms we don't double count
				sisi_bondlengths.push_back(distance);
				SiSi += 1.0;
			}
		}
	}
	SiSi_p = SiSi * pow(10, -21) / volume;

	//setting HH bond density
	for (unsigned i = 0; i < hydrogens.size(); i++) {
		for (unsigned j = 0; j < hydrogens.size(); j++) {
			distance = find_distance(hydrogens[i], hydrogens[j]);
			if (distance < HH_bond_length && distance > 0.1) {
				//hydrogens[i].add_partner(hydrogens[j]);
				//hydrogens[j].add_partner(hydrogens[i]);
				HH += 1.0;
			}
		}
	}
	HH_p = HH * pow(10, -21) / volume;
}

//pretty straightforward function that calculates the atomic percentages
void Contcar::get_atomic_percents() {
	int total_atoms = silicons.size() + carbons.size() + hydrogens.size();
	percentC = static_cast<double>(carbons.size()) / static_cast<double>(total_atoms);
	percentH = static_cast<double>(hydrogens.size()) / static_cast<double>(total_atoms);
	percentSi = static_cast<double>(silicons.size()) / static_cast<double>(total_atoms);
}

//fucntion that calculates the desnity
void Contcar::get_density() {
	double mass = 0;
	mass = mass + (silicons.size()*mass_Si) + (carbons.size()*mass_C) + (hydrogens.size()*mass_H);
	mass = mass*(1.66054 * pow(10.0, -24.0));
	density = mass / volume;
}

//function that calculates
void Contcar::get_cluster_sizes() {

	//finding clusters of silicon atoms
	for (unsigned i = 0; i < silicons.size(); i++) {
		if (silicons[i].get_color() == 1) {
			silicon_cluster_sizes.push_back(DFS(silicons[i]));
		}
	}

	for (int i = 0; i < 12; i++) {
		silicon_clusters[i] = 0;
	}

	for (unsigned i = 0; i < silicon_cluster_sizes.size(); i++ ) {
		int a = silicon_cluster_sizes[i];
		if (a > 11)
			a = 11;
		silicon_clusters[a]++;
	}

	//finding clusters of carbon atoms
	for (unsigned i = 0; i < carbons.size(); i++) {
		if (carbons[i].get_color() == 1) {
			carbon_cluster_sizes.push_back(DFS(carbons[i]));
		}
	}

	for (int i = 0; i < 12; i++) {
		carbon_clusters[i] = 0;
	}

	for (unsigned i = 0; i < carbon_cluster_sizes.size(); i++) {
		int b = carbon_cluster_sizes[i];
		if (b > 11)
			b = 11;
		carbon_clusters[b]++;
	}


}
void Contcar::get_parameters() {
	n_SiC = 1 - (SiC_p / SiC_crystalline);
	n_CH = CH_p / SiC_crystalline;
	n_SiH = SiH_p / SiC_crystalline;
	n_SiSi = SiSi_p / SiC_crystalline;
	n_CC = CC_p / SiC_crystalline;
}
void Contcar::get_band_gap() {
	std::ifstream fin;
	fin.open("EIGENVAL");
	std::string dummy;

	//ignores the first ten lines or so of the file because we dont really need those to calculate the band gap
	for (unsigned i = 0; i < 13; i++) {
		getline(fin, dummy);
	}

	double num = 0, energy1 = 0, energy2 = 0, occupation = -1;
	band_gap = 0;

	while (!fin.eof()) {
		energy1 = energy2;
		fin >> num >> energy2 >> occupation;
		if (occupation == 0.00000) {
			band_gap = energy2 - energy1;
			break;
		}
	}
}
void Contcar::write_silicon_cluster() {

	std::vector<Atom> temp;
	for (unsigned i = 0; i < silicons.size(); i++) {
		if (!silicons[i].partners.empty()) {
			temp.push_back(silicons[i]);
		}
		std::ofstream fout;
		fout.open("SILICON.vasp");
		fout << system_name << std::endl;
		fout << lattice_constant << std::endl;
		fout << a1[0] << " " << a1[1] << " " << a1[2] << std::endl;
		fout << a2[0] << " " << a2[1] << " " << a2[2] << std::endl;
		fout << a3[0] << " " << a3[1] << " " << a3[2] << std::endl;
		fout << "Si" << std::endl;
		fout << temp.size() << std::endl;
		fout << "Direct" << std::endl;
		fout << std::fixed;
		for (unsigned i = 0; i < temp.size(); i++) {
			fout << std::setprecision(8) << temp[i].x_direct << " " << temp[i].y_direct << " " << temp[i].z_direct << std::endl;
		}
		fout.close();
	}
}
void Contcar::write_carbon_cluster() {

	std::vector<Atom> temp;
	for (unsigned i = 0; i < carbons.size(); i++) {
		if (!carbons[i].partners.empty()) {
			temp.push_back(carbons[i]);
		}
		std::ofstream fout;
		fout.open("CARBONS.vasp");
		fout << system_name << std::endl;
		fout << lattice_constant << std::endl;
		fout << a1[0] << " " << a1[1] << " " << a1[2] << std::endl;
		fout << a2[0] << " " << a2[1] << " " << a2[2] << std::endl;
		fout << a3[0] << " " << a3[1] << " " << a3[2] << std::endl;
		fout << "C" << std::endl;
		fout << temp.size() << std::endl;
		fout << "Direct" << std::endl;
		fout << std::fixed;
		for (unsigned i = 0; i < temp.size(); i++) {
			fout << std::setprecision(8) << temp[i].x_direct << " " << temp[i].y_direct << " " << temp[i].z_direct << std::endl;
		}
		fout.close();
	}
}
void Contcar::get_free_energy() {
	std::ifstream fin;
	fin.open("OSZICAR");
	std::string var;

	while (fin >> var) {
		if (var == "E0=") {
			fin >> free_e;
		}
	}
}
void Contcar::write_data() {
	std::ofstream os;
	os.open("DATA.csv", std::ios_base::app);
	os << 0 << std::endl;

	//wrting out the actual numerical values
	os << system_name << ", " << band_gap << ", " << n_SiC << ", " << n_CH << ", " << n_SiH << ", " << n_SiSi << ", " << n_CC << ", " << SiC_p << ", " << CH_p << ", " << SiH_p << ", " << SiSi_p << ", " << CC_p << ", " << HH_p << ", " << mean_SiSi << ", " << MAD_SiSi << ", "  << STD_SiSi << ", " << mean_CC << ", " << MAD_CC << ", " << STD_CC << ", " << mean_SiC << ", " << MAD_SiC << ", " << STD_SiC << ", " << percentC << ", " << percentSi << ", " << percentH << ", " << density << ", " << volume << ", " << free_e << ", " << coordination_defect << ", ";
	for (unsigned i = 1; i < 12; i++) { os << carbon_clusters[i] << ", "; };
	for (unsigned j = 1; j < 12; j++) { os << silicon_clusters[j] << ", "; };

}

void Contcar::get_bond_lengths() {
	//for the SiSi bonds
	if (sisi_bondlengths.size() > 0) {
		double total = 0;
		double low = 5, high = 0;
		for (unsigned i = 0; i < sisi_bondlengths.size(); i++) {
			total += sisi_bondlengths[i];
			if (sisi_bondlengths[i] > high)
				high = sisi_bondlengths[i];
			if (sisi_bondlengths[i] < low)
				low = sisi_bondlengths[i];
		}
		mean_SiSi = total / sisi_bondlengths.size();

		double sum = 0;
		double num = 0;
		for (unsigned i = 0; i < sisi_bondlengths.size(); i++) {
		  num = sisi_bondlengths[i] - mean_SiSi;
		  if(num < 0)
		    num = num*(-1);
		  sum += num;
		}
		MAD_SiSi = sum / sisi_bondlengths.size();
		sum = 0;

		for (unsigned i = 0; i < sisi_bondlengths.size(); i++) {
			sum = sum + pow((sisi_bondlengths[i] - mean_SiSi), 2);
		}
		STD_SiSi = pow(sum / (sisi_bondlengths.size() - 1), 0.5);
	}

	/***********************************************************************************
	************************************************************************************
	************************************************************************************/
	//for the CC bonds
	if (cc_bondlengths.size() > 0) {
		double total = 0;
		double low = 5;
		double high = 0;
		for (unsigned i = 0; i < cc_bondlengths.size(); i++) {
			total += cc_bondlengths[i];
			if (cc_bondlengths[i] > high)
				high = cc_bondlengths[i];
			if (cc_bondlengths[i] < low)
				low = cc_bondlengths[i];
		}
		mean_CC = total / cc_bondlengths.size();
		double sum = 0;
		double num = 0;
		for (unsigned i = 0; i < cc_bondlengths.size(); i++) {
		  num = cc_bondlengths[i] - mean_CC;
		  if (num < 0) {
		    num = num*(-1);
		  }
		  sum += num;
		}
		MAD_CC = sum / cc_bondlengths.size();
		sum = 0;

		for (unsigned i = 0; i < cc_bondlengths.size(); i++) {
			sum += pow((cc_bondlengths[i] - mean_CC), 2);
		}
		STD_CC = pow(sum / (cc_bondlengths.size() - 1), 0.5);
	}
	
	
	/***********************************************************************************
	************************************************************************************
	************************************************************************************/
	//for the SiC bonds
	if (sic_bondlengths.size() > 0) {
		double total = 0;
		double low = 5;
		double high = 0;
		for (unsigned i = 0; i < sic_bondlengths.size(); i++) {
			total += sic_bondlengths[i];
			if (sic_bondlengths[i] > high)
				high = sic_bondlengths[i];
			if (sic_bondlengths[i] < low)
				low = sic_bondlengths[i];
		}
		mean_SiC = total / sic_bondlengths.size();
		double sum = 0;
		double num = 0;
		for (unsigned i = 0; i < sic_bondlengths.size(); i++) {
		  num = sic_bondlengths[i] - mean_SiC;
		  if (num < 0) {
		    num = num*(-1);
		  }
		  sum += num;
		}
		MAD_SiC = sum / sic_bondlengths.size();
		sum = 0;

		for (unsigned i = 0; i < sic_bondlengths.size(); i++) {
			sum += pow((sic_bondlengths[i] - mean_SiC), 2);
		}
		STD_SiC = pow(sum / (sic_bondlengths.size() - 1), 0.5);
	}
}

#endif
