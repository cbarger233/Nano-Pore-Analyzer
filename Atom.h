#include<vector>
#include<cmath>
#include<string>
#ifndef Atom_H
#define Atom_H

class Atom {
public:
	//made these functions public so we could access them from the Contcar class
	//set and get functions
	double get_xcartesian() { return x_cartesian; }
	double get_ycartesian() { return y_cartesian; }
	double get_zcartesian() { return z_cartesian; }
	void set_xcartesian(double x) { x_cartesian = x; }
	void set_ycartesian(double y) { y_cartesian = y; }
	void set_zcartesian(double z) { z_cartesian = z; }
	void set_cartesian(double a1[3], double a2[3], double a3[3]);
	double get_xdirect() { return x_direct; }
	double get_ydirect() { return y_direct; }
	double get_zdirect() { return z_direct; }
	void set_xdirect(double x) { x_direct = x; }
	void set_ydirect(double y) { y_direct = y; }
	void set_zdirect(double z) { z_direct = z; }
	int get_color() { return color; }
	void set_color(int x) { color = x; }
	void add_partner(Atom &a);
	int get_number() {return number;}
	void set_number(int x) {number = x;}
	void tally_bond() {bonds++;}

	//vector of Atoms that are the neighbors of a certain Atom
	std::vector<Atom*> partners;
	std::vector<double> bond_lengths;

	//color for the dfs cluster size algorithm
	int color = 1;
	int number = 0;
	int bonds = 0;
	std::string type;

	//coordinates for the atom
	double x_cartesian, x_direct;
	double y_cartesian, y_direct;
	double z_cartesian, z_direct;
};


void Atom::add_partner(Atom &a) {
	partners.push_back(&a);
}
//depth-first search recursive algorithm for quantifying silicon clusters
//algorithm works based off of different "colors"
//as described in https://www.bowdoin.edu/~ltoma/teaching/cs231/spring16/Lectures/10-bfsdfs/bfsdfs.pdf
int DFS(Atom &a) {
	int count = 0;
	a.color = 2;
	for (unsigned i = 0; i < a.partners.size(); i++) {
		if (a.partners[i]->color == 1) {
			count += DFS(*a.partners[i]);
		}
	}
	a.color = 3;
	count++;
	return count;
}
void Atom::set_cartesian(double a1[3], double a2[3], double a3[3]) {
	x_cartesian = (x_direct*(a1[0])) + (y_direct*a2[0]) + (z_direct*a3[0]);
	y_cartesian = (x_direct*(a1[1])) + (y_direct*a2[1]) + (z_direct*a3[1]);
	z_cartesian = (x_direct*(a1[2])) + (y_direct*a2[2]) + (z_direct*a3[2]);
}

#endif
