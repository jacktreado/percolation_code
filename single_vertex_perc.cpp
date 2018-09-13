// Single voronoi vertex percolation calc

#include "clustertree.h"
#include "voidcluster.h"
#include "voro++.hh"
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>

const double PI = 3.1415926;

int main(int argc, char *argv[]){
	// read in data
	string N_str = argv[1];			// number of particles
	string seed_str = argv[2];		// seed
	string statsf = argv[3];		// file to save stats
	string sf = argv[4];			// file to save cluster list

	// ofstream objects
	cout << "opening files to write...";
	ofstream statsobj(statsf.c_str());
	ofstream sobj(sf.c_str());

	if (!sobj.is_open() || !statsobj.is_open()){
		cout << endl << "output file could not open, ending program..." << endl;
		cout << "stat file: " << statsf << endl;
		cout << "cluster size list file: " << sf << endl;
		return 0;		
	}

	// get numerical values for input variables
	int N,seed;

	stringstream Nss(N_str);
	Nss >> N;

	stringstream s1ss(seed_str);
	s1ss >> seed;

	double poroL,poroH,aH,aL,epsilon;

	// set bounds for critical porosity
	poroH = 0.5;
	poroL = 0.0001;	

	aL = pow(abs(3*log(poroH)/(4*PI*N)),0.333333333333333);	
	aH = pow(abs(3*log(poroL)/(4*PI*N)),0.333333333333333);

	// set tolerance
	epsilon = 1e-8;

	// run voro vertex percolation
	voroperc vorobj(N);
	vorobj.rand_sphere_pos(seed);
	vorobj.get_voro();
	vorobj.voro_edge_perc(epsilon,seed,aH,aL);

	// print results to files
	vorobj.fprint_stats(statsobj,seed,1);
	vorobj.fprint_s(sobj);

	statsobj.close();
	sobj.close();

	// finished!
	return 0;
}