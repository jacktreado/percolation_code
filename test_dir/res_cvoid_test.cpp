// Find percolation point of single rcp packing of residues

#include "voidcluster.h"
#include "clustertree.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream>

#define NDIM 3

using namespace std;

int main(){
	// read in inputs
	string statstr = "stat.test";
	string sstr = "s.test";
	string pstr = "/Users/JackTreado/Jamming/ProteinVoids/cluster/res/rcp/config/res_rcp_config_N16_seed1.dat";

	// get values of input variables
	int N,NNN;
	double g;
	N = 16;
	NNN = 26;
	g = 0.3;

	// ofstream objects
	cout << "opening files to write...";
	ofstream statsobj(statstr.c_str());
	ofstream sobj(sstr.c_str());

	if (!sobj.is_open() || !statsobj.is_open()){
		cout << endl << "output file could not open, ending program..." << endl;
		cout << "stat file: " << statstr << endl;
		cout << "cluster size list file: " << sstr << endl;
		return 0;		
	}

	// first, instantiate test object to get box length in angstroms
	double B;
	int NA;
	int NS = 10;
	int NC = 3;
	
	voidcluster test_obj(pstr,NS,NDIM,NNN,NC);
	B = test_obj.get_B();
	NA = test_obj.get_NP();
	cout << "B = " << B << endl;
	cout << "NA =  = " << NA << endl;

	// get proper number of sites per side given grid spacing, instantiate second, "true" object
	NS = (int)floor(B/g);

	// also get proper number of cells
	if (NA < 500)
		NC = 3;
	else if (NA < 2000)
		NC = 4;
	else
		NC = 5;
	

	// setup percolation finder
	double aH,aL,aC,epsilon;
	aH = 5;
	aL = 0.001;
	aC = 0.5;
	epsilon = 1e-8;

	// print grid, initial info
	cout << "Beginning particle clustering algorithm for random, overlapping spheres on following grid: " << endl;
	cout << "-- N = " << N << endl;
	cout << "-- NDIM = " << NDIM << endl;
	cout << "-- NS = " << NS << endl;
	cout << "-- NC = " << NC << endl;
	cout << "-- g = " << g << endl;
	cout << "-- aH = " << aH << endl;
	cout << "-- aL = " << aL << endl;
	cout << "-- aC = " << aC << endl;

	// instantiate true object
	voidcluster respack(pstr,NS,NDIM,NNN,NC);

	// get NN list
	if (NDIM == 3){
		if (NNN==6)
			respack.cubic_lattice_F();
		else if (NNN==18)
			respack.cubic_lattice_FE();
		else if (NNN==26)
			respack.cubic_lattice_FEV();
		else{
			cout << "NNN = " << NNN << " not supported, setting NNN -> 6 and running normally." << endl;
			NNN = 6;
			respack.cubic_lattice_F();
		}
	}
	else if (NDIM == 2){
		if (NNN==4)
			respack.square_lattice_E();
		else if (NNN==8)
			respack.square_lattice_EV();
		else{
			cout << "NNN = " << NNN << " not supported, setting NNN -> 4 and running normally." << endl;
			NNN = 4;
			respack.square_lattice_E();
		}
	}

	// setup cell list
	cout << "getting cell positions..." << endl;
	respack.get_cell_positions();
	cout << "getting neighbor list for cells..." << endl;
	respack.get_cell_neighbors();
	cout << "making labels..." << endl;
	respack.label_cells();

	// get percolation
	cout << "get probe particle percolation..." << endl;
	respack.find_probe_perc(aH,aL,aC,epsilon);

	// save to file
	respack.fprint_stats(statsobj,1);
	respack.fprint_s(sobj);

	statsobj.close();
	sobj.close();
	return 0;
}