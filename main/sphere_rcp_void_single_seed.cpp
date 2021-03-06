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
#define sph 1

using namespace std;

int main(int argc, char *argv[]){
	// read in inputs
	string Nstr = argv[1];
	string NNNstr = argv[2];
	string gstr = argv[3];
	string pstr = argv[4];
	string statstr = argv[5];
	string sstr = argv[6];

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

	// get values of input variables
	int N,NNN;
	double g;
	stringstream Nss(Nstr);
	Nss >> N;

	stringstream NNNss(NNNstr);
	NNNss >> NNN;

	stringstream gss(gstr);
	gss >> g;

	// first, instantiate test object to get box length in angstroms
	double B;
	int NS = 10;
	int NC = 3;
	voidcluster test_obj(pstr,sph,NS,NDIM,NNN,NC);
	B = test_obj.get_B();
	cout << "B = " << B << endl;

	// get proper number of sites per side given grid spacing, instantiate second, "true" object
	NS = (int)floor(B/g);

	// also get proper number of cells
	int ncmin, ncmax;
	double rmax = 2;
	NC = (int)floor(B/rmax);
	ncmin = 4;
	ncmax = (int)round(pow(13*N,1/((double)NDIM)))-1;
	if (NC > ncmax)
		NC = ncmax;
	else if (NC < ncmin){
		NC = ncmin;
	}

	// setup percolation finder
	double aH,aL,aC,epsilon;
	aH = 1.3;
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
	voidcluster respack(pstr,sph,NS,NDIM,NNN,NC);

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