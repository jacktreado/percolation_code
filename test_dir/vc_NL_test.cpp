/*

	Get info for single seed

*/

#include "clustertree.h"
#include "voidcluster.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>

#define L 1
#define NL_DEBUG

const double PI = 3.1415269;
const double epsilon = 1e-8;

using namespace std;

int main(int argc, char *argv[]){
	// // read in data
	// string N_str = argv[1];			// number of particles
	// string NDIM_str = argv[2];		// dimension of space
	// string NNN_str = argv[3];		// number of nearest neighbors (decides lattice topology)
	// string g0_str = argv[4];		// grid spacing amplitude g0
	// string seed_str = argv[5];		// seed
	// string statsf = argv[6];		// file to save stats
	// string sf = argv[7];			// file to save cluster list

	// // ofstream objects
	// cout << "opening files to write...";
	// ofstream statsobj(statsf.c_str());
	// ofstream sobj(sf.c_str());

	// if (!sobj.is_open() || !statsobj.is_open()){
	// 	cout << endl << "output file could not open, ending program..." << endl;
	// 	cout << "stat file: " << statsf << endl;
	// 	cout << "cluster size list file: " << sf << endl;
	// 	return 0;		
	// }

	// // get numerical values for input variables
	// int N,NDIM,NNN,seed;
	// double g0;

	// stringstream Nss(N_str);
	// Nss >> N;

	// stringstream g0ss(g0_str);
	// g0ss >> g0;

	// stringstream NDIMss(NDIM_str);
	// NDIMss >> NDIM;

	// stringstream NNNss(NNN_str);
	// NNNss >> NNN;

	// stringstream s1ss(seed_str);
	// s1ss >> seed;

	int N = 1000;
	int NDIM = 3;
	int NNN = 26;
	double g0 = 0.1;
	int seed = 3;
	string statsf = "stat.dat";
	string sf = "clist.dat";
	string xyzf = "cell.xyz";

	// ofstream objects
	cout << "opening files to write...";
	ofstream statsobj(statsf.c_str());
	ofstream sobj(sf.c_str());
	ofstream xyzobj(xyzf.c_str());

	if (!sobj.is_open() || !statsobj.is_open()){
		cout << endl << "output file could not open, ending program..." << endl;
		cout << "stat file: " << statsf << endl;
		cout << "cluster size list file: " << sf << endl;
		return 0;		
	}


	// get bounds on radius
	double poroH,aH,poroL,aL,poroC,aC;
	poroH = 0.5;
	poroL = 0.005;	

	if (NDIM == 3){
		poroC = 0.03;

		aL = pow(abs(3*log(poroH)/(4*PI*N)),0.333333333333333);	
		aH = pow(abs(3*log(poroL)/(4*PI*N)),0.333333333333333);
		aC = pow(abs(3*log(poroC)/(4*PI*N)),0.333333333333333);
	}
	else if (NDIM == 2){
		poroC = 0.15;

		aL = sqrt(abs(log(poroH))/(PI*N));
		aH = sqrt(abs(log(poroL))/(PI*N));
		aC = sqrt(abs(log(poroC))/(PI*N));
	}
	else{
		cout << "NDIM = " << NDIM << " not yet supported, rerun program!" << endl;
		return 0;
	}

	// get grid size info
	double g;
	int NS;

	g = g0*aC;
	NS = round(L/g);

	// get number of cells
	// get number of cells
	int ncmin, ncmax, NC;
	double bl;
	ncmin = 4;
	ncmax = 20;
	bl = 1.5*aC;
	NC = floor(L/bl);
	if (NC < ncmin)
		NC = ncmin;
	else if (NC > ncmax)
		NC = ncmax;


	// print grid, initial info
	cout << "Beginning void clustering algorithm for random, overlapping spheres on following grid: " << endl;
	cout << "-- N = " << N << endl;
	cout << "-- NDIM = " << NDIM << endl;
	cout << "-- NS = " << NS << endl;
	cout << "-- NC = " << NC << endl;
	cout << "-- seed = " << seed << endl;
	cout << "-- g = " << g << endl;
	cout << "-- aH = " << aH << "; poroL = " << poroL << endl;
	cout << "-- aL = " << aL << "; poroH = " << poroH << endl;
	cout << "-- aC = " << aC << "; poroC = " << poroC << endl;

	// construct voicluster object
	voidcluster cp(N,NS,NDIM,NNN,NC); 

	// get NN list
	if (NDIM == 3){
		if (NNN==6)
			cp.cubic_lattice_F();
		else if (NNN==18)
			cp.cubic_lattice_FE();
		else if (NNN==26)
			cp.cubic_lattice_FEV();
		else{
			cout << "NNN = " << NNN << " not supported, setting NNN -> 6 and running normally." << endl;
			NNN = 6;
			cp.cubic_lattice_F();
		}
	}
	else if (NDIM == 2){
		if (NNN==4)
			cp.square_lattice_E();
		else if (NNN==8)
			cp.square_lattice_EV();
		else{
			cout << "NNN = " << NNN << " not supported, setting NNN -> 4 and running normally." << endl;
			NNN = 4;
			cp.square_lattice_E();
		}
	}
	

	// calc percolation threshold for single seed		

	// get rand sphere positions, radii
	cout << "* setting initial positions for N = " << N << " spheres with seed = " << seed << endl;
	cp.rand_sphere_pos(seed);
	cout << "getting cell positions..." << endl;
	cp.get_cell_positions();
	cout << "getting neighbor list for cells..." << endl;
	cp.get_cell_neighbors();
	cout << "making labels..." << endl;
	cp.label_cells();
	cout << "labelling lattice as void or not void..." << endl;
	cp.rand_sphere_void(aC);
	cp.print_pcell();
	cp.print_cellpos();	
	cp.print_cell_xyz(xyzobj);

	// loop over values, check for percolation
	cp.find_perc(aH,aL,epsilon);

	// once perc find completed, print info
	cp.print_void_xyz(xyzobj);
	cp.fprint_stats(statsobj,seed,1);
	cp.fprint_s(sobj);


	statsobj.close();
	sobj.close();
	xyzobj.close();
	return 0;
}