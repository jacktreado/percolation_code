/*

	Test voidcluster for N overlapping spheres in NDIM

*/

#include "clustertree.h"
#include "voidcluster.h"
#include <iostream>
#include <fstream>
#include <cmath>

#define L 1
#define NDIM 3
#define NNN 6
#define N 20
#define seed 1

const double PI = 3.1415269;


int main(){
	// percolation variables
	double poroH,aH,poroL,aL,poroC,aC,epsilon;
	poroH = 0.05;
	aL = pow(abs(3*log(poroH)/(4*PI*N)),0.333333333333333);

	poroL = 0.02;
	aH = pow(abs(3*log(poroL)/(4*PI*N)),0.333333333333333);

	poroC = 0.03;
	aC = pow(abs(3*log(poroC)/(4*PI*N)),0.333333333333333);

	epsilon = 1e-6;

	cout << "poroL = " << poroL << "; aH = " << aH << endl;
	cout << "poroH = " << poroH << "; aL = " << aL << endl;
	cout << "poroC = " << poroC << "; aC = " << aC << endl;


	double g0 = 0.1;
	double g = g0*aC;
	int NS;			

	// figure out number of lattice sites for given g
	NS = round(L/g);

	// double g;
	// int NS;

	// NS = 50;
	// g = (double)L/(NS-1);

	cout << "NS = " << NS << "; g = " << g << endl;

	// initialize voidcluster object
	voidcluster sv(N,NS,NDIM,NNN);

	// get percolation
	cout << "* setting initial positions for N = " << N << " spheres with seed = " << seed << endl;
	sv.cubic_lattice_FEV();
	sv.rand_sphere_pos(seed);
	sv.rand_sphere_void(aC);

	cout << "lattice sum before perc = " << sv.get_lattice_sum() << endl;

	// if (NS <= 30)
		// sv.print_labels();

	cout << "* finding percolation..." << endl;
	sv.find_perc(aH,aL,epsilon);

	// if (NS <= 30)
		// sv.print_labels();

	// sv.print_lattice();
	// sv.print_labels();
	int lattice_tot = sv.get_lattice_sum();
	// lattice_tot = (double)lattice_tot*g*g;
	double lsites = (double)lattice_tot;
	lsites = lsites*pow(g,NDIM);
	cout << "* outside, from lattice sum, lsum = " << lattice_tot << "; vctot = " << lsites << "; g^NDIM = " << pow(g,NDIM) << endl;

	// open output files, print
	ofstream sf("test_s.dat");
	ofstream statf("test_stat.dat");

	sv.fprint_s(sf);
	sv.fprint_stats(statf,seed,1);

	sf.close();
	statf.close();

	return 0;
}
