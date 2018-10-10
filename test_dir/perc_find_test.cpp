// Code to test percolation finder

#include <iostream>
#include <cmath>
#include "clustertree.h"

int main(){

	cout << "Testing perc search..." << endl;

	// declare box info
	int L = 50;
	int NDIM = 2;
	int NNN = 4;
	int i,j,perc,N;
	int Lhalf = round(0.5*L);
	double p = 0.56;

	// instantiate object of clustertree
	clustertree ct(L,NDIM,NNN);

	// setup NNN
	ct.square_lattice_E();

	// setup lattice
	N = ct.get_N();

	// XY perc
	cout << "==== TESTING SELF-WRAPPING perc..." << endl;

	ct.rand_site(p,1);
	for (i=0; i<L; i++){
		j = i*(L*L)+L*i+i;
		if (i==1){
			ct.set_lattice(j+1,1);		
			ct.set_lattice(j-1,1);
		}
		if (i >= Lhalf)
			ct.set_lattice(j-i+L-2,1);
		else
			ct.set_lattice(j,1);

		cout << "i = " << i << "; j = " << j << endl;
	}	

	// get whether or not percolation
	cout << "merging clusters..." << endl;
	ct.merge_clusters();
	cout << endl;
	cout << "max cluster = " << ct.get_pclus() << endl;
	ct.post_process();
	ct.print_rl_lattice();

	perc = ct.get_perc();

	cout << "perc = " << perc << endl;
	ct.print_cluster_stats(0,1);
	ct.reset_sys();

	return 0;
}