// Code to test percolation finder

#include <iostream>
#include <cmath>
#include "clustertree.h"

int main(){

	cout << "Testing perc search..." << endl;

	// declare box info
	int L = 7;
	int NDIM = 3;
	int NNN = 18;
	int i,j,perc,N;
	int Lhalf = round(0.5*L);

	// instantiate object of clustertree
	clustertree ct(L,NDIM,NNN);

	// setup NNN
	ct.cubic_lattice_FEV();

	// setup lattice
	N = ct.get_N();

	// XY perc
	cout << "==== TESTING XY perc..." << endl;

	ct.rand_site(0,1);
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
	ct.merge_clusters();
	cout << "max cluster = " << ct.get_pclus() << endl;
	ct.post_process();
	ct.print_rl_lattice();

	perc = ct.get_perc();

	cout << "perc = " << perc << endl;
	ct.print_cluster_stats(0,1);
	ct.reset_sys();


	// YZ perc
	cout << "==== TESTING XZ perc..." << endl;

	ct.rand_site(0,1);
	for (i=0; i<L; i++){
		j = L*i;
		if (i == 1){
			ct.set_lattice(j+1,1);
			ct.set_lattice(j+2,1);
		}
		if (i != Lhalf)
			ct.set_lattice(j,1);		
		cout << "i = " << i << "; j = " << j << endl;
	}	

	// get whether or not percolation
	ct.merge_clusters();
	cout << "max cluster = " << ct.get_pclus() << endl;
	ct.post_process();
	ct.print_rl_lattice();

	perc = ct.get_perc();

	cout << "perc = " << perc << endl;
	ct.print_cluster_stats(0,1);
	ct.reset_sys();


	// XZ perc
	cout << "==== TESTING YZ perc..." << endl;

	ct.rand_site(0,1);
	int xzi = Lhalf*(L*L) + 2*L + Lhalf;
	for (i=0; i<L; i++){
		j = xzi + i;
		ct.set_lattice(j,1);
		cout << "i = " << i << "; j = " << j << endl;
	}

	// get whether or not percolation
	ct.merge_clusters();
	cout << "max cluster = " << ct.get_pclus() << endl;
	ct.post_process();
	ct.print_rl_lattice();

	perc = ct.get_perc();

	cout << "perc = " << perc << endl;
	ct.print_cluster_stats(0,1);
	ct.reset_sys();

	

	return 0;
}