// Test voronoi void percolation

#include "voidcluster.h"
#include "voroperc.h"
#include "voro++.hh"
#include "clustertree.h"
#include <iostream>
#include <string>

using namespace std;

int main(){
	// initialize the only two important parameters
	int N = 1000;
	int seed = 12;
	string xyzstr = "voro_perc.xyz";

	// instantiate object
	voroperc vpo(N);
	vpo.open_xyz(xyzstr);
	
	// initialize spheres with random positions
	vpo.rand_sphere_pos(seed);

	// get voronoi network (vertices & edges in global frame)
	vpo.get_voro();

	// find percolation
	double epsilon = 1e-8;
	double aH = 1;
	double aL = 1e-12;
	vpo.voro_edge_perc(epsilon,seed,aH,aL);


	// test on single radius
	// double a,poro;
	// double PI = 3.1415296;
	// poro = 0.035;
	// a = pow(abs(3*log(poro)/(4*PI*N)),0.333333333333333);

	// cout << "setting up system..." << endl;
	// vpo.reset_sys();
	// vpo.reset_ptr();
	// vpo.set_radius(r);

	// cout << "get particle/cpa intersections..." << endl;
	// vpo.set_lattice_cpa_intersect();

	// that's all folks
	return 0;
}