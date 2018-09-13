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
	int N = 500;
	int seed = 1;
	string xyzstr = "voro_perc.xyz";

	// instantiate object
	voroperc vpo(N);
	vpo.open_xyz(xyzstr);
	
	// initialize spheres with random positions
	vpo.rand_sphere_pos(seed);

	// get voronoi network (vertices & edges in global frame)
	vpo.get_voro();

	// find percolation
	double epsilon = 1e-6;
	double aH = 1;
	double aL = 1e-12;
	vpo.voro_edge_perc(epsilon,seed,aH,aL);

	// that's all folks
	return 0;
}