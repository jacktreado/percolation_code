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
	int N = 20;
	int seed = 1;
	string xyzstr = "voro_perc.xyz";

	// instantiate object
	voroperc vpo(N);
	vpo.open_xyz(xyzstr);
	
	// initialize spheres with random positions
	vpo.rand_sphere_pos(seed);

	// get voronoi network (vertices & edges in global frame)
	vpo.get_voro(1);

	// setup lattice for edge percolation
	vpo.setup_edge_perc_lattice();

	// 

	// that's all folks
	return 0;
}