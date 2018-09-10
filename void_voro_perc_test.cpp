// Test voronoi void percolation

#include "voidcluster.h"
#include "voroperc.h"
#include "voro++.hh"
#include "clustertree.h"
#include <iostream>

int main(){
	// initialize the only two important parameters
	int N = 3;
	int seed = 1;

	// instantiate object
	voroperc vpo(N);
	
	// initialize spheres with random positions
	vpo.rand_sphere_pos(seed);

	// get voronoi, output v cell information
	vpo.get_voro(1);

	// that's all folks
	return 0;
}