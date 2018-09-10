/*

	Methods for VOROPERC class

	BY Jack Treado, 09/07/2018

*/

#include "voro++.hh"
using namespace voro;

#include "voroperc.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

const double PI = 3.1415926;

voroperc::voroperc(int np) : voidcluster(np,1,3,6){
	NVp = nullptr;
	NEp = nullptr;
	NFp = nullptr;

	this->init_NVp(np);
	this->init_NEp(np);
	this->init_NFp(np);
}


voroperc::~voroperc(){
	// clear contents of vectors
	vpx.clear();
	vpy.clear();
	vpz.clear();

	epx.clear();
	epy.clear();
	epz.clear();

	cpax.clear();
	cpay.clear();
	cpaz.clear();

	// delete dynamically allocated memory
	delete [] NVp;
	delete [] NEp;
	delete [] NFp;

	// close stream objects
	if (xyzobj.is_open())
		xyzobj.close();
}




void voroperc::get_voro(int printit){
	// get variables from clustertree
	int NDIM;
	NDIM = this->get_NDIM();

	/******************************

		Generate Voronoi Diagram
		using voro++
	
	*******************************/

	// initialize voro++ loop variables
	// NOTE: naming conventions based on 
	// example file polygons.cc
	int i,j,id,n,init;
	double x,y,z;
	bool pbc;
	voronoicell_neighbor c;
	vector<double> v;

	// instantiate main object of container class
	n = 2;
	init = 5;
	pbc = true;
	container con(0,B[0],0,B[1],0,B[2],n,n,n,pbc,pbc,pbc,init);

	for (i=0; i<NP; i++)
		con.put(i,pos[i][0],pos[i][1],pos[i][2]);

	// loop over all particles, get information about cells
	c_loop_all cl(con);
	if (cl.start()) do if(con.compute_cell(c,cl)) {	
		cl.pos(x,y,z);
		id = cl.pid();	

		// Gather information about the computed voronoi cell
		c.vertices(x,y,z,v);

		// get number of vertices per particles
		this->set_NVp(id,v.size()/3);
		this->set_NEp(id,c.number_of_edges());
		this->set_NFp(id,c.number_of_faces());

		if (printit == 1){
			cout << endl << endl;
			cout << "on particle id = " << id;
			cout << "; x = " << x << ", pos[" << id << "][0] = " << pos[id][0];
			cout << "; y = " << y << ", pos[" << id << "][1] = " << pos[id][1];
			cout << "; z = " << z << ", pos[" << id << "][2] = " << pos[id][2];
			cout << endl;
			output_vpp_stats(c,x,y,z);
			cout << endl << endl;

			if (xyzobj.is_open())
				this->print_vertices_xyz(id,c);
		}

		// merge vertices from vertex vector to global vertex list
		// this->merge_vertices(id,c);

	} while(cl.inc());
}


double voroperc::merge_vertices(int id, voronoicell_neighbor& c){
	// get variables from clustertree
	int NDIM;
	NDIM = this->get_NDIM();

	// local variables
	vector<int> neigh,f_vert;
	vector<double> v;
	int i,k,nv;

	// get information about computed voronoi cell
	c.neighbors(neigh);
	c.face_vertices(f_vert);
	c.vertices(pos[id][0],pos[id][1],pos[id][2],v);

	// add unique vertices to list of vertices
	nv = v.size()/NDIM;
	k = 0;
	for (i=0; i<nv; i++){
		// get position of vertex
		xi = NDIM*i;	vx = v.at(xi);
		yi = NDIM*i+1;	vy = v.at(yi);
		zi = NDIM*i+2;	vz = v.at(zi);
		if (vpx.size()==0)
			this->add_vertex(vx,vy,vz);	
		else
			k = this->check_vertex(vx,vy,vz);

		// map particle and vertex number to k value
		k2p.at(k) = id;
		k2v.at(k) = i;
	}

}

void voroperc::add_vertex(double vx, double vy, double vz){
	vpx.push_back(vx);
	vpy.push_back(vy);
	vpz.push_back(vz);
}


double voroperc::check_vertex(double vx, double vy, double vz){
	// get variables from clustertree
	int NDIM;
	NDIM = this->get_NDIM();

	// local variables
	int i,k;
	double dx,dy,dz,dr;
	double v_ep = 2e-2;




}










// PRINTING INFO

// Voronoi vertex access
void output_vpp_stats(voronoicell_neighbor& v, double x, double y, double z){
	int nvertices = 0;
	vector<double> vec;
	v.vertices(x,y,z,vec);
	nvertices = vec.size();

	// Output vertex-based statistics
	cout << endl;
	cout << "Vertex global positions from vector : ";
	for (int i=0; i<nvertices; i+=3){
		cout << "(";
		cout << vec.at(i)-x << ", ";
		cout << vec.at(i+1)-y << ", ";
		cout << vec.at(i+2)-z << ") ";
	}
	cout << endl;
	cout << endl;
	printf("Total vertices      		: %d\n",nvertices);	
	printf("Vertex local positions    	: ");v.output_vertices();puts("");
	printf("Vertex orders      		  	: ");v.output_vertex_orders();puts("");
	printf("Max rad. sq. vertex 		: %g\n\n",0.25*v.max_radius_squared());

	// Output edge-based statistics
	printf("Total edges         		: %d\n",v.number_of_edges());
	printf("Total edge distance 		: %g\n",v.total_edge_distance());
	printf("Face perimeters     		: ");v.output_face_perimeters();puts("\n");

	// Output face-based statistics
	printf("Total faces         		: %d\n",v.number_of_faces());
	printf("Surface area        		: %g\n",v.surface_area());
	printf("Face freq. table    		: ");v.output_face_freq_table();puts("");
	printf("Face orders         		: ");v.output_face_orders();puts("");
	printf("Face areas          		: ");v.output_face_areas();puts("");
 	printf("Face normals        		: ");v.output_normals();puts("");
	printf("Face vertices       		: ");v.output_face_vertices();puts("\n");

	v.centroid(x,y,z);
	printf("Volume              : %g\n"
		   "Centroid vector     : (%g,%g,%g)\n",v.volume(),x,y,z);
}


// Print system to xyz
void voroperc::print_vertices_xyz(int id, voronoicell_neighbor& c){
	int i,d,w;
	int NDIM = this->get_NDIM();
	vector<double> v;
	w = 16;

	// get vertex positions for cell id
	c.vertices(pos[id][0],pos[id][1],pos[id][2],v);

	// print info to xyz string
	xyzobj << NVp[id] << endl;
	xyzobj << "Lattice=\"" << B[0] << " 0.0 0.0";
	xyzobj << " 0.0 " << B[1] << " 0.0";
	xyzobj << " 0.0 0.0 " << B[2] << "\"";
	xyzobj << '\t';
	xyzobj << "Properties=species:S:1:pos:R:3:radius:R:1" << endl;
	for (i=0; i<NVp[id]; i++){
		xyzobj << setw(w) << "X" << id;
		for (d=0; d<NDIM; d++)
			xyzobj << setw(w) << v.at(NDIM*i+d);
		xyzobj << setw(w) << 0.01;
		xyzobj << endl;
	}
}










