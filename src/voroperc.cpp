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

	// initialize size arrays
	this->init_NVp(np);
	this->init_NEp(np);
	this->init_NFp(np);

	// initialize vector arrays
	facen = new vector<int>[np];
	vface = new vector<int>*[np];
	vpx = new vector<double>[np];
	vpy = new vector<double>[np];
	vpz = new vector<double>[np];

	for (int i=0; i<NP; i++){
		NVp[i] = 0;
		NEp[i] = 0;
		NFp[i] = 0;
		vface[i] = nullptr;
	}
}


voroperc::~voroperc(){
	// clear contents of vectors
	vx.clear();
	vy.clear();
	vz.clear();

	ex.clear();
	ey.clear();
	ez.clear();

	cpax.clear();
	cpay.clear();
	cpaz.clear();
	
	// delete vector arrays
	int i,f;
	for (i=0; i<NP; i++){
		facen[i].clear();
		for (f=0; f<NFp[i]; f++)
			vface[i][f].clear();
		delete [] vface[i];
		vpx[i].clear();
		vpy[i].clear();
		vpz[i].clear();
	}
	delete [] facen;
	delete [] vface;
	delete [] vpx;
	delete [] vpy;
	delete [] vpz;
	delete [] NVp;
	delete [] NEp;
	delete [] NFp;

	// close stream objects
	if (xyzobj.is_open())
		xyzobj.close();
}


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
	int i,j,id,nvtmp,nftmp,netmp,n,init;
	double x,y,z;
	bool pbc;
	voronoicell_neighbor c;
	vector<int> neigh,f_vert;
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
		c.neighbors(neigh);
		c.face_vertices(f_vert);
		c.vertices(x,y,z,v);

		// get number of vertices/edges/faces per particles
		nvtmp = v.size()/3;
		netmp = c.number_of_edges();
		nftmp = c.number_of_faces();
		this->set_NVp(id,nvtmp);
		this->set_NEp(id,netmp);
		this->set_NFp(id,nftmp);

		if (printit == 1){
			cout << endl << endl;
			cout << "on particle id = " << id;
			cout << "; x = " << x << ", pos[" << id << "][0] = " << pos[id][0];
			cout << "; y = " << y << ", pos[" << id << "][1] = " << pos[id][1];
			cout << "; z = " << z << ", pos[" << id << "][2] = " << pos[id][2];
			cout << endl;
			cout << "printing neigh vector: " << endl;
			for (i=0; i<neigh.size(); i++)
				cout << setw(8) << neigh[i];
			cout << endl;
			cout << "printing f_vert vector: " << endl;
			for (i=0; i<f_vert.size(); i++)
				cout << setw(8) << f_vert[i];
			cout << endl;
			output_vpp_stats(c,x,y,z);
			cout << endl << endl;

			if (xyzobj.is_open())
				this->print_vertices_xyz(id,c);
		}

		// store particle vertex/face information
		cout << "creating face vectors..." << endl;		
		this->store_face_neighbors(id,neigh);
		this->store_face_vertices(id,f_vert);
		this->store_particle_vertices(id,v);
	} while(cl.inc());

	if (printit == 1){
		cout << endl << endl;
		this->print_face_vectors();
	}
}


// auxilliary functions for get_voro
void voroperc::store_face_neighbors(int id, vector<int>& neigh){
	facen[id].resize(neigh.size());
	facen[id] = neigh;
}

void voroperc::store_face_vertices(int id, vector<int>& f_vert){
	// local variables
	int i,j,f,v,k;
	int nf = f_vert.size();
	int nvf = 0;
	int vtmp = 0;
	int kmax = 2*nf;

	// initialize vface[id] to have NFp[id] faces 
	this->init_vface(id,NFp[id]);

	// loop over f_vert vector, store faces for each vertex v
	i = 0;
	f = 0;
	k = 0;
	while(i < nf && k < kmax){		
		// get number of vertices for face f		
		nvf = f_vert[i];
		cout << "nvf = " << nvf << endl;
		cout << "i = " << i << endl << endl;

		// loop over vertices on face f
		for (j=0; j<nvf; j++){
			// increment f_vert index
			i++;

			// get vertex number
			v = f_vert[i];

			// store vertex number in face position
			vface[id][f].push_back(v);
		}
		// increment i by 1 to get to next face
		i++;

		// increment face index
		f++;

		// increment failsafe k
		k++;
	}
	if (k == kmax){
		cout << "ERROR: store_face_vertices loop did not execute properly, exiting..." << endl;
		throw;
	}
}

void voroperc::store_particle_vertices(int id, vector<double>& v){
	// get variables from clustertree
	int NDIM;
	NDIM = this->get_NDIM();

	// local variables
	int k,nvtmp,kxi,kyi,kzi;
	double vxtmp,vytmp,vztmp;
	nvtmp = v.size()/NDIM;

	// store neighbors for each face on each particle
	for (k=0; k<nvtmp; k++){
		kxi = NDIM*k;   vxtmp = v[kxi];
		kyi = NDIM*k+1;	vytmp = v[kyi];
		kzi = NDIM*k+2;	vztmp = v[kzi];
	
		if ( (k == 0 && vpx[id].size() == 0) || k > 0 ){
			vpx[id].push_back(vxtmp);
			vpy[id].push_back(vytmp);
			vpz[id].push_back(vztmp);
		}
		else{
			cout << "trying to push_back vertices on pre-allocated vector vpx, ending..." << endl;
			throw;
		}
	}
}


// PRINTING INFO


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

// Print face vector info
void voroperc::print_face_vectors(){
	int i,j,k;

	cout << endl;
	cout << "Face neighbors: " << endl;
	for (i=0; i<NP; i++){
		cout << "facen[" << i << "] = ";
		for (j=0; j<facen[i].size(); j++)
			cout << setw(6) << facen[i][j];
		cout << endl;
	}
	cout << endl;
	cout << "Vertex face info: " << endl;
	for (i=0; i<NP; i++){
		cout << "NFp[" << i << "] = " << NFp[i] << ", vface[" << i  << "]:" << endl;
		for (j=0; j<NFp[i]; j++){
			cout << "** face = " << j << ": ";
			for (k=0; k<vface[i][j].size(); k++)
				cout << setw(6) << vface[i][j][k];
			cout << endl;
		}
		cout << endl;
	}
	cout << endl;
	cout << "Vertices per particle: " << endl;
	for (i=0; i<NP; i++){
		cout << endl;
		cout << "NVp[" << i << "] = " << NVp[i] << endl;
		cout << "vpx[" << i <<  "]:";
		for (j=0; j<NVp[i]; j++)
			cout << setw(12) << vpx[i][j];
		cout << endl;
		cout << "vpy[" << i <<  "]:";
		for (j=0; j<NVp[i]; j++)
			cout << setw(12) << vpy[i][j];
		cout << endl;
		cout << "vpz[" << i <<  "]:";
		for (j=0; j<NVp[i]; j++)
			cout << setw(12) << vpz[i][j];
		cout << endl;
		cout << endl;
	}
}










