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

	// set other arrays to nullptr
	vneigh = nullptr;
	vp = nullptr;
	eneigh = nullptr;
	ep = nullptr;
}


voroperc::~voroperc(){
	cout << "in voroperc constructor" << endl;
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

	// delete face pair map
	int nfp = vmap.size();
	for (i=0; i<nfp; i++){
		delete [] vmap[i];
		vmap[i] = nullptr;
	}
	vmap.clear();

	// delete v neighbor map
	for (i=0; i<NV; i++){
		vn_map[i].clear();
		delete [] eij[i];
		vneigh[i].clear();
		vp[i].clear();
	}
	vn_map.clear();
	delete [] eij;
	delete [] vneigh;
	delete [] vp;

	// delete edge objects
	for (i=0; i<NE; i++){
		eneigh[i].clear();
		ep[i].clear();
	}
	delete [] eneigh;
	delete [] ep;

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

	// merge vertices based on faces
	this->merge_vertices();
	this->update_NV();
	if (printit == 1){
		cout << endl << endl;
		this->print_global_vertices();
	}

	// initialize vertex neighbor array, populate using vn_map	
	this->update_vertex_neighbors();
	if (printit == 1){
		cout << endl << endl;
		this->print_vertex_neighbors();

		cout << "Printing vp: " << endl;
		for (i=0; i<NV; i++){
			cout << "vp[" << i << "] : ";
			for (j=0; j<vp[i].size(); j++)
				cout << setw(8) << vp[i][j];
			cout << endl;
		}
	}

	// get edge positions and neighbors based on vertex neighbors
	this->get_edges();
	if (printit == 1){
		cout << endl << endl;
		this->print_edge_neighbors();

		cout << "Printing ep: " << endl;
		for (i=0; i<NE; i++){
			cout << "ep[" << i << "] : ";
			for (j=0; j<ep[i].size(); j++)
				cout << setw(8) << ep[i][j];
			cout << endl;
		}
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


// MERGE VERTICES FUNCTION
void voroperc::merge_vertices(){
	// local variables
	int i,f,g;
	int fnp,fnf,red;

	// loop over particles, add vertices to global list
	for (i=0; i<NP; i++){
		// loop over faces, get neighboring faces
		for (f=0; f<NFp[i]; f++){
			// merge face f with face fnf on neighboring particle fnp

			// check whether particle i shares multiple faces with particle g
			g = this->check_multiple_faces(i,f);
			if (g == -1){
				fnp = facen[i][f];
				fnf = this->get_neighboring_face(i,f,fnp);				
			}
			else if (g >= 0){
				fnp = g;
				fnf = this->get_closest_face(i,f,fnp);
			}
			else {
				cout << "g returned negative in merge_vertices at i = " << i << ", f = " << f << ", ending..." << endl;
				throw;
			}

			// check if (i,f) & (fnp,fnf) pairing is redundant
			if (vmap.size() == 0)
				red = 0;
			else
				red = this->check_redundancy(i,f,fnp,fnf);

			// if not redundant, then combine vertices on (i,f) and (fnp,fnf)
			if (red == 0)
				this->combine_faces(i,f,fnp,fnf);
			else
				this->add_vertex_neighbors(i,f);
		}	
	}
}

// auxilliary functions for merge_vertices
int voroperc::check_multiple_faces(int i, int f){
	// local variables
	int j,g,gtmp;

	// let g be the neighbor of face f
	g = facen[i][f];

	// loop over faces on i, check to see if any other faces adjacent to g
	for (j=0; j<NFp[i]; j++){
		if (j != f){
			gtmp = facen[i][j];

			// if gtmp = g, at least duplicate, so return
			if (gtmp == g)
				return g;
		}
	}

	// if you get this far, no duplicates
	g = -1;
	return g;
}

int voroperc::get_neighboring_face(int i, int f, int fnp){
	// local variables
	int k;

	// loop over faces on fnp, check adjacency to particle i			
	for (k=0; k<NFp[fnp]; k++){
		if (facen[fnp][k] == i)
			return k;
	}

	// if you get this far, throw an error
	cout << "ERROR: no adjacency found on face opposite i = " << i << ", f = " << f << "...ending." << endl;
	throw;
}

int voroperc::get_closest_face(int i, int f, int fnp){
	// there are multiple adjacencies for particle i that are on particle fnp,
	// so loop over multiples and check which one is closest to face (i,f)
	vector<double> f1com,f2com;
	int j,nadj,mindistj;
	double dx,dy,dz,dr,mindist,mindist0;

	// get com of f1
	this->calc_face_com(i,f,f1com);

	// loop over adjacents to faces on fnp, determine min distance to (i,f)
	mindist0 = 10*B[0];
	mindist = mindist0;
	nadj = facen[fnp].size();
	for (j=0; j<nadj; j++){
		if (facen[fnp][j]==i){			
			// get face com of (fnp,j)
			this->calc_face_com(fnp,j,f2com);

			// get x distance
			dx = f2com[0]-f1com[0];
			dx = dx - B[0]*round(dx/B[0]);

			// get x distance
			dy = f2com[1]-f1com[1];
			dy = dy - B[1]*round(dy/B[1]);

			// get z distance
			dz = f2com[2]-f1com[2];
			dz = dz - B[2]*round(dz/B[2]);

			// get total distance
			dr = sqrt(dx*dx + dy*dy + dz*dz);

			// if distance is less than min distance, update min distance
			if (dr < mindist){
				mindist = dr;
				mindistj = j;
			}
		}
	}

	if (mindist < mindist0)
		return mindistj;
	else{
		cout << "ERROR: mindist not updated in get_closest_face, ending..." << endl;
		throw;
	}
}

void voroperc::calc_face_com(int i, int f, vector<double>& com){
	// get variables from clustertree
	int NDIM;
	NDIM = this->get_NDIM();

	// resize com
	com.resize(NDIM);

	// local variables
	int j,vi,nv;
	double xsum,ysum,zsum;

	// calculate avg x,y,z pos
	xsum = 0;
	ysum = 0;
	zsum = 0;
	nv = vface[i][f].size();
	for (j=0; j<nv; j++){
		vi = vface[i][f][j];
		xsum += vpx[i][vi];
		ysum += vpy[i][vi];
		zsum += vpz[i][vi];
	}
	xsum /= nv;
	ysum /= nv;
	zsum /= nv;

	// store com
	com[0] = xsum;
	com[1] = ysum;
	com[2] = zsum;
}

int voroperc::check_redundancy(int i, int f, int fnp, int fnf){
	int nfp,itmp,ftmp,jtmp,gtmp;
	bool revchecked, thischecked;

	// initialize booleans to false
	revchecked = false;
	thischecked = false;

	// loop over prior matches, check for repeats
	nfp = vmap.size();
	for (i=0; i<nfp; i++){
		// get map values
		itmp = vmap[i][0];
		ftmp = vmap[i][1];
		jtmp = vmap[i][2];
		gtmp = vmap[i][3];

		// check redundancy
		revchecked = (i==jtmp && f==gtmp && fnp==itmp && fnf==ftmp);
		thischecked = (i==itmp && f==ftmp && fnp==jtmp && fnf==gtmp);

		if (revchecked)
			return 1;
		else if (thischecked)
			return 1;
	}

	// if you get this far, no redundancies
	return 0;
}

void voroperc::combine_faces(int i, int f, int fnp, int fnf){
	// loop over vertices, add vertices to global list, make note
	// of which vertices are comprised of which faces
	int k,nv,nv1,nv2,nfp,vi,merge_vertex,vm,vtrue;
	int* p;
	vector<int> vtmp;

	// verify that number of vertices are same
	nv1 = vface[i][f].size();
	nv2 = vface[fnp][fnf].size();

	if (nv1 != nv2){
		cout << "number of vertices do not match between i = " << i << ", f = " << f;
		cout << " and fnp = " << fnp << ", fnf = " << fnf << endl;
		cout << "Ending..." << endl;
		throw;
	}

	// count number of merged vertices
	vm = 0;

	// if numbers are the same, add to vertex list from face (i,f)	
	for (k=0; k<nv1; k++){
		vi = vface[i][f][k];
		merge_vertex = this->check_vertex_history(i,f,vi);
		if (merge_vertex == 1){
			// check to see if already on master list
			vtrue = this->get_vtrue(i,vi);

			if (vtrue == -1){
				// add vertices to global list
				vx.push_back(vpx[i][vi]);
				vy.push_back(vpy[i][vi]);
				vz.push_back(vpz[i][vi]);
				vtrue = vx.size()-1;

				// incrememnt vm
				vm++;

				// add another row to vn_map
				vn_map.push_back(vtmp);			

				// add vertex neighbors
				this->add_vertex_neighbors(vtrue,vi,i,f);

				if (vn_map[vtrue].size() == 0){
					cout << "* neighbors not added at vertex " << vtrue << " : "; 
					cout << "i=" << i << ", f=" << f << ", v=" << vi << endl;
					throw;
				}
			}
			else
				this->add_vertex_neighbors(vtrue,vi,i,f);
		}
		else{
			vtrue = this->get_vtrue(i,vi);
			this->add_vertex_neighbors(vtrue,vi,i,f);
		}
	}

	// add vertex neighbors from other face
	this->add_vertex_neighbors(fnp,fnf);

	// add info about faces if at least 1 vertex was merged
	if (vm > 0){
		nfp = vmap.size();			// number of face pairs already in vector
		vmap.push_back(p);			// push pointer back to vector, increase size by nfp+1
		vmap[nfp] = new int[4];		// because index from 0, vmap[nfp] needs to be initialized
		vmap[nfp][0] = i;			// first column: particle i
		vmap[nfp][1] = f;			// second column: face f
		vmap[nfp][2] = fnp;			// third column: particle fnp
		vmap[nfp][3] = fnf;			// fourth column: face fnf
	}
}

int voroperc::check_vertex_history(int i, int f, int vi){
	// local variables
	int j,k,nv,nfp,ftmp;
	bool already_merged;
	vector<int> facelist;

	// get all faces that vi is a part of on particle i
	for (j=0; j<NFp[i]; j++){
		// do not check f
		if (j != f){
			// loop over vertices on face j
			nv = vface[i][j].size();			
			for (k=0; k<nv; k++){
				// if vfacen[i][j][k] = vi, add to face list
				if (vface[i][j][k] == vi)
					facelist.push_back(j);
			}
		}
	}

	// throw error if facelist not populated
	if (facelist.size() == 0){
		cout << "error populating list of faces with vertex " << vi << "on (" << i << "," << f << ") ending..." << endl;
		throw;
	}

	// check facelist to see if any faces have been merged previously
	nv = facelist.size();
	for (k=0; k<nv; k++){
		// get face k that vi is a part of
		ftmp = facelist[k];

		// get number of face pairs in map
		nfp = vmap.size();

		// loop over map, check for (i,ftmp) pairs in either half of matrix
		for (j=0; j<nfp; j++){
			already_merged = (vmap[j][0] == i && vmap[j][1] == ftmp);
			already_merged = already_merged || (vmap[j][2] == i && vmap[j][3] == ftmp);
			if (already_merged)
				return 0;
		}
	}

	// if you made it this far, none of the other faces with vertex vi have been merged,
	// so return 1 to merge
	return 1;
}

int voroperc::get_vtrue(int i, int vi){
	// local variables
	double vix,viy,viz,dx,dy,dz,dr,mindist0,mindist;
	int j,vtrue,nv;

	vix = vpx[i][vi];
	viy = vpy[i][vi];
	viz = vpz[i][vi];

	// loop over true vertices, get min dist
	mindist0 = 10*B[0];
	mindist = mindist0;
	nv = vx.size();
	for (j=0; j<nv; j++){
		// get x distance
		dx = vx[j]-vix;
		dx = dx - B[0]*round(dx/B[0]);

		// get x distance
		dy = vy[j]-viy;
		dy = dy - B[1]*round(dy/B[1]);

		// get z distance
		dz = vz[j]-viz;
		dz = dz - B[2]*round(dz/B[2]);

		// get total distance
		dr = sqrt(dx*dx + dy*dy + dz*dz);

		if (dr < mindist){
			mindist = dr;
			vtrue = j;
		}
	}

	if (abs(mindist-mindist0) < 1e-15){
		// cout << "ERROR: true vertex not found for i = " << i << ", vi = " << vi << "; check that it exists." << endl;
		return -1;
	}
	else if (mindist > 1e-14){
		// cout << "ERROR: true vertex not very close to to i = " << i << ", vi = " << vi << "; mindist = " << mindist << endl;
		return -1;
	}
	else{
		return vtrue;
	}
}

void voroperc::add_vertex_neighbors(int vtrue, int vi, int i, int f){
	// local variables
	int j,nfv,fwd,bwd,vfwd,vbwd;

	// loop over face (i,f), add neighbors of vi on (i,f)
	// to vn_map[vtrue]
	nfv = vface[i][f].size();
	for (j=0; j<nfv; j++){
		if (vface[i][f][j]==vi){
			// get fwd and bwd indices
			fwd = j + 1;
			if (fwd == nfv)
				fwd = 0;
			bwd = j - 1;
			if (bwd == -1)
				bwd = nfv-1;

			// get neighboring vertices
			vfwd = vface[i][f][fwd];
			vbwd = vface[i][f][bwd];

			// add to vn_map
			vn_map[vtrue].push_back(i);
			vn_map[vtrue].push_back(vbwd);
			vn_map[vtrue].push_back(i);
			vn_map[vtrue].push_back(vfwd);
			break;
		}
	}
}

void voroperc::add_vertex_neighbors(int i, int f){
	int j,nv,vj,vtrue;

	// loop over all vertices on a given face, get vtrue, add neighbors to vn_map
	nv = vface[i][f].size();
	for (j=0; j<nv; j++){
		vj = vface[i][f][j];
		vtrue = get_vtrue(i,vj);
		this->add_vertex_neighbors(vtrue,vj,i,f);
	}
}


// UPDATE VERTEX NEIGHBOR FUNCTION
void voroperc::update_vertex_neighbors(){
	// intialize array vneigh
	this->init_vneigh();

	// initialize array vp
	this->init_vp();

	// local variables
	int i,j,nv,cell,vertex,vtrue;
	vector<int> vlist;

	// loop over vn_map, turn face-frame vertices into global-frame vertices
	for (i=0; i<NV; i++){
		// loop over vertices connected to i	
		vlist = vn_map[i];
		nv = vlist.size();
		for (j=0; j<nv; j+=2){
			// get cell and vertex connected to i
			cell = vlist[j];
			vertex = vlist[j+1];

			// get true vertex
			vtrue = this->get_vtrue(cell,vertex);

			// decide to add to vneigh or not
			this->add_unique_to_vneigh(i,vtrue);

			// decide to add to vp list or not
			this->add_unique_to_vp(cell,vtrue);
		}
	}
}

void voroperc::add_unique_to_vneigh(int i, int vtrue){
	int j,nv,skip;
	nv = vneigh[i].size();

	skip = 0;
	if (nv == 0)
		// if vneigh[i] empty, add vtrue
		vneigh[i].push_back(vtrue);
	else{
		// check if vtrue already there
		for (j=0; j<nv; j++){
			if (vneigh[i][j] == vtrue){
				skip = 1;
				break;
			}
		}

		// if skip = 0, add vtrue
		if (skip == 0)
			vneigh[i].push_back(vtrue);
	}
}

void voroperc::add_unique_to_vp(int cell, int vtrue){
	int j,np,skip;
	np = vp[vtrue].size();

	skip = 0;
	if (np == 0)
		// if vp[cell] empty, add vtrue
		vp[vtrue].push_back(cell);
	else{
		// check if vtrue already there
		for (j=0; j<np; j++){
			if (vp[vtrue][j] == cell){
				skip = 1;
				break;
			}
		}

		// if skip = 0, add vtrue
		if (skip == 0)
			vp[vtrue].push_back(cell);
	}
}


// GET EDGE POSITIONS AND NEIGHBORS
void voroperc::get_edges(){
	int i,j,v,e,nvn;

	// initialize edge def array
	eij = new int*[NV];
	for (i=0; i<NV; i++){
		eij[i] = new int[NV];
		for (j=0; j<NV; j++)
			eij[i][j] = -1;
	}

	e = 0;
	for (i=0; i<NV; i++){
		// get number of vertex neighbors
		nvn = vneigh[i].size();

		// loop over neighbors
		for (j=0; j<nvn; j++){
			v = vneigh[i][j];		
			if (v > i){
				// ex.push_back(0.5*(vx[i]+vx[j]));
				// ey.push_back(0.5*(vy[i]+vy[j]));
				// ez.push_back(0.5*(vz[i]+vz[j]));	
				// calculate edge position (relative to i)
				this->get_edge_positions(i,v);
				eij[i][v] = e;
				eij[v][i] = e;
				e++;
			}
		}
	}

	// get NE
	NE = ex.size();
	eneigh = new vector<int>[NE];	
	ep = new vector<int>[NE];

	// get edge neighbors
	for (i=0; i<NV; i++){
		nvn = vneigh[i].size();
		for (j=0; j<nvn; j++){
			v = vneigh[i][j];
			if (v > i){
				e = eij[i][v];
				this->get_edge_neighbors(e,i,v);
				this->get_ep(e,i,v);
			}
		}
	}
}

void voroperc::get_edge_positions(int vi, int vj){
	double dx,dy,dz;
	double vix,viy,viz;
	double vjx,vjy,vjz;

	// get edge position taking PBCs into account
	vix = vx[vi];	vjx = vx[vj];
	viy = vy[vi];	vjy = vy[vj];
	viz = vz[vi];	vjz = vy[vj];

	// check if vertices are too far away

	// x direction 
	if (abs(vix-vjx)>0.5*B[0]){
		if (vjx > vix)
			vjx -= B[0];
		else
			vjx += B[0];
	}

	// x direction
	if (abs(viy-vjy)>0.5*B[1]){
		if (vjy > viy)
			vjy -= B[1];
		else
			vjy += B[1];
	}

	// x direction
	if (abs(viz-vjz)>0.5*B[2]){
		if (vjz > viz)
			vjz -= B[2];
		else
			vjz += B[2];
	}

	// get edge positions from translated vertices
	ex.push_back(0.5*(vix+vjx));
	ey.push_back(0.5*(viy+vjy));
	ez.push_back(0.5*(viz+vjz));
}

void voroperc::get_edge_neighbors(int e, int i, int j){
	int k,n;

	// run over neighbors of i and j, get vertices, add to eneigh[e]
	for (k=0; k<vneigh[i].size(); k++){
		n = vneigh[i][k];
		if (eij[i][n] != e)
			eneigh[e].push_back(eij[i][n]);
	}
	for (k=0; k<vneigh[j].size(); k++){
		n = vneigh[j][k];
		if (eij[j][n] != e)
			eneigh[e].push_back(eij[j][n]);
	}
}

void voroperc::get_ep(int e, int i, int j){
	int k,l,p1,p2;

	// get the union of particles adjacent to vertices i and j
	for (k=0; k<vp[i].size(); k++){
		p1 = vp[i][k];
		for (l=0; l<vp[j].size(); l++){
			p2 = vp[j][l];
			if (p2 == p1){
				ep[e].push_back(p1);
				break;
			}
		}
	}
}


// LOOP OVER RADII, CHECK FOR PERCOLATION
void voroperc::voro_edge_perc(double epsilon, int seed, double aH, double aL){
	// set lattice up for edge percolation
	this->setup_edge_perc_lattice();

	// get NDIM from clustertree
	int NDIM = this->get_NDIM();

	// bisection variables
	double check,check_new,r;
	int k,kmax;
	kmax = 1e4;
	k = 0;
	check = 10*epsilon;
	check_new = check;
	r = 0.5*(aH+aL);
	int perc = 0;
	string percdir = "NA";

	// get vector edge_nn, for cluster merging
	vector<int> edge_nn(NE);
	for (int e=0; e<NE; e++)
		edge_nn[e] = eneigh[e].size();	

	// output information to console
	int w = 15;
	int prc = 5;
	double poro; 
	cout << "Looping over sphere radii, looking for edge percolation" << endl;
	cout << "aH = " << aH << endl;
	cout << "aL = " << aL << endl;
	cout << "a = " << r << endl;
	cout << setw(w) << "k";
	cout << setw(w) << "r";
	cout << setw(w) << "aH";
	cout << setw(w) << "aL";
	cout << setw(w) << "check";
	cout << setw(w) << "~poro";
	cout << setw(w) << "perc";
	cout << setw(w) << "dir";
	cout << setw(w) << "cnum";
	cout << setw(w) << "smax";
	cout << setw(w) << "lsites";
	cout << setw(w) << "fcalls";
	cout << endl;
	for (int s=0; s<12*w; s++)
		cout << "=";
	cout << endl;
	
	while( ((check > epsilon) || perc == 0) && k < kmax){
		check = r;

		this->reset_sys();
		this->reset_ptr();
		this->set_radius(r);				
		percdir = "NA";

		// check intersection by edge midpoint
		this->set_lattice_midpoint_intersect();

		// check point of closest approach intersection
		// NOTE: NEED TO CALCULATE CPA FIRST!
		// this->check_cpa_intersect();

		// merge clusters		
		this->merge_clusters(edge_nn);
		this->post_process_voro();
		this->check_voro_percolation(percdir);		

		// update check based on perc or not
		if (perc == 1)
			aL = r;
		else
			aH = r;

		r = 0.5*(aH+aL);
		check = abs(check-r)/r;
		check_new = check;
		k++;

		// print info to console
		perc = this->get_perc();
		poro = exp(-NP*(1.33333333333)*PI*pow(r,NDIM));
		cout << setprecision(prc) << setw(w) << k;
		cout << setprecision(prc) << setw(w) << r;
		cout << setprecision(prc) << setw(w) << aH;
		cout << setprecision(prc) << setw(w) << aL;
		cout << setprecision(prc) << setw(w) << check_new;
		cout << setprecision(prc) << setw(w) << poro;
		cout << setprecision(prc) << setw(w) << perc;
		cout << setprecision(prc) << setw(w) << percdir;
		cout << setprecision(prc) << setw(w) << this->get_cnum();
		cout << setprecision(prc) << setw(w) << this->get_smax();
		cout << setprecision(prc) << setw(w) << this->get_lattice_sum();
		cout << setprecision(prc) << setw(w) << this->get_fcalls();
		cout << endl;
	}
}


// SET LATTICE UP FOR EDGE PERCOLATION
void voroperc::setup_edge_perc_lattice(){
	// local variables
	int e,f,enn;

	// Reset lattice with total number of edges
	cout << "resetting percolation lattice with NE = " << NE << " edges..." << endl;
 	this->reset_lattice(NE);

 	// Get nearest neighbor info for each lattice
 	cout << "getting nearest neighbor info..." << endl;
 	for (e=0; e<NE; e++){
 		enn = eneigh[e].size();
 		this->initialize_nn(e,enn);
 		for (f=0; f<enn; f++)
 			this->set_nn(e,f,eneigh[e][f]);
 	}

 	// Also get edge com, for perc checking
 	ecomx = 0;
 	ecomy = 0;
 	ecomz = 0;
 	for (e=0; e<NE; e++){
 		ecomx += ex[e];
 		ecomy += ey[e];
 		ecomz += ez[e];
 	}
 	ecomx /= NE;
 	ecomy /= NE;
 	ecomz /= NE;
}

void voroperc::set_lattice_midpoint_intersect(){
	double dx,dy,dz,dr;
	int e,p,pn,edgeoverlap;

	// loop over edges, see if overlapping with any particles
	for (e=0; e<NE; e++){
		dr = 0;
		edgeoverlap = 0;

		// loop over particles adjacent to e
		for (p=0; p<ep[e].size(); p++){
			// get particle
			pn = ep[e][p];

			// check distance
			dx = ex[e] - pos[pn][0];
			dx = dx - round(dx/B[0])*B[0];

			dy = ey[e] - pos[pn][1];
			dy = dy - round(dy/B[1])*B[1];

			dz = ez[e] - pos[pn][2];
			dz = dz - round(dz/B[2])*B[2];

			dr = sqrt(dx*dx + dy*dy + dz*dz);
			if (dr < rad[pn]){
				edgeoverlap = 1;
				break;
			}
		}

		if (edgeoverlap == 1)
			this->set_lattice(e,0);
		else
			this->set_lattice(e,1);
	}
}

void voroperc::check_voro_percolation(string& percdir){
	// local variables
	int e;
	double ex0,ey0,ez0;

	// perc check variables
	int cf = 1;
	double subdiv = 10.0;
	double loc,dloc;
	int perc = 0;

	// check percolation in X direction (post-process)
	perc = 0;
	cf = 1;
	loc = 0;
	dloc = B[0]/subdiv;
	while(cf == 1 && loc < B[0]-dloc){
		cf = 0;

		// check for all points between loc & loc+dloc
		for (e=0; e<NE; e++){

			// get edge positions (with translated edges)
			ex0 = ex[e]-ecomx+0.5;
			if (ex0 < 0)
				ex0 += B[0];
			else if (ex0 > B[0])
				ex0 -= B[0];

			// check for ex0 in x bin
			if (ex0 > loc && ex0 <= loc+dloc){
				// check if edge is occupied with biggest cluster
				if (this->findroot(e)==this->get_pclus()){					
					cf = 1;
					break;
				}	
			}
		}
		if (cf == 1)
			loc += dloc;			
	}
	if (cf == 1){
		perc = 1;
		percdir = "x";
	}

	// check percolation in Y direction (post-process)
	if (perc == 0){
		cf = 1;
		loc = 0;
		dloc = B[1]/subdiv;
		while(cf == 1 && loc < B[1]-dloc){
			cf = 0;

			// check for all points between loc & loc+dloc
			for (e=0; e<NE; e++){

				// get edge positions (with translated edges)
				ey0 = ey[e]-ecomy+0.5;
				if (ey0 < 0)
					ey0 += B[1];
				else if (ey0 > B[1])
					ey0 -= B[1];

				// check for ex0 in x bin
				if (ey0 > loc && ey0 <= loc+dloc){
					// check if edge is occupied with biggest cluster
					if (this->findroot(e)==this->get_pclus()){					
						cf = 1;
						break;
					}	
				}
			}
			if (cf == 1)
				loc += dloc;			
		}
		if (cf == 1){
			perc = 1;
			percdir = "y";
		}
	}

	// check percolation in Z direction (post-process)		
	if (perc == 0){
		cf = 1;
		loc = 0;
		dloc = B[2]/subdiv;
		while(cf == 1 && loc < B[2]-dloc){
			cf = 0;

			// check for all points between loc & loc+dloc
			for (e=0; e<NE; e++){

				// get edge positions (with translated edges)
				ez0 = ez[e]-ecomz+0.5;
				if (ez0 < 0)
					ez0 += B[2];
				else if (ez0 > B[2])
					ez0 -= B[2];

				// check for ex0 in x bin
				if (ez0 > loc && ez0 <= loc+dloc){
					// check if edge is occupied with biggest cluster
					if (this->findroot(e)==this->get_pclus()){					
						cf = 1;
						break;
					}	
				}
			}
			if (cf == 1)
				loc += dloc;			
		}
		if (cf == 1){
			perc = 1;
			percdir = "z";
		}
	}

	// set percolation value
	this->set_perc(perc);
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
		xyzobj << "Type_" << id;
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

// Print global vertex info
void voroperc::print_global_vertices(){
	int i,j,dup;
	double dx,dy,dz,dr;

	cout << "Printing global vertex information: " << endl;
	cout << "NV = " << NV << endl << endl;
	for (i=0; i<NV; i++){
		cout << "vertex " << i << " :";
		cout << setw(16) << vx[i];
		cout << setw(16) << vy[i];
		cout << setw(16) << vz[i];
		cout << endl;
	}

	// check for duplicates
	dup = 0;
	for (i=0; i<NV; i++){
		for (j=i+1; j<NV; j++){
			// get x distance
			dx = vx[j]-vx[i];
			dx = dx - B[0]*round(dx/B[0]);

			// get x distance
			dy = vy[j]-vy[i];
			dy = dy - B[1]*round(dy/B[1]);

			// get z distance
			dz = vz[j]-vz[i];
			dz = dz - B[2]*round(dz/B[2]);

			// get total distance
			dr = sqrt(dx*dx + dy*dy + dz*dz);

			if (dr < 1e-8){
				cout << "duplicate found at i=" << i << "j=" << j;
				cout << "; dr = " << dr << endl;
				dup++;
			}
		}
	}

	if (dup == 0)
		cout << "No duplicates found in global vertex positions" << endl;
	cout << endl << endl;
}

// Print vertex neighbor info
void voroperc::print_vertex_neighbors(){
	int i,j,nv;

	cout << "Printing vertex neighbors" << endl;
	for (i=0; i<NV; i++){
		cout << "vertex i = " << i << " neighbors : ";
		nv = vneigh[i].size();
		for (j=0; j<nv; j++)
			cout << setw(8) << vneigh[i][j];
		cout << endl;
	}
	cout << endl << endl;
}

void voroperc::print_edge_neighbors(){
	int i,j,ne;

	for (i=0; i<NE; i++){
		cout << "edge i = " << i << " neighbors : ";
		ne = eneigh[i].size();
		for (j=0; j<ne; j++)
			cout << setw(8) << eneigh[i][j];
		cout << endl;
	}
	cout << endl << endl;
}





