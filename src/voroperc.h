/*

	VOROPERC class

	BY Jack Treado, 09/07/2018

*/

#ifndef VOROPERC
#define VOROPERC

#include "voidcluster.h"
#include "voro++.hh"
#include <vector>
#include <fstream>
#include <string>

class voroperc : public voidcluster
{
private:
	// Voronoi network information
	int* NVp; 					// number of vertices per particle
	int* NEp;					// number of edges per particle
	int* NFp; 					// number of faces per particle
	int NV; 					// total number of vertices
	int NE; 					// total number of edges

	// particle-to-global mapping info
	std::vector<int>* facen;
	std::vector<int>** vface;
	std::vector<double>* vpx;
	std::vector<double>* vpy;
	std::vector<double>* vpz;
	std::vector<int*> vmap;
	std::vector< std::vector<int> > vn_map;
	std::vector<int>* vneigh;

	// vertex positions in global coodinate system
	std::vector<double> vx;
	std::vector<double> vy;
	std::vector<double> vz;
	std::vector<int>* vp; 		// list of particles associated with each vertex

	// edge positions in global coordinate system
	std::vector<double> ex;
	std::vector<double> ey;
	std::vector<double> ez;
	double ecomx,ecomy,ecomz; 	// edge center of mass!

	// edge neighbors, vertex map, list of particles
	std::vector<int>* eneigh;
	int** eij;
	std::vector<int>* ep;


	// point of closest approach in global coordinate system
	std::vector<double> cpax;
	std::vector<double> cpay;
	std::vector<double> cpaz;

	// ofstream object
	std::ofstream xyzobj;
public:
	// constructors & destructors
	voroperc(int np);
	~voroperc();

	// getters
	int get_NV(){ return NV; };
	int get_NE(){ return NE; };
	void update_NV(){ NV = vx.size(); }

	// setters
	void set_NV(int val) { NV = val; };
	void set_NE(int val) { NE = val; };
	void set_NVp(int id, int val) { if (NVp != nullptr) {NVp[id] = val;} else {cout << "cannot set, nullptr\n"; throw;} };
	void set_NEp(int id, int val) { if (NEp != nullptr) {NEp[id] = val;} else {cout << "cannot set, nullptr\n"; throw;} };
	void set_NFp(int id, int val) { if (NFp != nullptr) {NFp[id] = val;} else {cout << "cannot set, nullptr\n"; throw;} };
	
	// ptr initialization
	void init_NVp(int val) { if (NVp == nullptr) {NVp = new int[val];} else {cout << "cannot init, not nullptr\n"; throw;} };
	void init_NEp(int val) { if (NEp == nullptr) {NEp = new int[val];} else {cout << "cannot init, not nullptr\n"; throw;} };
	void init_NFp(int val) { if (NFp == nullptr) {NFp = new int[val];} else {cout << "cannot init, not nullptr\n"; throw;} };
	void init_vface(int id, int val) { if (vface[id] == nullptr) {vface[id] = new std::vector<int>[val];} else {cout << "cannot init, not nullptr\n"; throw;} };
	void init_vneigh() { if (vneigh == nullptr) {vneigh = new std::vector<int>[NV];} else {cout << "cannot init, not nullptr\n"; throw;} };
	void init_vp() { if (vp == nullptr) {vp = new std::vector<int>[NV];} else {cout << "cannot init, not nullptr\n"; throw;} };

	// get voronoi vertex information
	void get_voro();
	void store_particle_vertices(int id, std::vector<double>& v);
	void store_face_neighbors(int id, std::vector<int>& neigh);
	void store_face_vertices(int id, std::vector<int>& f_vert);

	// merge vertices from face information
	void merge_vertices();
	int check_multiple_faces(int i, int f);
	int get_neighboring_face(int i, int f, int fnp);
	int get_closest_face(int i, int f, int fnp);
	void calc_face_com(int i, int f, std::vector<double>& com);
	int check_redundancy(int i, int f, int fnp, int fnf);
	void combine_faces(int i, int f, int fnp, int fnf);
	int check_vertex_history(int i, int f, int vi);
	int get_vtrue(int i, int vi);
	void add_vertex_neighbors(int vtrue, int vi, int i, int f);
	void add_vertex_neighbors(int i, int f);

	// use info from vertex merge to create connectivity list (vneigh)
	void update_vertex_neighbors();
	void add_unique_to_vneigh(int i, int vtrue);
	void add_unique_to_vp(int i, int vtrue);

	// get edges
	void get_edges();
	void get_edge_positions(int vi, int vj);
	void get_edge_neighbors(int e, int vi, int vj);
	void get_ep(int e, int i, int j);
	int get_edge_v1(int e);
	int get_edge_v2(int e);
	void get_edge_vertices(int e, int& v1, int& v2);
	void get_cpa();
	void get_particle_positions(int e, std::vector<double>& p0, std::vector<double>& p1, std::vector<double>& p2);

	// voronoi network edge percolation
	void voro_edge_perc(double epsilon, int seed, double aH, double aL);
	void setup_edge_perc_lattice();
	void set_lattice_midpoint_intersect();
	void set_lattice_cpa_intersect();
	void check_voro_percolation(std::string& str);
	int check_voro_perc_x(std::string& percdir);
	int check_voro_perc_y(std::string& percdir);
	int check_voro_perc_z(std::string& percdir);

	// print information: lattice, network, xyz
	void open_xyz(std::string& str){
		xyzobj.open(str.c_str());
		if (!xyzobj.is_open()){
			std::cout << "ERROR: xyzobj could not open..." << std::endl;
			std::cout << "ERROR: file str = " << str << std::endl;
			throw;
		}
	}
	void print_vertices_xyz(int id, voro::voronoicell_neighbor& c);
	void print_edge_intersects_xyz(int e, int v0, int v1, vector<double>& p0, vector<double>& p1, vector<double>& p2);
	void print_perc_cluster_xyz();
	void print_face_vectors();
	void print_global_vertices();
	void print_vertex_neighbors();
	void print_edge_neighbors();
};




#endif