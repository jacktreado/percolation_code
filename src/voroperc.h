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
	// std::vector< std::vector<int> >* facen;
	std::vector<int>* vface;
	std::vector<double>* vpx;
	std::vector<double>* vpy;
	std::vector<double>* vpz;

	// vertex positions in global coodinate system
	std::vector<double> vx;
	std::vector<double> vy;
	std::vector<double> vz;

	// edge positions in global coordinate system
	std::vector<double> ex;
	std::vector<double> ey;
	std::vector<double> ez;

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

	// setters
	void set_NV(int val) { NV = val; };
	void set_NE(int val) { NE = val; };
	void set_NVp(int id, int val) { if (NVp != nullptr) {NVp[id] = val;} else {cout << "cannot set, not nullptr\n"; throw;} };
	void set_NEp(int id, int val) { if (NEp != nullptr) {NEp[id] = val;} else {cout << "cannot set, not nullptr\n"; throw;} };
	void set_NFp(int id, int val) { if (NFp != nullptr) {NFp[id] = val;} else {cout << "cannot set, not nullptr\n"; throw;} };

	
	// ptr initialization
	void init_NVp(int val){ if (NVp == nullptr) {NVp = new int[val];} else {cout << "cannot init, not nullptr\n"; throw;} };
	void init_NEp(int val){ if (NEp == nullptr) {NEp = new int[val];} else {cout << "cannot init, not nullptr\n"; throw;} };
	void init_NFp(int val){ if (NFp == nullptr) {NFp = new int[val];} else {cout << "cannot init, not nullptr\n"; throw;} };

	// get voronoi vertex information
	void get_voro(int printit);
	void store_particle_vertices(int id, std::vector<double>& v);
	void store_face_neighbors(int id, std::vector<int>& neigh);
	void store_face_vertices(int id, std::vector<int>& f_vert);

	// setup lattice for percolation

	// get edge percolation

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
};




#endif