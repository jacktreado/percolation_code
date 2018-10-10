#ifndef VOIDCLUSTER_H
#define VOIDCLUSTER_H

#include <iostream>
#include "clustertree.h"
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

class voidcluster : public clustertree
{
protected:
	// scalars
	double g;					// grid spacing		
	double ac;					// critical radius/probe size
	double vctot;				// total void volume at percolation
	int NP;						// # of particles

	// array of  box lengths, particle positions, radii
	double* B;					// NDIM array of box lengths
	double** pos;				// particle positions
	double* rad;				// particle radii

	// cell info, for faster void updating
	int NCL;					// number of cells along one dimension
	int NCELLS;					// total number of cells
	int NCN;					// number of cell neighbors
	double* sp;					// array of cell dimensions
	int** cellneighbors;		// array of cell neighbors
	double** cellpos;			// cell positions
	int* pclabel;				// cell label for each particle
	int* sclabel;				// cell label for each lattice site	
	std::vector<int>* pcell;	// array of vectors: indices of atoms in cell m
	std::vector<int>* scell;	// array of vectors: indices of sites in cell m
public:
	voidcluster(int np, int l, int ndim, int nnn);
	voidcluster(int np, int l, int ndim, int nnn, int larr[]);
	voidcluster(int np, int l, int ndim, int nnn, int nc);
	voidcluster(string& pstr, int l, int ndim, int nnn, int nc);
	voidcluster(string& spherestr, int sph, int l, int ndim, int nnn, int nc);
	~voidcluster();

	// setters/getters
	void initialize_box(double b);
	void initialize_particles();
	void initialize_cells();
	void set_radius(double a) {for(int i=0; i<NP; i++) {rad[i] = a;} };
	double get_vctot(){return pow(g,this->get_NDIM())*(double)this->get_lattice_sum();};
	double get_B(){return B[0];};
	double get_NP(){return NP;};

	// particle overlap check
	int check_rad_overlap(int site);
	int check_probe_overlap(int site, double a);

	// set particle positions
	void rand_sphere_pos(int seed);

	// lattice initialization 
	void rand_sphere_void(double a);
	void rand_sphere_pperc(double a);
	void probe_void(double a);

	// percolation threshold finder
	void find_perc(double aH, double aL, double epsilon);
	void find_particle_perc(double aH, double aL, double aC, double epsilon);
	void find_probe_perc(double aH, double aL, double aC, double epsilon);

	// voronoi vertex percolation finder
	void find_voro_perc(double epsilon, double seed, double aH, double aL);

	// cell list info
	void add_pcell(int c, int p);
	void add_scell(int c, int s);
	void add_pcell_neighbors(int c, int p);
	void clear_cells();
	void get_cell_positions();
	void get_cell_neighbors();	
	int get_pcell(int p);
	int get_scell(int s);
	void label_cells();

	// file printing
	void fprint_stats(ofstream& obj, int seed, int header);
	void fprint_stats(ofstream& obj, int header);
	void fprint_s(ofstream& obj);

	// test/extra
	void print_lattice() {this->print_std_lattice();};
	void print_labels() {this->print_rl_lattice();};
	void print_lattice_positions();
	void print_void_xyz(ofstream &obj);
	void print_res_xyz(ofstream &obj);
	void print_pcell();
	void print_cellpos();
	void print_cell_xyz(ofstream &obj);
	void load_zhe_pos(string &zstr);
};


#endif