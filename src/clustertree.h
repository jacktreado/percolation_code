/*

	CLUSTERTREE class

	BY Jack Treado, 06/17/2018

*/


#ifndef CLUSTERTREE_H
#define CLUSTERTREE_H

#include <iostream>
#include <vector>
using namespace std;

class clustertree{
private:
	// integers
	int L; 						// number of sites per side
	long long int N;			// number of lattice sites total	
	int NDIM;					// number of dimensions	
	int NNN;					// number of nearest neighbors

	// pointers (dyn. mem. alloc.)			
	int** nn; 					// number of nearest neighbors
	int* ptr; 					// array of pointers to root site	
	int* lattice; 				// occupied sites (or not) lattice

	// perc info	
	long long int pclus;					// percolated cluster label
	long long int smax;						// max cluster size
	long long int smean;					// mean cluster size
	long long int fcalls;					// number of function calls
	long long int perc;						// whether or not cluster percolations
	long long int cnum;						// number of distinct clusters
	vector<long long int> s;				// cluster size list
public:
	// constructor/destructor
	clustertree(int l, int ndim, int nnn); 		// constructor, initializes all pointers dynamically
	~clustertree();
	void reset_lattice(long long int n);
	void initialize_nn(long long int i, int n){nn[i] = new int[n];};	

	// setters/getters (for inherited classes)
	void set_lattice(long long int site, int val) {lattice[site] = val;};
	void set_ptr(long long int site, int val) {ptr[site] = val;};
	void set_nn(long long int i, int j, int val){nn[i][j] = val;};
	void set_perc(int val) {perc = val;};
	void reset_sys();
	void reset_ptr();
	int get_L() {return L;};
	long long int get_N() {return N;};
	int get_NDIM() {return NDIM;};
	int get_NNN() {return NNN;};
	long long int get_pclus() {return pclus;};
	long long int get_smax() {return smax;};
	long long int get_smean() {return smean;};
	long long int get_fcalls() {return fcalls;};
	long long int get_perc() {return perc;};	
	long long int get_cnum() {return cnum;};
	long long int get_lattice_site(int site) {return lattice[site];};
	long long int get_sum_s();
	long long int get_lattice_sum();

	// nearest neighbor initialization
	void square_lattice_E();					// nearest-neighbor (NN) square lattice (square edges)
	void square_lattice_EV();					// next-to-nearest-neighbor (NNN) square lattice (sqaure edges, vertices)
	void cubic_lattice_F();						// NN cubic lattice (cube faces)
	void cubic_lattice_FE();					// NNN cubic lattice (cube faces, edges)
	void cubic_lattice_FEV();					// NNNN cubic lattice (cube faces, edges, vertices)
	void print_nn();

	// initialization
	void rand_site(double p, int seed);	

	// percolation methods
	long long int findroot(long long int i); 						// finds root of ptr
	long long int findroot(long long int i, int &kf);				// finds root, keeps track of function calls
	void merge_clusters();						// merge clusters
	void merge_clusters(vector<int>& NNNvec);
	void merge_clusters_edge_perc(vector<int>& NNNvec,vector<int> ev_0[], vector<double>& vx_0, 
		vector<double>& vy_0, vector<double>& vz_0, double B_0[]);
	int get_site_distance(int s1, int s2);
	double get_site_distance(int s1, int s2, vector<int> ev_0[], vector<double>& vx_0, vector<double>& vy_0, vector<double>& vz_0);
	void merge_boundary_pairs(vector< vector<int> >& boundpairs, long long int& big, long long int& bigr);
	int check_spanning(vector<int> ev[], vector<double>& vx, vector<double>& vy, vector<double>& vz, double B[]);
	int span_check(vector<int> ev[], vector<double>& vtmp, double B[]);

	void post_process();						// check percolation, get cluster stats
	void post_process_voro();
	int perc_search_XY();
	int perc_search_XZ();
	int perc_search_YZ();
	int last_frame_XY(long long int i);
	int last_frame_XZ(long long int i);
	int last_frame_YZ(long long int i);

	// print to console/file
	void print_cluster_stats(double p, int header);
	void print_cluster_stats(ofstream& fobj, double p, int header);	

	void print_cluster_sizes();
	void print_cluster_sizes(ofstream& fobj);

	void print_std_lattice();
	void print_rl_lattice();

	// extra/test		
	void merge_clusters_test();
	void print_all_vals();	
};

#endif