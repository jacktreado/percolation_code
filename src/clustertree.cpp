/*

	Methods implementation 
	for clustertree class

	BY Jack Treado, 06/17/2018

*/


#include "clustertree.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

#define DEBUG_OFF




/* 
==================================

	CONSTRUCTORS & DESTRUCTORS		 

================================== 
*/


/* 
	Main constructor: initialize box length, box dimension, 
	number of nearest neighbors. 

	ptr -> -1
	lattice -> -1
	nn -> -1

	smax -> 0
	smean -> 0
	fcalls -> -1
	perc -> -1
	cnum -> 0
	s -> empty vector, no initialization needed
*/

clustertree::clustertree(int l, int ndim, int nnn){
	/* SAVE PRIVATE VARIABLES */
	L = l;
	NDIM = ndim;
	N = pow(L,NDIM);
	NNN = nnn;

	long long int EMPTY = -1*N-1;
	long long int i;

	/* INITIAlIZE POINTERS DYNAMICALLY */
	nn = nullptr;
	ptr = nullptr;
	lattice = nullptr;

	cout << "initial pointers:" << endl;
	cout << "nn = " << nn << endl;
	cout << "lattice = " << lattice << endl;
	cout << "ptr = " << ptr << endl;

	// nearest neighbors: rank = 2
	nn = new int*[N];

	// lattice ptr array: rank = 1
	ptr = new int[N];

	// lattice array: rank = 1
	lattice = new int[N];

	/* LOOP OVER POINTS, SET VALUES TO 0 */
	for (i=0; i<N; i++){
		// initialize nn row to nullptr
		nn[i] = nullptr;

		// allocate NNN to row, initialize to 0
		nn[i] = new int[nnn];
		for (int j=0; j<nnn; j++)
			nn[i][j] = 0;

		// initialize lattice ptr to 0
		ptr[i] = -1;

		// initialize lattice ptr to 0
		lattice[i] = 0;
	}

	/* INITIALIZE INTEGER PRIVATE VARIABLES */
	pclus = -1;
	smax = 0;
	smean = 0;
	fcalls = -1;
	perc = -1;
	cnum = 0;
}

// Destructor: delete all dynamically allocated memory
clustertree::~clustertree(){
	cout << "destroying clustertree ptrs ..." << endl;
	delete [] ptr;
	delete [] lattice;
	for (int i=0; i<N; i++){
		delete [] nn[i];
	}
	delete [] nn;

	ptr = nullptr;
	lattice = nullptr;
	nn = nullptr;

	cout << "== clustertree POINTERS AFTER clustertree DESTRUCTOR CALL:" << endl;
	cout << "ptr = " << ptr << endl;
	cout << "lattice = " << lattice << endl;
	cout << "nn = " << nn << endl;
}

void clustertree::reset_lattice(long long int n){
	// set irrelevant private variables to -1
	L = -1;
	NNN = -1;

	// delete allocated memory with old value of N
	delete [] ptr;
	delete [] lattice;
	for (int i=0; i<N; i++){
		delete [] nn[i];
	}
	delete [] nn;

	ptr = nullptr;
	lattice = nullptr;
	nn = nullptr;

	// get new value of N
	N = n;

	// allocate new memory
	ptr = new int[N];
	lattice = new int[N];
	nn = new int*[N];

	int i;
	for (i=0; i<N; i++){
		// set ptr to -1 (no clusters), lattice to 0 (no populated edges)
		ptr[i] = -1;
		lattice[i] = 0;

		// dont initialize nn yet, needs to be done in voronoi vertex percolation part
	}
}


/* 
==================================

	  	 SETTERS/GETTERS		 

================================== 
*/

long long int clustertree::get_sum_s(){
	long long int c,ssum;
	ssum = 0;
	for(c=0; c<cnum; c++)
		ssum += s.at(c);

	return ssum;
}

long long int clustertree::get_lattice_sum(){
	long long int i,lsum;
	lsum = 0;
	for (i=0; i<N; i++)		
		lsum += lattice[i];

	return lsum;
}

void clustertree::reset_sys(){
	// reset list s
	s.clear();

	// reset stats
	smax = 0;
	smean = 0;
	cnum = 0;
	fcalls = 0;
	perc = 0;
	
}

void clustertree::reset_ptr(){
	for (int i=0; i<N; i++)
		ptr[i] = -1;
}



/* 
==================================

	  LATTICE INITIALIZATION		 

================================== 
*/


// set NN lists for a square lattice in NDIM = 2 
// with NN = 4 (edges of square centered on i)
void clustertree::square_lattice_E(){
	long long int i;
	for(i=0; i<N; i++){		
		nn[i][0] = (i+N-1) % N; 	// left neighbor (i-1)
		nn[i][1] = (i+N-L) % N;		// bottom neighbor (j-1)
		nn[i][2] = (i+1) % N; 		// right neighbor (i+1)
		nn[i][3] = (i+L) % N;		// top neighbor (j+1)
		if (i % L == 0)
			nn[i][0] = i+L-1;
		if ((i+1) % L == 0)
			nn[i][2] = i-L+1;
	}
}

// set NN lists for a square lattice in NDIM = 2 
// with NN = 8 (next-to-NN, edges & vertices of square centered on i)
void clustertree::square_lattice_EV(){
	long long int i;
	for (i=0; i<N; i++){
		// edges
		nn[i][0] = (i+N-1) % N; 	// left neighbor (i-1)
		nn[i][1] = (i+N-L) % N;		// bottom neighbor (j-1)
		nn[i][2] = (i+1) % N; 		// right neighbor (i+1)
		nn[i][3] = (i+L) % N;		// top neighbor (j+1)
		if (i % L == 0)
			nn[i][0] = i+L-1;
		if ((i+1) % L == 0)
			nn[i][2] = i-L+1;

		// vertices
		nn[i][4] = nn[i][1] - 1;	// bottom-left corner (i-1,j-1)
		nn[i][5] = nn[i][1] + 1;	// bottom-right corner (i+1,j-1)
		nn[i][6] = nn[i][3] - 1;	// top-left corner (i-1,j+1)
		nn[i][7] = nn[i][3] + 1;	// top-right corner (i+1,j+1)
		if (i % L == 0){
			nn[i][4] = i+L-1;
			nn[i][6] = i+L-1;
		}
		if ((i+1) % L == 0){
			nn[i][5] = i-L+1;
			nn[i][7] = i-L+1;
		}
	}
}


/* VARIATIONS ON CUBIC LATTICES: NDIM = 3 */

// cubic lattice, check faces of cube centered on i
void clustertree::cubic_lattice_F(){
	long long int i,zi;
	for (i=0; i<N; i++){
		zi = i/(L*L);

		// faces
		nn[i][0] = (i+N-1) % N; 			// left neighbor (i-1)
		nn[i][1] = i-L;						// bottom neighbor (j-1)
		nn[i][2] = (i+1) % N; 				// right neighbor (i+1)
		nn[i][3] = i+L;						// top neighbor (j+1)
		nn[i][4] = (i+N-L*L) % N;			// backward neighbor (k-1)
		nn[i][5] = (i+L*L) % N;				// forward neighbor (k+1)



		if (i % L == 0)
			nn[i][0] = i+L-1;
		if ((i+1) % L == 0)
			nn[i][2] = i-L+1;	
		if (i-zi*L*L < L)
			nn[i][1] = i-L+L*L;	
		if ((i+L)/(L*L) > zi)
			nn[i][3] = i-L*L+L;		
	}
}

// cubic lattice, check faces,edges of cube centered on i
void clustertree::cubic_lattice_FE(){
	long long int i,zi;
	for (i=0; i<N; i++){
		zi = i/(L*L);

		// faces
		nn[i][0] = (i+N-1) % N; 			// left neighbor (i-1)
		nn[i][1] = i-L;						// bottom neighbor (j-1)
		nn[i][2] = (i+1) % N; 				// right neighbor (i+1)
		nn[i][3] = i+L;						// top neighbor (j+1)
		nn[i][4] = (i+N-L*L) % N;			// backward neighbor (k-1)
		nn[i][5] = (i+L*L) % N;				// forward neighbor (k+1)		

		// y-direction bc
		if (i-zi*L*L < L)
			nn[i][1] = i-L+L*L;	
		if ((i+L)/(L*L) > zi)
			nn[i][3] = i-L*L+L;	

		// edges

		// * xy plane
		nn[i][6] = nn[i][1] - 1;	// bottom-left neightbor (i-1,j-1)
		nn[i][7] = nn[i][1] + 1;	// bottom-right neighbor (i+1,j-1)
		nn[i][8] = nn[i][3] - 1; 	// top-left neighbor (i-1,j+1)
		nn[i][9] = nn[i][3] + 1;	// top-right neighbor (i+1.j+1)

		// * xz plane
		nn[i][10] = (nn[i][4] - 1) % N;	// back-left neighbor (i-1,k-1)
		nn[i][11] = (nn[i][4] + 1) % N;	// back-right neighbor (i+1,k-1)
		nn[i][12] = (nn[i][5] - 1) % N;	// front-left neighbor (i-1,k+1)		
		nn[i][13] = (nn[i][5] + 1) % N;	// back-front neighbor (i+1,k+1)

		// * yz plane
		nn[i][14] = nn[i][4] - L;	// back-bottom neighbor (j-1,k-1)
		nn[i][15] = nn[i][4] + L;	// back-top neighbor (j+1,k-1)
		nn[i][16] = nn[i][5] - L;	// front-bottom neighbor (j-1,k+1)
		nn[i][17] = nn[i][5] + L; 	// front-top neighbor (j+1,k+1)

		// y-direction bc
		if (i-zi*L*L < L){
			nn[i][14] = nn[i][4]-L+L*L;
			nn[i][16] = nn[i][5]-L+L*L;	
		}
		if ((i+L)/(L*L) > zi){
			nn[i][15] = nn[i][4]-L*L+L;
			nn[i][17] = nn[i][5]-L*L+L;			
		}


		if (i % L == 0){
			// left BC
			nn[i][0] = i+L-1;
			nn[i][6] = nn[i][1]+L-1;
			nn[i][8] = nn[i][3]+L-1;
			nn[i][10] = nn[i][4]+L-1;
			nn[i][12] = nn[i][5]+L-1;			
		}
		if ((i+1) % L == 0){
			// right BC
			nn[i][2] = i-L+1;
			nn[i][7] = nn[i][1]-L+1;
			nn[i][9] = nn[i][3]-L+1;
			nn[i][11] = nn[i][4]-L+1;
			nn[i][13] = nn[i][5]-L+1;
		}		
	}
}


// cubic lattice, check faces,edges,vertices of cube centered on i
void clustertree::cubic_lattice_FEV(){
	long long int i,zi;
	for (i=0; i<N; i++){
		zi = i/(L*L);		
		
		// faces
		nn[i][0] = (i+N-1) % N; 			// left neighbor (i-1)
		nn[i][1] = i-L;						// bottom neighbor (j-1)
		nn[i][2] = (i+1) % N; 				// right neighbor (i+1)
		nn[i][3] = i+L;						// top neighbor (j+1)
		nn[i][4] = (i+N-L*L) % N;			// backward neighbor (k-1)
		nn[i][5] = (i+L*L) % N;				// forward neighbor (k+1)		

		// y-direction bc
		if (i-zi*L*L < L)
			nn[i][1] = i-L+L*L;	
		if ((i+L)/(L*L) > zi)
			nn[i][3] = i-L*L+L;	

		// edges

		// * xy plane
		nn[i][6] = nn[i][1] - 1;	// bottom-left neightbor (i-1,j-1)
		nn[i][7] = nn[i][1] + 1;	// bottom-right neighbor (i+1,j-1)
		nn[i][8] = nn[i][3] - 1; 	// top-left neighbor (i-1,j+1)
		nn[i][9] = nn[i][3] + 1;	// top-right neighbor (i+1,j+1)

		// * xz plane
		nn[i][10] = (nn[i][4] - 1) % N;	// back-left neighbor (i-1,k-1)
		nn[i][11] = (nn[i][4] + 1) % N;	// back-right neighbor (i+1,k-1)
		nn[i][12] = (nn[i][5] - 1) % N;	// front-left neighbor (i-1,k+1)		
		nn[i][13] = (nn[i][5] + 1) % N;	// back-front neighbor (i+1,k+1)

		// * yz plane
		nn[i][14] = nn[i][4] - L;	// back-bottom neighbor (j-1,k-1)
		nn[i][15] = nn[i][4] + L;	// back-top neighbor (j+1,k-1)
		nn[i][16] = nn[i][5] - L;	// front-bottom neighbor (j-1,k+1)
		nn[i][17] = nn[i][5] + L; 	// front-top neighbor (j+1,k+1)

		// y-direction bc
		if (i-zi*L*L < L){
			nn[i][14] = nn[i][4]-L+L*L;
			nn[i][16] = nn[i][5]-L+L*L;	
		}
		if ((i+L)/(L*L) > zi){
			nn[i][15] = nn[i][4]-L*L+L;
			nn[i][17] = nn[i][5]-L*L+L;			
		}
		


		// cubic vertices
		nn[i][18] = (nn[i][14] - 1) % N; // back-bottom-left neighbor (i-1,j-1,k-1)
		nn[i][19] = (nn[i][14] + 1) % N; // back-bottom-right neighbor (i+1,j-1,k-1)
		nn[i][20] = (nn[i][15] - 1) % N; // back-top-left neighbor (i-1,j+1,k-1)
		nn[i][21] = (nn[i][15] + 1) % N; // back-top-right neighbor (i+1,j+1,k-1)
		nn[i][22] = (nn[i][16] - 1) % N; // front-bottom-left neighbor (i-1,j-1,k+1)
		nn[i][23] = (nn[i][16] + 1) % N; // front-bottom-right neighbor (i+1,j-1,k+1)
		nn[i][24] = (nn[i][17] - 1) % N; // front-top-left neighbor (i-1,j+1,k+1)
		nn[i][25] = (nn[i][17] + 1) % N; // front-top-right neighbor (i+1,j+1,k+1)



		if (i % L == 0){
			// left BC
			nn[i][0] = i+L-1;
			nn[i][6] = nn[i][1]+L-1;
			nn[i][8] = nn[i][3]+L-1;
			nn[i][10] = nn[i][4]+L-1;
			nn[i][12] = nn[i][5]+L-1;	
			nn[i][18] = nn[i][14]+L-1;
			nn[i][20] = nn[i][15]+L-1;
			nn[i][22] = nn[i][16]+L-1;
			nn[i][24] = nn[i][17]+L-1;
		}
		if ((i+1) % L == 0){
			// right BC
			nn[i][2] = i-L+1;
			nn[i][7] = nn[i][1]-L+1;
			nn[i][9] = nn[i][3]-L+1;
			nn[i][11] = nn[i][4]-L+1;
			nn[i][13] = nn[i][5]-L+1;
			nn[i][19] = nn[i][14]-L+1;
			nn[i][21] = nn[i][15]-L+1;
			nn[i][23] = nn[i][16]-L+1;
			nn[i][25] = nn[i][17]-L+1;
		}
	}	
}

void clustertree::print_nn(){
	long long int i,n;
	cout << "Printing nearest neighbors..." << endl;
	for (i=0; i<N; i++){
		cout << "i = " << i << "; ";
		for (n=0; n<NNN; n++)
			cout << setw(5) << nn[i][n];
		cout << endl;
	}
	cout << endl << endl;

}


// Randomly populate lattice as a function of p
void clustertree::rand_site(double p, int seed){
	this->reset_sys();

	srand48(seed);
	long long int i;
	double r;
	for (i=0; i<N; i++){
		r = drand48();
		if (r <= p)
			lattice[i] = 1;					
		else
			lattice[i] = 0;
		ptr[i] = -1;
	}
}








/*
============================

	PERCOLATION METHODS

============================
*/


// get root site of given site i
long long int clustertree::findroot(long long int i){
	// cout << "&&& in findroot(), i = " << i << endl;
	if (ptr[i]<0)
		return i;
	else
		return ptr[i] = findroot(ptr[i]);
}


// get root site of given site i, keep track of function evaluations
long long int clustertree::findroot(long long int i, int &kf){
	kf++;
	if (ptr[i]<0)
		return i;
	else
		return ptr[i] = findroot(ptr[i],kf);
}


/*
	UNION/FIND ALGORITHM:
		1. Find root of site i
		2. Check all populated lattice sites adjacent to i
		3. For all populated lattice sites r2, get root of r2
		4. IF r2 != r1, then clusters are merged using weighted
			merging.
		5. Largest cluster size is monitored.	
*/
void clustertree::merge_clusters(){
	long long int i,j;
	long long int s1,s2;
	long long int r1,r2;
	long long int big = 0;
	long long int bigr = 0;
	long long int nn_tmp;

	// # of function calls
	int kf = 0;

	for (i=0; i<N; i++){				
		// only continue if on an occupied lattice site
		if (lattice[i] == 1){
			s1 = i;
			r1 = this->findroot(s1,kf);

			for (j=0; j<NNN; j++){
				s2 = nn[i][j];

				// skip nn if lattice is empty
				if (lattice[s2] == 1){
					// get root of adj pt
					r2 = this->findroot(s2,kf);

					// if diff roots, then need to merge (if same, already merged!)
					if (r2 != r1){

						// if r2 > r1 (mind minus sign), merge 1 -> 2
						if (ptr[r1] > ptr[r2]){
							ptr[r2] += ptr[r1];
							ptr[r1] = r2;
							r1 = r2;
						}

						// else, merge r2 -> r1
						else{
							ptr[r1] += ptr[r2];
							ptr[r2] = r1;
						}

						// if new cluster is max, increase max!
						if (-ptr[r1] > big){							
							big = -ptr[r1];
							bigr = r1;
						}
					}
				}
			}
		}
	}

	smax = big;
	pclus = bigr;
	fcalls = kf;
}

void clustertree::merge_clusters(vector<int>& NNNvec){
	long long int i,j;
	long long int s1,s2;
	long long int r1,r2;
	long long int big = 0;
	long long int bigr = 0;
	long long int nn_tmp;

	// # of function calls
	int kf = 0;

	for (i=0; i<N; i++){				
		// only continue if on an occupied lattice site
		if (lattice[i] == 1){
			s1 = i;
			r1 = this->findroot(s1,kf);

			for (j=0; j<NNNvec.at(i); j++){
				s2 = nn[i][j];

				// skip nn if lattice is empty
				if (lattice[s2] == 1){
					// get root of adj pt
					r2 = this->findroot(s2,kf);

					// if diff roots, then need to merge (if same, already merged!)
					if (r2 != r1){

						// if r2 > r1 (mind minus sign), merge 1 -> 2
						if (ptr[r1] > ptr[r2]){
							ptr[r2] += ptr[r1];
							ptr[r1] = r2;
							r1 = r2;
						}

						// else, merge r2 -> r1
						else{
							ptr[r1] += ptr[r2];
							ptr[r2] = r1;
						}

						// if new cluster is max, increase max!
						if (-ptr[r1] > big){							
							big = -ptr[r1];
							bigr = r1;
						}
					}
				}
			}
		}
	}

	smax = big;
	pclus = bigr;
	fcalls = kf;
}


void clustertree::merge_clusters_edge_perc(vector<int>& NNNvec, vector<int> ev_0[], vector<double>& vx_0, vector<double>& vy_0, vector<double>& vz_0, double B_0[]){
	long long int ii,i,j,irand1,irand2,itmp;
	long long int s1,s2;
	long long int r1,r2;
	long long int big = 0;
	long long int bigr = 0;
	long long int nn_tmp;
	int span;
	// long long int v1,v2,vr;
	// double dr1x,dr1y,dr1z,dr1;
	// double dr2x,dr2y,dr2z,dr2;
	// double dx1,dy1,dz1,dv;
	// double dx2,dy2,dz2,d;
	double dx,dy,dz,ds;

	// vector of possible cross-boundary pairs
	vector< vector<int> > boundpairs;
	vector<int> vtmp(4);

	// # of function calls
	int kf = 0;

	// percolated or not
	perc = 0;

	// create random list of indices
	vector<long long int> uq_inds(N);
	for (i=0; i<N; i++){
		uq_inds[i] = i;
	}

	
	int NRM = 1e3;
	for (i=0; i<NRM; i++){
		irand1 = round((N-1)*drand48());
		irand2 = round((N-1)*drand48());
		itmp = uq_inds[irand1];
		uq_inds[irand1] = uq_inds[irand2];
		uq_inds[irand2] = itmp;
	}
	

	for (ii=0; ii<N; ii++){
		// get random index
		i = uq_inds[ii];

		// only continue if on an occupied lattice site
		if (lattice[i] == 1){
			s1 = i;
			r1 = this->findroot(s1,kf);

			for (j=0; j<NNNvec.at(i); j++){
				s2 = nn[i][j];				

				// merge if lattice site is occupied
				if (lattice[s2] == 1){
					// get distance between 2 sites (specifically for edge percolation)
					dx = 0;
					dy = 0;
					dz = 0;
					this->get_site_distance(s1,s2,ev_0,vx_0,vy_0,vz_0,dx,dy,dz);

					// get root of adj pt
					r2 = this->findroot(s2,kf);

					// if ds > (1/2) box length, then boundary points, so do not merge yet
					if (abs(dx) > 0.5*B_0[0]){
						vtmp[0] = s1;
						vtmp[1] = s2;
						vtmp[2] = 0;
						vtmp[3] = round(dx);		
						boundpairs.push_back(vtmp);
						continue;
					}
					else if (abs(dy) > 0.5*B_0[1]){
						vtmp[0] = s1;
						vtmp[1] = s2;
						vtmp[2] = 1;
						vtmp[3] = round(dy);			
						boundpairs.push_back(vtmp);
						continue;
					}
					else if (abs(dz) > 0.5*B_0[2]){
						vtmp[0] = s1;
						vtmp[1] = s2;
						vtmp[2] = 2;
						vtmp[3] = round(dz);			
						boundpairs.push_back(vtmp);
						continue;
					}


					// if diff roots, then need to merge (if same, already merged!)
					if (r2 != r1){
						// if r2 > r1 (mind minus sign), merge 1 -> 2
						if (ptr[r1] > ptr[r2]){
							ptr[r2] += ptr[r1];
							ptr[r1] = r2;
							r1 = r2;
						}

						// else, merge r2 -> r1
						else{
							ptr[r1] += ptr[r2];
							ptr[r2] = r1;
						}

						// if new cluster is max, increase max!
						if (-ptr[r1] > big){							
							big = -ptr[r1];
							bigr = r1;
						}
					}
				}				
			}
		}
	}

	// detect spanning cluster
	pclus = bigr;
	span = this->check_spanning(boundpairs);
	// cout << "span = " << span << endl;

	// if spanning cluster found, check boundary pairs
	if (span == 1){
		// perc = 1;
		this->merge_boundary_pairs(boundpairs,big,bigr);
		cout << endl;
	}	

	smax = big;
	pclus = bigr;
	fcalls = kf;
}


void clustertree::get_site_distance(int s1, int s2, vector<int> ev[], vector<double>& vx, vector<double>& vy, vector<double>& vz, double& dx, double& dy, double& dz){
	// get non-pbc distance between edge s1 and s2
	int v1,v2;

	// get principle vertices for s1 & s2
	v1 = ev[s1][0];
	v2 = ev[s2][0];

	// get distances
	dx = vx[v2]-vx[v1];
	dy = vy[v2]-vy[v1];
	dz = vz[v2]-vz[v1];
}

int clustertree::check_spanning(vector< vector<int> >& boundpairs){
	int i,s1,s2,r1,r2,ds,ptrGS1,ptrGS2,rgs1,rgs2,sz,pclussz,d,s,nbp,span;
	int pcf = -1;
	int new_pclus;
	span = 0;

	// span = this->span_check(ev,vx,B);
	// if (span == 0)
	// 	span = this->span_check(ev,vy,B);
	// if (span == 0)
	// 	span = this->span_check(ev,vz,B);

	// merge with boundary ghost points
	int ghostp[NDIM][2];
	for (d=0; d<NDIM; d++){
		ghostp[d][0] = -1;
		ghostp[d][1] = -1;
	}

	// loop over boundary points, decide which boundary each is on, merge ghost point 
	// to largest cluster on given boundary
	nbp = boundpairs.size();
	for (i=0; i<nbp; i++){
		// get sites, crossing boundary
		s1 = boundpairs[i][0];
		s2 = boundpairs[i][1];
		d = boundpairs[i][2];
		ds = boundpairs[i][3];

		r1 = this->findroot(s1);
		r2 = this->findroot(s2);
		// cout << "s1 = " << s1 << ", s2 = " << s2 << ", d = " << d << ", ds = " << ds << ", r1 = " << r1 << ", r2 = " << r2 << endl;

		// get boundary side s of s1: s = 0 -> 0, s = 1 -> L
		if (ds > 0)
			s = 0;
		else if (ds < 0)
			s = 1;
		else{
			cout << "ds = 0 in spanning check, throwing error..." << endl;
			throw;
		}


		// get current root site of gs on side of s1
		if (ghostp[d][s]>=0)
			ptrGS1 = ghostp[d][s];
		else{
			ptrGS1 = s1;
			ghostp[d][s] = ptrGS1;
		}

		// get current root site of gs on side of s2		
		if (ghostp[d][1-s]>=0)
			ptrGS2 = ghostp[d][1-s];
		else{
			ptrGS2 = s2;
			ghostp[d][1-s] = ptrGS2;
		}

		// get root sites of s1, s2, and ghost points		
		rgs1 = this->findroot(ptrGS1);
		rgs2 = this->findroot(ptrGS2);

		// point ghostp[d][s] to larger cluster if cluster with root site r1 is larger
		// than cluster 
		if (ptr[r1] < ptr[rgs1])
			ghostp[d][s] = s1;
		// else, either ghostp points to s1, OR ghostp points to larger cluster. In either
		// case, do nothing

		// check same for s2
		if (ptr[r2] < ptr[rgs2])
			ghostp[d][1-s] = s2;

	}

	// once ghost points are merged with largest cluster on each boundary, check across boundary
	// to see if pair ghost point has same root cluster. IF SO, then spanning cluster found!
	for (d=0; d<NDIM; d++){
		if (ghostp[d][0] >= 0 && ghostp[d][1] >= 0){
			rgs1 = this->findroot(ghostp[d][0]);
			rgs2 = this->findroot(ghostp[d][1]);
			if (rgs1 == rgs2){
				if (span == 0){
					cout << "spanning cluster detected! pclus = " << pclus;
					span = 1;
				}				

				// check if cluster > pclus
				sz = -ptr[rgs1];
				pclussz = -ptr[pclus];
				if (sz < pclussz){
					cout << "; d = " << d << ", rgs1 = " << rgs1 << ", spanning cluster smaller than pclus found...";
					if (pcf == -1){
						pcf = 0;
						new_pclus = rgs1;
					}
				}
				else if (sz == pclussz){
					cout << "; d = " << d << ", rgs1 = " << rgs1 << ", which is the pclus! ";
					pcf = 1;
				}
				else{
					cout << "; d = " << d << ", rgs1 = " << rgs1 << ", spanning cluster found that is larger than pclus, so making pclus spanning cluster...";					
					pclus = rgs1;
				}
			}
		}
	}

	if (pcf == 0){
		cout << "pclus not spanning, so setting pclus to spanning cluster! ";
		pclus = new_pclus;
	}
	
	// return value of span
	return span;
}

int clustertree::span_check(vector<int> ev[], vector<double>& vtmp, double B[]){
	// local variables
	int e,span,boxedge,cf,v1;
	double loc,loc0,loc1,dloc,subdiv,p;

	// check percolation in X direction
	span = 0;
	cf = 1;
	subdiv = 10.0;
	loc0 = -2*B[0];
	loc1 = 0;
	loc = 0;
	dloc = B[0]/subdiv;
	boxedge = 0;

	while(cf == 1 && boxedge == 0){
		cf = 0;

		// check which bin to test	
		// if loc is 0, check left edge and out to -inf
		if (abs(loc) < 1e-14){	
			loc0 = -10*B[0];
			loc1 = loc+dloc;
		}
		// if loc is close to box edge, check from edge to +inf
		else if (loc >= B[0]-dloc-1e-8){
			loc0 = loc;
			loc1 = loc+10*B[0];
			boxedge = 1;
		}
		else{
			loc0 = loc;
			loc1 = loc+dloc;
		}		

		// check for all points between loc & loc+dloc
		for (e=0; e<N; e++){
			if (lattice[e]==0)
				continue;
					
			v1 = ev[e][0];
			p = vtmp[v1];

			// check for v pos in bin
			if (p > loc0 && p <= loc1){
				// check if edge is occupied with biggest cluster && 
				if (this->findroot(e)==pclus){					
					cf = 1;
					break;
				}	
			}
		}
		if (cf == 1)
			loc += dloc;

	}
	if (cf == 1 && boxedge == 1)
		span = 1;

	return span;
}


void clustertree::merge_boundary_pairs(vector< vector<int> >& boundpairs, long long int& big, long long int& bigr){
	// THIS FORM GIVES CORRECT nu VALUE! UPDATE pclus IN SECOND LOOP!!

	int i,nbp,s1,s2,r1,r2;

	nbp = boundpairs.size();
	// unite non-spanning pairs
	for (i=0; i<nbp; i++){
		s1 = boundpairs[i][0];
		s2 = boundpairs[i][1];

		r1 = this->findroot(s1);
		r2 = this->findroot(s2);	
		
		if (r1 != r2 && r1 != pclus && r2 != pclus){
			// if r2 > r1 (mind minus sign), merge 1 -> 2
			if (ptr[r1] > ptr[r2]){
				ptr[r2] += ptr[r1];
				ptr[r1] = r2;
				r1 = r2;
			}

			// else, merge r2 -> r1
			else{
				ptr[r1] += ptr[r2];
				ptr[r2] = r1;
			}

			// if new cluster is max, increase max!
			if (-ptr[r1] > big){							
				big = -ptr[r1];
				bigr = r1;
			}
		}
	}
	if (-ptr[pclus] < big){		
		cout << "; nonspanning found that is larger than old pclus!";
	}

	// unite non-spanning with spanning
	for (i=0; i<nbp; i++){
		s1 = boundpairs[i][0];
		s2 = boundpairs[i][1];

		r1 = this->findroot(s1);
		r2 = this->findroot(s2);	
		
		if (r1 != r2 && ( (r1 != pclus && r2 == pclus) || (r1 == pclus && r2 != pclus) ) ) {
			// if r2 > r1 (mind minus sign), merge 1 -> 2
			if (ptr[r1] > ptr[r2]){
				ptr[r2] += ptr[r1];
				ptr[r1] = r2;
				r1 = r2;
			}

			// else, merge r2 -> r1
			else{
				ptr[r1] += ptr[r2];
				ptr[r2] = r1;
			}

			// if new cluster is max, increase max!
			if (-ptr[r1] > big){							
				big = -ptr[r1];
				bigr = r1;
			}
			pclus = bigr;
		}		
	}	
	if (-ptr[pclus] < big && bigr != pclus){		
		cout << "; pclus merged into larger cluster!";
		pclus = bigr;
	}
	else if (-ptr[pclus] < big && bigr == pclus)
		cout << "; pclus still the largest cluster...";

	// check if 
	for (i=0; i<nbp; i++){
		s1 = boundpairs[i][0];
		s2 = boundpairs[i][1];

		r1 = this->findroot(s1);
		r2 = this->findroot(s2);	
		
		if (r1 == r2 && r1 == pclus){
			perc = 1;
		}
	}
}



















// perc search - XY
int clustertree::perc_search_XY(){
	long long int i,j;
	long long int loc,slice,tmp;
	int percXY,cf;

	slice = pow(L,NDIM-1);

	// check XY slices for largest cluster, if found then advance to next slice
	loc = 0;
	cf = 1;

	while(cf==1 && loc < N){
		// set cf = 0			
		cf = 0;	

		// check NDIM-1 slice to see if cmax is there
		for (i=0; i<slice; i++){
			// map to location in array
			j = loc+i;

			// find root of tmp location, break if cmax found
			tmp = findroot(j);
			if (tmp == pclus){
				// cout << "clus found at i = " << i << endl;
				if (NNN == 26){
					cf = this->last_frame_XY(j);
					if (cf == 1)
						break;
					else
						cout << "does not connect to last frame!\n";
				}
				else{
					cf = 1;
					break;
				}
			}
		}

		// if cf still = 1, increase loc, check next slice
		if (cf == 1){
			loc += slice;
		}
	}

	if (cf == 1)
		percXY = 1;
	else
		percXY = 0;

	return percXY;

}

// check that frame i connects to frame i-1
int clustertree::last_frame_XY(long long int i){
	int* lflist;
	int NLFN = 9;
	int j,ind,neighbor;

	lflist = new int[NLFN];
	lflist[0] = nn[i][4];
	lflist[1] = nn[i][10];
	lflist[2] = nn[i][11];
	lflist[3] = nn[i][14];
	lflist[4] = nn[i][15];
	lflist[5] = nn[i][18];
	lflist[6] = nn[i][19];
	lflist[7] = nn[i][20];
	lflist[8] = nn[i][21];

	for (j=0; j<NLFN; j++){
		ind = lflist[j];
		neighbor = findroot(ind);
		if (neighbor == pclus){
			delete [] lflist;
			return 1;
		}
	}

	delete [] lflist;
	return 0;
}




// perc search - XY
int clustertree::perc_search_XZ(){
	long long int i,k,xx,zz;
	long long int orig,tmp;
	int percXZ,cf;

	// keeps track of which slice
	k = 0;

	// check XY slices for largest cluster, if found then advance to next slice
	cf = 1;

	while(k < L){
		// set cf = 0			
		cf = 0;	

		// set origin
		orig = L*k;

		// loop over XZ slice, incrementing by L^2 ever L sites
		for (zz=0; zz<L; zz++){				
			for (xx=0; xx<L; xx++){	
				// increment site		
				i = orig+xx;

				// find root of site index
				tmp = findroot(i);
				if (tmp == pclus){
					// cout << "clus found at i = " << i << endl;
					if (NNN == 26){
						cf = this->last_frame_XZ(i);
						if (cf == 1)
							break;
						else
							cout << "does not connect to last frame!\n";
					}
					else{
						cf = 1;
						break;
					}
				}
			}
			// break zz loop if max clus found
			if (cf == 1)
				break;
			else
				// increment to next x-row in z
				orig += L*L;
		}

		if (cf == 0)
			break;

		// increment y direction
		k++;
	}

	// if cf = 1 at the end of while loop, then max clus found in each slice!
	if (cf == 1)
		percXZ = 1;
	// else, gap in max cluster and no percolation
	else
		percXZ = 0;		

	return percXZ;
}

// check that frame i connects to frame i-1
int clustertree::last_frame_XZ(long long int i){
	int* lflist;
	int NLFN = 9;
	int j,ind,neighbor;

	lflist = new int[NLFN];
	lflist[0] = nn[i][1];
	lflist[1] = nn[i][6];
	lflist[2] = nn[i][7];
	lflist[3] = nn[i][14];
	lflist[4] = nn[i][16];
	lflist[5] = nn[i][18];
	lflist[6] = nn[i][19];
	lflist[7] = nn[i][22];
	lflist[8] = nn[i][23];

	for (j=0; j<NLFN; j++){
		ind = lflist[j];
		neighbor = findroot(ind);
		if (neighbor == pclus){
			delete [] lflist;
			return 1;
		}
	}

	return 0;
}





// perc search - XY
int clustertree::perc_search_YZ(){
	long long int i,k,yy,zz;
	long long int orig,tmp;
	int percYZ,cf;

	// keeps track of which slice
	k = 0;

	// check YZ slices for largest cluster, if found then advance to next slice
	cf = 1;

	while(k < L){
		// set cf = 0			
		cf = 0;	

		// set origin
		orig = k;

		// loop over XZ slice, incrementing by L^2 ever L sites
		for (zz=0; zz<L; zz++){				
			for (yy=0; yy<L; yy++){	
				// increment site		
				i = orig+L*yy;				

				// find root of site index
				tmp = findroot(i);				
				if (tmp == pclus){
					// cout << "clus found at i = " << i << endl;
					if (NNN == 26){
						cf = this->last_frame_YZ(i);
						if (cf == 1)
							break;
						else
							cout << "does not connect to last frame!\n";
					}
					else{
						cf = 1;
						break;
					}
				}
			}
			// break zz loop if max clus found
			if (cf >= 1)
				break;
			// else, increment to next y-row in z direction
			else					
				orig += L*L;
		}

		if (cf == 0)
			break;

		// increment y direction
		k++;
	}

	// if cf = 1 at the end of while loop, then max clus found in each slice!
	if (cf == 1)
		percYZ = 1;
	// else, gap in max cluster and no percolation
	else
		percYZ = 0;


	return percYZ;

}

// check that frame i connects to frame i-1
int clustertree::last_frame_YZ(long long int i){
	int* lflist;
	int NLFN = 9;
	int j,ind,neighbor;

	lflist = new int[NLFN];
	lflist[0] = nn[i][0];
	lflist[1] = nn[i][6];
	lflist[2] = nn[i][8];
	lflist[3] = nn[i][10];
	lflist[4] = nn[i][12];
	lflist[5] = nn[i][18];
	lflist[6] = nn[i][20];
	lflist[7] = nn[i][22];
	lflist[8] = nn[i][24];

	for (j=0; j<NLFN; j++){
		ind = lflist[j];
		neighbor = findroot(ind);
		if (neighbor == pclus){
			delete [] lflist;
			return 1;
		}
	}

	delete [] lflist;
	return 0;
}

void clustertree::post_process_voro(){
	// get number of root clusters at least > L in 

	long long int i;
	long long int tmp;
	int nocs = 0;

	// set cnum = 0, smax = 0;
	cnum = 0;
	for (i=0; i<N; i++){
		// check if lattice is occupied or not
		if (lattice[i]==1){
			// increase number of occupied sites
			nocs++;
			tmp = ptr[i];			
			if (tmp < 0){
				// increase number of clusters
				cnum++;

				// save cluster to vector s
				s.push_back(-1*tmp);
			}
		}		
	}

	// set mean cluster size, private variable smean
	if (cnum > 0)
		smean = round((double)nocs/(double)cnum);
	else{
		cout << "cnum == 0, no clusters found...." << endl;
		smean = 0;
		smax = 0;
	}
}

void clustertree::post_process(){
	// get number of root clusters at least > L in 

	long long int i,j;
	long long int tmp,max;
	int nocs = 0;

	max = L-1;

	// set cnum = 0, smax = 0;
	cnum = 0;
	for (i=0; i<N; i++){
		// check if lattice is occupied or not
		if (lattice[i]==1){
			// increase number of occupied sites
			nocs++;
			tmp = ptr[i];			
			if (tmp < 0){
				// increase number of clusters
				cnum++;

				// save cluster to vector s
				s.push_back(-1*tmp);
			}
		}		
	}

	// set mean cluster size, private variable smean
	if (cnum > 0)
		smean = round((double)nocs/(double)cnum);
	else{
		// cout << "cnum == 0, no clusters found...." << endl;
		smean = 0;
		smax = 0;
	}

	// check for percolation
	if (smax > max){		
		perc = this->perc_search_XY();
		if (perc == 1)
			cout << "XY";

		if (NDIM == 3){
			if (perc == 0){			
				// cout << "checking XZ perc....";
				perc = this->perc_search_XZ();	
				if (perc == 1)
					cout << "XZ";		

				if (perc == 0){
					// cout << "checking YZ perc....";
					perc = this->perc_search_YZ();
					if (perc == 1)
						cout << "YZ";
				}
			}
		}
		else if(NDIM > 3){
			cout << "NDIM > 3 not supported at this time." << endl;
			throw "NDIM > 3 not yet supported.\n";
		}
	}
	else		
		perc = 0;

	if (perc == 0)
		pclus = -1;
}


/*
============================

	  PRINTER METHODS

============================
*/

// print stats to console
void clustertree::print_cluster_stats(double p, int header){
	// Print header
	int w = 12;
	if (header == 1){
		// Print header
		cout << setw(w) << "p";
		cout << setw(w) << "#clusters";
		cout << setw(w) << "perc";
		cout << setw(w) << "smax";
		cout << setw(w) << "smean";
		cout << setw(w) << "fcalls";
		cout << endl;
		for (int i=0; i<6*w; i++)
			cout << "=";
		cout << endl;
	}

	cout << setw(w) << setprecision(6) << p;
	cout << setw(w) << cnum;
	cout << setw(w) << perc;
	cout << setw(w) << smax;
	cout << setw(w) << smean;
	cout << setw(w) << fcalls;
	cout << endl; 
}


// print stats to file
void clustertree::print_cluster_stats(ofstream& fobj, double p, int header){
	cout << "Printing cluster stats to ofstream object." << endl;

	int w = 12;	
	if (header == 1){
		// Print header
		fobj << setw(w) << "p";
		fobj << setw(w) << "#clusters";
		fobj << setw(w) << "perc";
		fobj << setw(w) << "smax";
		fobj << setw(w) << "smean";
		fobj << setw(w) << "fcalls";
		fobj << endl;
		for (int i=0; i<6*w; i++)
			fobj << "=";
		fobj << endl;
	}

	fobj << setw(w) << setprecision(6) << p;
	fobj << setw(w) << cnum;
	fobj << setw(w) << perc;
	fobj << setw(w) << smax;
	fobj << setw(w) << smean;
	fobj << setw(w) << fcalls;
	fobj << endl; 
}


// print cluster size distribution to console
void clustertree::print_cluster_sizes(){
	cout << "Printing cluster sizes to console: " << endl;
	long long int i;

	for (i=0; i<cnum; i++){
		cout << s.at(i) << endl;
	}

	cout << "cnum = " << cnum << ", s size = " << s.size() << endl;
}


// print cluster size distribution to file
void clustertree::print_cluster_sizes(ofstream &fobj){
	cout << "Printing cluster sizes to ofstream object." << endl;
	int i;
	for (i=0; i<cnum; i++) {fobj << s.at(i) << endl;}
}



// print lattice to console, unlabelled
void clustertree::print_std_lattice(){
	long long int i;
	for (i = 0; i < N; i++){
		if (i % L == 0)
			cout << endl << endl;
		if (i % (L*L) == 0)
			cout << endl << endl << "===" << endl << endl;
		if (lattice[i] == 1)
			cout << setw(5) << "*";
		else
			cout << setw(5) << "-";			
	}
	cout << endl << endl;
}


// print relabelled lattice to console
void clustertree::print_rl_lattice(){
	long long int i,r,d,test;
	const int SUB = NDIM - 1;
	for (i = 0; i < N; i++){

		// print lines for different higher dimensional 2D slices (# of lines = NDIM-2)
		d = 0;
		while (d < SUB){
			d++;
			test = pow(L,d);
			if (i % test == 0){
				cout << endl << endl;
				for(int dd = 1; dd < d; dd++){
					cout << "------------------------" << endl;
					cout << endl << endl;					
				}
			}
		}

		if (lattice[i] == 1){
			r = findroot(i);
			cout << setw(5) << r;		
		}
		else
			cout << setw(5) << "-";

	}
	cout << endl << endl;
}



/* 
=========================

	      TEST 				

========================= 
*/

// Print values of all variables to console
void clustertree::print_all_vals(){
	cout << "PRINTING ALL MEMBER VARIABLES" << endl;

	cout << "== INTEGERS:" << endl;
	cout << "N = " << N << endl;
	cout << "NDIM = " << NDIM << endl;
	cout << "NNN = " << NNN << endl;


	cout << endl << "== POINTERS:" << endl;
	cout << "lattice ptr:" << endl;
	for (long long int i=0; i<N; i++)
		cout << "ptr[" << i << "] = " << ptr[i] << endl;

	cout << "lattice ptr:" << endl;
	for (long long int i=0; i<N; i++)
		cout << "lattice[" << i << "] = " << lattice[i] << endl;

	cout << "nn ptr:" << endl;
	for (long long int i=0; i<N; i++){
		cout << "nn[" << i << "][:] = ";
		for (int j=0; j<NNN; j++)
			cout << nn[i][j] << " ";
		cout << endl;
	}

	cout << endl << endl << "std. lattice:";
	this->print_std_lattice();
	cout << endl << endl << "relabelled lattice:"; 
	this->print_rl_lattice();
}


// test cluster merging algorithm, with diagnostic outputs
void clustertree::merge_clusters_test(){
	long long int i,j;
	long long int s1,s2;
	long long int r1,r2;
	int big = 0;
	int nn_tmp;

	// # of function calls
	int kf = 0;

	for (i=0; i<N; i++){	
		// only continue if on an occupied lattice site
		if (lattice[i] == 1){
			s1 = i;
			r1 = this->findroot(s1,kf);

			for (j=0; j<NNN; j++){
				s2 = nn[i][j];

				// skip nn if lattice is empty
				if (lattice[s2] == 1){
					// get root of adj pt
					r2 = this->findroot(s2,kf);					

					// if diff roots, then need to merge (if same, already merged!)
					if (r2 != r1){
						cout << "rl lattice, i = " << i;
						cout << ", s2 = " << s2;
						cout << ", r1 = " << r1;
						cout << ", r2 = " << r2;
						cout << ". ptr[r1] = " << ptr[r1];
						cout << ". ptr[r2] = " << ptr[r2];

						// if r2 > r1 (mind minus sign), merge 1 -> 2
						if (ptr[r1] > ptr[r2]){
							ptr[r2] += ptr[r1];
							ptr[r1] = r2;
							r1 = r2;
						}																		
						else{
							// else, merge r2 -> r1
							ptr[r1] += ptr[r2];
							ptr[r2] = r1;
						}

						// if new cluster is max, increase max!
						if (-ptr[r1] > big)
							big = -ptr[r1];

						cout << ", big = " << big;
						cout << ": " << endl;
						this->print_rl_lattice();
					}
				}
			}
		}
	}

	smax = big;
	fcalls = kf;
}

