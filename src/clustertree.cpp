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

