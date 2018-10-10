/*

	Methods implementation 
	for voidcluster class

	BY Jack Treado

*/

#include "clustertree.h"
#include "voidcluster.h"
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

const double PI = 3.1415926;


/* 
==================================

	CONSTRUCTORS & DESTRUCTORS		 

================================== 
*/

voidcluster::voidcluster(int np, int l, int ndim, int nnn) : clustertree(l, ndim, nnn){
	#ifdef DEBUG
		cout << "in constructor, L = 1..." << endl;
	#endif

	// dynamically allocate B,pos,rad
	B = new double[ndim];
	pos = new double*[np];
	rad = new double[np];
	
	int i,d;

	// if no input arguments, set box length to 1		
	for(d=0; d<ndim; d++)
		B[d] = 1.0;

	// initialize pos, rad to zero
	for (i=0; i<np; i++){
		rad[i] = 0.0;
		pos[i] = nullptr;
		pos[i] = new double[ndim];
		for (d=0; d<ndim; d++)
			pos[i][d] = 0.0;		
	}

	// set grid spacing
	g = 1.0/l;

	// set other private variables
	NP = np;

	#ifdef DEBUG
		cout << "leaving constructor." << endl << endl;
	#endif

	// no cell, so set all else to nothing
	NCL = -1;
	NCELLS = -1;
	sp = nullptr;
	cellneighbors = nullptr;
	cellpos = nullptr;
	pclabel = nullptr;
	sclabel = nullptr;	
	pcell = nullptr;
	scell = nullptr;
}



voidcluster::voidcluster(int np, int l, int ndim, int nnn, int larr[]) : clustertree(l, ndim, nnn){	
	// dynamically allocate L
	B = new double[ndim];
	pos = new double*[np];
	rad = new double[np];

	int i,d;

	// set box length to larr values	
	for(d=0; d<ndim; d++)
		B[d] = larr[d];

	// initialize pos, rad to zero
	for (i=0; i<np; i++){
		rad[i] = 0.0;
		pos[i] = nullptr;
		pos[i] = new double[ndim];
		for (d=0; d<ndim; d++)
			pos[i][d] = 0.0;		
	}

	// set grid spacing
	g = 1.0/l;

	// set other private variables
	NP = np;

	// no cell, so set all else to nothing
	NCL = -1;
	NCELLS = -1;
	sp = nullptr;
	cellneighbors = nullptr;
	cellpos = nullptr;
	pclabel = nullptr;
	sclabel = nullptr;	
	pcell = nullptr;
	scell = nullptr;
}

// CELL LIST CONSTRUCTOR
voidcluster::voidcluster(int np, int l, int ndim, int nnn, int nc) : clustertree(l, ndim, nnn){	
	long long int N = this->get_N();

	// dynamically allocate L
	B = new double[ndim];
	pos = new double*[np];
	rad = new double[np];

	// set other private variables
	NP = np;

	// initialize cell info
	NCL = nc;
	NCELLS = pow(nc,ndim);
	sp = new double[ndim];
	cellneighbors = new int*[NCELLS];
	cellpos = new double*[NCELLS];
	pclabel = new int[np];
	sclabel = new int[N];	
	pcell = new vector<int>[NCELLS];
	scell = new vector<int>[NCELLS];

	long long int i;
	int d;

	if (ndim == 2)
		NCN = 8;
	else if (ndim == 3)
		NCN = 26;
	else{
		cout << "ERROR: NDIM = " << ndim << " not yet supported...." << endl;
		throw "NDIM not supported\n";
	}

	// set box length to larr values	
	for(d=0; d<ndim; d++){		
		B[d] = 1.0;
		sp[d] = B[d]/nc;
	}

	// initialize pos, rad to zero
	for (i=0; i<np; i++){
		rad[i] = 0.0;
		pos[i] = nullptr;
		pos[i] = new double[ndim];
		for (d=0; d<ndim; d++)
			pos[i][d] = 0.0;

		pclabel[i] = -1;	
	}

	// initilialize site label	
	for (i=0; i<N; i++)
		sclabel[i] = -1;	

	for (i=0; i<NCELLS; i++){
		cellpos[i] = new double[ndim];
		for (d=0; d<ndim; d++)
			cellpos[i][d] = -1;

		cellneighbors[i] = new int[NCN];
		for (d=0; d<NCN; d++)
			cellneighbors[i][d] = -1;
	}

	// set grid spacing
	g = 1.0/l;	
}

// Read in particle positions with cell list
voidcluster::voidcluster(string& spherestr, int sph, int l, int ndim, int nnn, int nc) : clustertree(l, ndim, nnn){	
	long long int N = this->get_N();

	// load in particle info
	ifstream obj(spherestr.c_str());
	if (!obj.is_open()){
		cout << "particle data file " << spherestr << " does not open, ending." << endl;
		throw "file not open!\n";
	}
	cout << "loading particle info from " << spherestr << endl;

	// get number of particles
	int nres = 0;
	obj >> nres;
	cout << "nres = " << nres << endl;

	// get sys size
	double Btmp = 1.0;
	obj >> Btmp >> Btmp >> Btmp;

	// get initial phi
	double phitmp = 0.0;
	obj >> phitmp;
	cout << "phi = " << phitmp << endl;

	// get rid of header
	char ctr;
	obj >> ctr;
	string header(" ");
	getline(obj,header);
	cout << "header #1: " << header << endl;
	getline(obj,header);
	cout << "header #2: " << header << endl;

	// now loading in particle information
	int ind,j,p;
	int k,kmax;
	double tr,radtmp1,xtmp1,ytmp1,ztmp1,mtmp,minrad,scale;
	vector<double> radtmp;
	vector<double> xtmp;
	vector<double> ytmp;
	vector<double> ztmp;

	// loop through file
	p = 0;
	k = 0;
	kmax = 1e6;
	minrad = Btmp;
	cout << "Reading in values:" << endl;
	while(!obj.eof() && k < kmax){
		obj >> ind;
		obj >> radtmp1;
		obj >> xtmp1;
		obj >> ytmp1;
		obj >> ztmp1;
		obj >> mtmp;

		if (p > ind){
			cout << "found end of file, breaking..." << endl;
			break;
		}

		if (radtmp1 < minrad)
			minrad = radtmp1;

		radtmp.push_back(radtmp1);
		xtmp.push_back(xtmp1);
		ytmp.push_back(ytmp1);
		ztmp.push_back(ztmp1);

		// ouput pushback values
		cout << setw(13) << ind;
		cout << setw(13) << radtmp1;
		cout << setw(13) << xtmp1;
		cout << setw(13) << ytmp1;
		cout << setw(13) << ztmp1;
		cout << endl;

		// increment number of particles, number of iterations
		p++;
		k++;
	}
	if (k == kmax){
		cout << "ERROR: infinite loop in reading in from particle file, ending..." << endl;
		throw "infinite loop\n";
	}
	else
		cout << "Read successful! # of atoms = " << p << endl;

	// save number of particles
	NP = p;

	// set scale, scale all lengths system
	scale = minrad;
	cout << "scale = " << scale << endl;
	Btmp = Btmp/scale;
	
	if (ndim == 3)
		NCN = 26;
	else{
		cout << "ERROR: NDIM = " << ndim << " not yet supported for cell neighbor list when reading in particle info...." << endl;
		throw "NDIM not supported\n";
	}

	// set box length to L values	
	this->initialize_box(Btmp);

	// initialize pos and rad
	this->initialize_particles();

	// initialize pos, rad with particle info
	long long int i;
	for (i=0; i<NP; i++){
		// input particle radius
		rad[i] = radtmp.at(i)/scale;

		// input particle position (ONLY WORKS FOR NDIM = 3)
		pos[i][0] = xtmp.at(i)/scale;
		pos[i][1] = ytmp.at(i)/scale;
		pos[i][2] = ztmp.at(i)/scale;
	}

	// get cell info
	NCL = nc;
	NCELLS = pow(nc,ndim);
	this->initialize_cells();
	
	// set grid spacing
	g = Btmp/l;	

	// clear vectors
	xtmp.clear();
	ytmp.clear();
	ztmp.clear();

	// close fstream object
	obj.close();
}


// Read in particle positions with cell list
voidcluster::voidcluster(string& pstr, int l, int ndim, int nnn, int nc) : clustertree(l, ndim, nnn){	
	long long int N = this->get_N();

	// load in particle info
	ifstream obj(pstr.c_str());
	if (!obj.is_open()){
		cout << "particle data file " << pstr << " does not open, ending." << endl;
		throw "file not open!\n";
	}
	cout << "loading particle info from " << pstr << endl;

	// get header
	// get number of particles
	int nres = 0;
	obj >> nres;
	cout << "nres = " << nres << endl;

	// get sys size
	double Btmp = 1.0;
	obj >> Btmp >> Btmp >> Btmp;
	cout << "L = " << Btmp << endl;

	// get initial phi
	double phitmp = 0.0;
	obj >> phitmp;
	cout << "phi = " << phitmp << endl;

	// get rid of header
	char ctr;
	obj >> ctr;
	string header(" ");
	getline(obj,header);
	cout << "header 1: " << header << endl;
	getline(obj,header);
	cout << "header 2: " << header << endl;

	// now loading in particle information
	int ind,j,p,NEXTRA;
	int k,kmax;
	double tr,radtmp1,xtmp1,ytmp1,ztmp1,minrad,scale;
	vector<double> radtmp;
	vector<double> xtmp;
	vector<double> ytmp;
	vector<double> ztmp;

	// loop through file
	p = 0;
	NEXTRA = 12;		// number of extra entries we don't care about
	k = 0;
	kmax = 1e6;
	minrad = Btmp;
	cout << "Reading in values:" << endl;
	while(!obj.eof() && k < kmax){
		obj >> ind;
		obj >> radtmp1;
		obj >> xtmp1;
		obj >> ytmp1;
		obj >> ztmp1;
		for (j=0; j<NEXTRA; j++)
			obj >> tr;

		if (radtmp1 < minrad)
			minrad = radtmp1;

		radtmp.push_back(radtmp1);
		xtmp.push_back(xtmp1);
		ytmp.push_back(ytmp1);
		ztmp.push_back(ztmp1);

		// ouput pushback values
		cout << setw(13) << ind;
		cout << setw(13) << radtmp1;
		cout << setw(13) << xtmp1;
		cout << setw(13) << ytmp1;
		cout << setw(13) << ztmp1;
		cout << endl;

		// increment number of particles, number of iterations
		p++;
		k++;
	}
	if (k == kmax){
		cout << "ERROR: infinite loop in reading in from particle file, ending..." << endl;
		throw "infinite loop\n";
	}
	else
		cout << "Read successful! # of atoms = " << p << endl;

	// save number of particles
	NP = p;

	// set scale, scale all lengths system
	scale = minrad;
	cout << "scale = " << scale << endl;
	Btmp = Btmp/scale;
	
	if (ndim == 3)
		NCN = 26;
	else{
		cout << "ERROR: NDIM = " << ndim << " not yet supported for cell neighbor list when reading in particle info...." << endl;
		throw "NDIM not supported\n";
	}

	// set box length to L values	
	this->initialize_box(Btmp);

	// initialize pos and rad
	this->initialize_particles();

	// initialize pos, rad with particle info
	long long int i;
	for (i=0; i<NP; i++){
		// input particle radius
		rad[i] = radtmp.at(i)/scale;

		// input particle position (ONLY WORKS FOR NDIM = 3)
		pos[i][0] = xtmp.at(i)/scale;
		pos[i][1] = ytmp.at(i)/scale;
		pos[i][2] = ztmp.at(i)/scale;
	}

	// get cell info
	NCL = nc;
	NCELLS = pow(nc,ndim);
	this->initialize_cells();
	
	// set grid spacing
	g = Btmp/l;	

	// clear vectors
	xtmp.clear();
	ytmp.clear();
	ztmp.clear();

	// close fstream object
	obj.close();
}

voidcluster::~voidcluster(){
	int i;

	cout << "destroying voidcluster ptrs ..." << endl;
	delete [] B;
	delete [] rad;
	for (i=0; i<NP; i++)
		delete [] pos[i];

	delete [] pos;
	delete [] pclabel;
	delete [] sclabel;
	delete [] sp;

	for (i=0; i<NCELLS; i++){
		delete [] cellpos[i];
		delete [] cellneighbors[i];
	}
	delete [] cellpos;
	delete [] cellneighbors;

	// clear cell arrays
	this->clear_cells();
	delete [] pcell;
	delete [] scell;

	B = nullptr;
	rad = nullptr;
	pos = nullptr;
	sp = nullptr;
	cellpos = nullptr;
	cellneighbors = nullptr;
	pclabel = nullptr;
	sclabel = nullptr;
	pcell = nullptr;
	scell = nullptr;

	cout << "== voidcluster POINTERS AFTER voidcluster DESTRUCTOR CALL:" << endl;
	cout << "B = " << B << endl;
	cout << "rad = " << rad << endl;
	cout << "pos = " << pos << endl;
}


/* 
==================================

		SETTERS & GETTERS		 

================================== 
*/

void voidcluster::initialize_box(double b){
	int d,NDIM;
	NDIM = this->get_NDIM();

	B = new double[NDIM];
	for(d=0; d<NDIM; d++)		
		B[d] = b;
}

void voidcluster::initialize_particles(){
	int i,d,NDIM;
	NDIM = this->get_NDIM();

	rad = new double[NP];
	pos = new double*[NP];

	for (i=0; i<NP; i++){
		rad[i] = 0.0;
		pos[i] = new double[NDIM];
		for (d=0; d<NDIM; d++)
			pos[i][d] = 0.0;
	}
}


void voidcluster::initialize_cells(){
	int i,d,N,NDIM;
	NDIM = this->get_NDIM();
	N = this->get_N();

	// initialize cell info	
	sp = new double[NDIM];
	cellneighbors = new int*[NCELLS];
	cellpos = new double*[NCELLS];
	pclabel = new int[NP];
	sclabel = new int[N];	
	pcell = new vector<int>[NCELLS];
	scell = new vector<int>[NCELLS];

	// initilialize site label	
	for (i=0; i<N; i++)
		sclabel[i] = -1;	

	for (i=0; i<NCELLS; i++){
		cellpos[i] = new double[NDIM];
		for (d=0; d<NDIM; d++)
			cellpos[i][d] = -1;

		cellneighbors[i] = new int[NCN];
		for (d=0; d<NCN; d++)
			cellneighbors[i][d] = -1;
	}

	for (i=0; i<NP; i++)
		pclabel[i] = -1;

	// set box length to larr values	
	for(d=0; d<NDIM; d++)
		sp[d] = B[d]/NCL;
}






/* 
==================================

	VOID PERCOLATION METHODS		 

================================== 
*/


// randomly place NP spheres into NDIM box, get lattice values
void voidcluster::rand_sphere_pos(int seed){
	// get important clustertree variables
	int NDIM = this->get_NDIM();

	// declare local variables
	int i,d;

	// set initial seed
	srand48(seed);

	// get random initial positions, save to private variable
	for (i=0; i<NP; i++){
		// random positions
		for (d=0; d<NDIM; d++)
			pos[i][d] = B[d]*drand48();
	}
}


// set radii of particles, get void lattice
void voidcluster::rand_sphere_void(double a){
	#ifdef DEBUG
		cout << "void space for random spheres of radius a = " << a << "..." << endl;
	#endif

	// reset system
	this->reset_sys();

	// get important clustertree variables
	int NDIM = this->get_NDIM();
	int L = this->get_L();
	long long int N = this->get_N();	

	// declare local variables
	long long int i;
	int ov;

	// set particle radius a
	this->set_radius(a);

	// loop over void lattice
	for (i=0; i<N; i++){
		// reset ptr to -1
		this->set_ptr(i,-1);

		// check if i overlaps with atoms
		ov = this->check_rad_overlap(i);
		if (ov == 0)
			this->set_lattice(i,1);
		else
			this->set_lattice(i,0);
	}
}

// set radii of particles, get particle lattice
void voidcluster::rand_sphere_pperc(double a){
	// reset system
	this->reset_sys();

	// get important clustertree variables
	int NDIM = this->get_NDIM();
	int L = this->get_L();
	long long int N = this->get_N();	

	// declare local variables
	long long int i;
	int ov;

	// set particle radius a
	this->set_radius(a);

	// loop over void lattice
	for (i=0; i<N; i++){
		// reset ptr to -1
		this->set_ptr(i,-1);

		// check if i overlaps with atoms
		ov = this->check_rad_overlap(i);
		if (ov == 1)
			this->set_lattice(i,1);
		else
			this->set_lattice(i,0);
	}
}

void voidcluster::probe_void(double a){
	// reset system
	this->reset_sys();

	// get important clustertree variables
	int NDIM = this->get_NDIM();
	int L = this->get_L();
	long long int N = this->get_N();	

	// declare local variables
	long long int i;
	int ov;

	// loop over void lattice
	for (i=0; i<N; i++){
		// reset ptr to -1
		this->set_ptr(i,-1);

		// check if i overlaps with atoms
		ov = this->check_probe_overlap(i,a);
		if (ov == 0)
			this->set_lattice(i,1);
		else
			this->set_lattice(i,0);
	}
}



// check for radius overlap with particles in system
int voidcluster::check_rad_overlap(int site){
	// get important clustertree variables
	int NDIM = this->get_NDIM();
	int L = this->get_L();

	// declare local variables
	int d,pp,p,NPit,c,i;
	double dr,h;

	// get position in space
	double sposd;

	// get vector of particles to check
	vector<int> pcheck;

	// check distances to every particle
	if (NCL > 0){
		// get number of particles in cell of lattice site "site"
		c = sclabel[site];
		NPit = pcell[c].size();
		for (i=0; i<NPit; i++)
			pcheck.push_back(pcell[c].at(i));
	}
	else{
		NPit = NP;
		for (i=0; i<NPit; i++)
			pcheck.push_back(i);
	}

	for (pp=0; pp<NPit; pp++){
		h = 0;
		p = pcheck.at(pp);

		// use pythagorean thm. to get distance
		for(d=0; d<NDIM; d++){
			sposd = g*floor((site % (int)pow(L,d+1))/(pow(L,d)));
			dr = pos[p][d] - sposd;
			dr = dr - B[d]*round(dr/B[d]);
			h += dr*dr;
		}		
		h = sqrt(h);

		// if h < radius of p, then inside particle (ov = 1), and can return value to function
		if (h < rad[p]){
			pcheck.clear();
			return 1;
		}
		
	}
	// clear vector just in case
	pcheck.clear();

	// if you got this far, then no particle overlaps (ov = 0), and can return value to function
	return 0;	
}

// check for radius overlap with particles in system
int voidcluster::check_probe_overlap(int site, double a){
	// get important clustertree variables
	int NDIM = this->get_NDIM();
	int L = this->get_L();

	// declare local variables
	int d,pp,p,NPit,c,i,s1,s2;
	double dr,h,prober;

	// get position in space
	double sposd;

	// get vector of particles to check
	vector<int> pcheck;

	// check distances to every particle
	if (NCL > 0){
		// get number of particles in cell of lattice site "site"
		c = sclabel[site];
		NPit = pcell[c].size();
		for (i=0; i<NPit; i++)
			pcheck.push_back(pcell[c].at(i));		
	}
	else{
		NPit = NP;
		for (i=0; i<NPit; i++)
			pcheck.push_back(i);
	}

	for (pp=0; pp<NPit; pp++){
		h = 0;
		p = pcheck.at(pp);

		// use pythagorean thm. to get distance
		for(d=0; d<NDIM; d++){
			s1 = pow(L,d+1);
			s2 = pow(L,d);
			sposd = g*floor((site % s1)/s2);

			dr = pos[p][d] - sposd;
			dr = dr - B[d]*round(dr/B[d]);
			h += dr*dr;
		}		
		h = sqrt(h);

		// get probe radius
		prober = rad[p]+a;

		// if h < probe radius, then inside particle (ov = 1), and can return value to function
		if (h < prober){
			pcheck.clear();
			return 1;
		}
		
	}
	// clear vector just in case
	pcheck.clear();

	// if you got this far, then no particle overlaps (ov = 0), and can return value to function
	return 0;	
}




// use binary root search to find percolation
void voidcluster::find_perc(double aH, double aL, double epsilon){
	#ifdef DEBUG
		cout << "starting binary search to get ac..." << endl;
	#endif

	// get important clustertree variables
	int NDIM = this->get_NDIM();

	// declare local variables
	int k,kmax,perc;
	long long int lsum;
	double a,ap,check,check_new;	

	// while loop iterator
	k = 0;
	kmax = 1e4;

	// termination criterion
	a = 0.5*(aH+aL);
	check = 10.0*epsilon;
	check_new = check;

	// output information to console
	int w = 15;
	int prc = 3;
	double poro; 
	cout << setw(w) << "k";
	cout << setw(w) << "a";
	cout << setw(w) << "aH";
	cout << setw(w) << "aL";
	cout << setw(w) << "check";
	cout << setw(w) << "~poro";
	cout << setw(w) << "vctot";
	cout << setw(w) << "perc";
	cout << setw(w) << "cnum";
	cout << setw(w) << "smax";
	cout << setw(w) << "lsites";
	cout << setw(w) << "fcalls";
	cout << endl;
	for (int e=0; e<12*w; e++)
		cout << "=";
	cout << endl;	

	// loop until converged on percolation probability and make sure cluster is percolated
	while (check >= epsilon && k < kmax){
		k++;
		check = a;		

		// get lattice with new radius
		this->rand_sphere_void(a);		

		// merge clusters, check for percolation
		this->merge_clusters();
		this->post_process();
		perc = this->get_perc();		

		// output information to console
		poro = exp(-NP*(1.33333333333)*PI*pow(a,NDIM));
		lsum = this->get_sum_s();

		cout << setprecision(prc) << setw(w) << k;
		cout << setprecision(prc) << setw(w) << a;
		cout << setprecision(prc) << setw(w) << aH;
		cout << setprecision(prc) << setw(w) << aL;
		cout << setprecision(prc) << setw(w) << check_new;
		cout << setprecision(prc) << setw(w) << poro;
		cout << setprecision(prc) << setw(w) << pow(g,NDIM)*lsum;
		cout << setprecision(prc) << setw(w) << perc;
		cout << setprecision(prc) << setw(w) << this->get_cnum();
		cout << setprecision(prc) << setw(w) << this->get_smax();
		cout << setprecision(prc) << setw(w) << this->get_lattice_sum();
		cout << setprecision(prc) << setw(w) << this->get_fcalls();
		cout << endl;	

		// check percolation, update aH or aL
		if (perc ==	1){
			aL = a;
			ap = a;
		}
		else
			aH = a;

		// update check, a
		a = 0.5*(aL+aH);
		check = abs(check-a)/a;	
		check_new = check;	
	}
	if (k < kmax && check < epsilon){
		cout << "percolation found!";
		ac = ap;
		this->rand_sphere_void(ac);		

		// merge clusters, check for percolation
		this->merge_clusters();
		this->post_process();
		perc = this->get_perc();
		lsum = this->get_sum_s();
		cout << "; lsum = " << lsum;
		vctot = (double)lsum;
		vctot *= pow(g,NDIM);
		cout << "; ac = " << ac << "; vctot = " << vctot << "; g = " << g << endl;
	}
	else{
		cout << "percolation not found in k < kmax :(" << endl;
		ac = -1;
		vctot = -1;
	}
}

void voidcluster::find_particle_perc(double aH, double aL, double aC, double epsilon){
	// get important clustertree variables
	int NDIM = this->get_NDIM();

	// declare local variables
	int k,kmax,perc;
	long long int lsum;
	double a,ap,check,check_new;	

	// while loop iterator
	k = 0;
	kmax = 1e4;

	// termination criterion
	a = aC;
	check = 10.0*epsilon;
	check_new = check;

	// output information to console
	int w = 15;
	int prc = 3;
	double pack; 
	cout << setw(w) << "k";
	cout << setw(w) << "a";
	cout << setw(w) << "aH";
	cout << setw(w) << "aL";
	cout << setw(w) << "check";
	cout << setw(w) << "~pack";
	cout << setw(w) << "vctot";
	cout << setw(w) << "perc";
	cout << setw(w) << "cnum";
	cout << setw(w) << "smax";
	cout << setw(w) << "lsites";
	cout << setw(w) << "fcalls";
	cout << endl;
	for (int e=0; e<12*w; e++)
		cout << "=";
	cout << endl;

	// loop until converged on percolation probability and make sure cluster is percolated
	while (check >= epsilon && k < kmax){
		k++;
		check = a;		

		// get lattice with new radius		
		this->rand_sphere_pperc(a);		

		// merge clusters, check for percolation
		this->merge_clusters();
		this->post_process();
		perc = this->get_perc();


		// output information to console
		pack = 1-exp(-NP*(1.33333333333)*PI*pow(a,NDIM));
		lsum = this->get_sum_s();

		cout << setprecision(prc) << setw(w) << k;
		cout << setprecision(prc) << setw(w) << a;
		cout << setprecision(prc) << setw(w) << aH;
		cout << setprecision(prc) << setw(w) << aL;
		cout << setprecision(prc) << setw(w) << check_new;
		cout << setprecision(prc) << setw(w) << pack;
		cout << setprecision(prc) << setw(w) << pow(g,NDIM)*lsum;
		cout << setprecision(prc) << setw(w) << perc;
		cout << setprecision(prc) << setw(w) << this->get_cnum();
		cout << setprecision(prc) << setw(w) << this->get_smax();
		cout << setprecision(prc) << setw(w) << this->get_lattice_sum();
		cout << setprecision(prc) << setw(w) << this->get_fcalls();
		cout << endl;	

		// check percolation, update aH or aL
		if (perc ==	1){
			aH = a;
			ap = a;
		}
		else
			aL = a;

		// update check, a
		a = 0.5*(aL+aH);
		check = abs(check-a)/a;	
		check_new = check;	
	}
	if (k < kmax && check < epsilon){
		cout << "percolation found!";
		ac = ap;
		this->rand_sphere_pperc(ac);		

		// merge clusters, check for percolation
		this->merge_clusters();
		this->post_process();
		perc = this->get_perc();
		lsum = this->get_sum_s();
		cout << "; lsum = " << lsum;
		vctot = (double)lsum;
		vctot *= pow(g,NDIM);
		cout << "; ac = " << ac << "; vctot = " << vctot << "; g = " << g << endl;
	}
	else{
		cout << "percolation not found in k < kmax :(" << endl;
		ac = -1;
		vctot = -1;
	}
}

void voidcluster::find_probe_perc(double aH, double aL, double aC, double epsilon){
	// get important clustertree variables
	int NDIM = this->get_NDIM();

	// declare local variables
	int k,kmax,perc;
	long long int lsum;
	double a,ap,check,check_new;	

	// while loop iterator
	k = 0;
	kmax = 1e4;

	// termination criterion
	a = aC;
	check = 10.0*epsilon;
	check_new = check;
	perc = 0;

	// output information to console
	int w = 15;
	int prc = 3;
	cout << setw(w) << "k";
	cout << setw(w) << "a";
	cout << setw(w) << "aH";
	cout << setw(w) << "aL";
	cout << setw(w) << "check";
	cout << setw(w) << "vctot";
	cout << setw(w) << "perc";
	cout << setw(w) << "cnum";
	cout << setw(w) << "smax";
	cout << setw(w) << "lsites";
	cout << setw(w) << "fcalls";
	cout << endl;
	for (int e=0; e<11*w; e++)
		cout << "=";
	cout << endl;	

	// loop until converged on percolation probability and make sure cluster is percolated
	while( ((check > epsilon) || perc == 0) && k < kmax){
		k++;
		check = a;		

		// get lattice with new radius
		this->reset_sys();
		this->reset_ptr();
		this->probe_void(a);		

		// merge clusters, check for percolation
		this->merge_clusters();
		this->post_process();
		perc = this->get_perc();

		// output information to console
		lsum = this->get_sum_s();

		cout << setprecision(prc) << setw(w) << k;
		cout << setprecision(prc) << setw(w) << a;
		cout << setprecision(prc) << setw(w) << aH;
		cout << setprecision(prc) << setw(w) << aL;
		cout << setprecision(prc) << setw(w) << check_new;
		cout << setprecision(prc) << setw(w) << pow(g,NDIM)*lsum/(B[0]*B[1]*B[2]);
		cout << setprecision(prc) << setw(w) << perc;
		cout << setprecision(prc) << setw(w) << this->get_cnum();
		cout << setprecision(prc) << setw(w) << this->get_smax();
		cout << setprecision(prc) << setw(w) << this->get_lattice_sum();
		cout << setprecision(prc) << setw(w) << this->get_fcalls();
		cout << endl;	

		// check percolation, update aH or aL
		if (perc ==	0)
			aH = a;
		else{
			aL = a;
			ap = a;
		}

		// update check, a
		a = 0.5*(aL+aH);
		check = abs(check-a)/a;	
		check_new = check;	
	}
	if (k < kmax && check < epsilon){
		cout << "percolation found!";
		ac = ap;

		lsum = this->get_sum_s();
		cout << "; lsum = " << lsum;
		vctot = (double)lsum;
		vctot *= pow(g,NDIM);
		vctot /= B[0]*B[1]*B[2];
		cout << "; ac = " << ac << "; vctot = " << vctot << "; g = " << g << endl;

		// plot if uncommented
		// string xyzf = "res_cvoids_final.xyz";
		// ofstream xyzobj(xyzf.c_str());
		// this->print_res_xyz(xyzobj);
		// xyzobj.close();
	}
	else{
		cout << "percolation not found in k < kmax :(" << endl;
		ac = -1;
		vctot = -1;
	}
}



void voidcluster::fprint_stats(ofstream& obj, int seed, int header){
	int w = 16;

	// output optional header
	if (header == 1){
		obj << setw(w) << "# seed";
		obj << setw(w) << "ac";
		obj << setw(w) << "vctot";
		obj << setw(w) << "smax";
		obj << setw(w) << "smean";
		obj << setw(w) << "cnum";
		obj << setw(w) << "fcalls";
		obj << endl;
		for (int e=0; e<7*w; e++)
			obj << "=";
		obj << endl; 		
	}

	// output variables
	obj << setw(w) << seed;
	obj << setw(w) << ac;
	obj << setw(w) << vctot;
	obj << setw(w) << this->get_smax();
	obj << setw(w) << this->get_smean();
	obj << setw(w) << this->get_cnum();
	obj << setw(w) << this->get_fcalls();
	obj << endl;
}

void voidcluster::fprint_stats(ofstream& obj, int header){
	int w = 16;

	// output optional header
	if (header == 1){
		obj << setw(w) << "ac";
		obj << setw(w) << "vctot";
		obj << setw(w) << "smax";
		obj << setw(w) << "smean";
		obj << setw(w) << "cnum";
		obj << setw(w) << "fcalls";
		obj << endl;
		for (int e=0; e<6*w; e++)
			obj << "=";
		obj << endl; 		
	}

	// output variables
	obj << setw(w) << ac;
	obj << setw(w) << vctot;
	obj << setw(w) << this->get_smax();
	obj << setw(w) << this->get_smean();
	obj << setw(w) << this->get_cnum();
	obj << setw(w) << this->get_fcalls();
	obj << endl;
}


void voidcluster::fprint_s(ofstream& obj){
	this->print_cluster_sizes(obj);
}


void voidcluster::print_lattice_positions(){
	cout << "PRINTING LATTICE POSITIONS MAPPED TO INDEX SPACE" << endl;
	long long int N = this->get_N();
	int NDIM = this->get_NDIM();
	int L = this->get_L();

	cout << "L = " << L << "; N = " << N << "; g = " << g << endl << endl;


	cout << setw(10) << "i";
	cout << setw(10) << "x";
	cout << setw(10) << "y";
	cout << setw(10) << "z";
	cout << endl;
	cout << "==============================================";
	cout << endl;

	double site_pos[NDIM];
	for (long long int i=0; i<N; i++){
		cout << setprecision(3) << setw(10) << i;
		for(int d=0; d<NDIM; d++){
			if (d == 0)
				site_pos[d] = g*(i % L);
			else
				site_pos[d] = g*(floor(i/pow(L,d)));

			cout << setprecision(3) << setw(10) << site_pos[d];
		}		
		cout << endl;
	}
}


void voidcluster::print_void_xyz(ofstream &obj){
	cout << "Printing both particle and void to .xyz file, read in Ovito..." << endl;

	// get relevant variables from base class
	long long int N = this->get_N();
	int NDIM = this->get_NDIM();
	int L = this->get_L();

	// local variables
	int d,dd,w,p;
	long long int i,Ntot;

	// output column width and precision
	w = 15;
	p = 5;

	// get total number of void points and particles
	Ntot = (long long int)NP + this->get_lattice_sum();

	// print .xyz header
	obj << Ntot << endl;
	obj << "Lattice=\"";
	for (d=0; d<NDIM; d++){
		for(dd=0; dd<NDIM; dd++){
			if (dd==d)
				obj << B[d];
			else
				obj << " 0.0 ";
		}
	}
	obj << "\" ";
	obj << '\t';
	obj << "Properties=species:S:1:pos:R:" <<  NDIM << ":radius:R:1" << endl;

	// print particle positions
	char sp = 'C';
	for (i=0; i<NP; i++){
		obj << setw(w) << sp;
		for (d=0; d<NDIM; d++)
			obj << setw(w) << setprecision(p) << pos[i][d];
		obj << setw(w) << setprecision(p) << rad[i];
		obj << endl;
	}

	// print lattice positions
	sp = 'O';
	double la = 0.05*rad[0];
	double sposd;
	int pclus,root;

	// if percolated cluster, label as white
	pclus = this->get_pclus();
	for (i=0; i<N; i++){

		if ((this->get_lattice_site(i))==1){
			// get location of site			
			
			if (pclus > 0){
				root = this->findroot(i);
				if (root == pclus)
					sp = 'H';
				else
					sp = 'O';
			}

			// print lattice info 
			obj << setw(w) << sp;
			for (d=0; d<NDIM; d++){
				sposd = g*floor((i % (int)pow(L,d+1))/(pow(L,d)));
				obj << setw(w) << setprecision(p) << sposd;
			}
			obj << setw(w) << setprecision(p) << la;
			obj << endl;
		}		
	}	
}



void voidcluster::print_res_xyz(ofstream &obj){
	cout << "Printing both residue and probe void to .xyz file, read in Ovito..." << endl;

	// get relevant variables from base class
	long long int N = this->get_N();
	int NDIM = this->get_NDIM();
	int L = this->get_L();

	// local variables
	int d,dd,w,p;
	long long int i,Ntot;

	// output column width and precision
	w = 15;
	p = 5;

	// get total number of void points and particles
	Ntot = (long long int)NP + this->get_lattice_sum();

	// print .xyz header
	obj << Ntot << endl;
	obj << "Lattice=\"";
	for (d=0; d<NDIM; d++){
		for(dd=0; dd<NDIM; dd++){
			if (dd==d)
				obj << B[d];
			else
				obj << " 0.0 ";
		}
	}
	obj << "\" ";
	obj << '\t';
	obj << "Properties=species:S:1:pos:R:" <<  NDIM << ":radius:R:1" << endl;

	// print particle positions
	char sp = 'C';
	int Nfnd = 0;
	double rtmp = 0.0;
	for (i=0; i<NP; i++){
		// determine particle type
		rtmp = rad[i];

		if (rtmp < 1.3)
			sp = 'H';
		else if (rtmp < 1.31){
			if (Nfnd == 0){
				sp = 'N';
				Nfnd = 1;
			}
			else{
				sp = 'C';
				Nfnd = 0;
			}
		}
		else if (rtmp < 1.44)
			sp = 'O';
		else if (rtmp < 1.74)
			sp = 'C';
		else
			sp = 'S';


		obj << setw(w) << sp;
		for (d=0; d<NDIM; d++)
			obj << setw(w) << setprecision(p) << pos[i][d];
		obj << setw(w) << setprecision(p) << rtmp;
		obj << endl;
	}

	// print lattice positions
	sp = 'O';
	double la = ac;
	double sposd;
	int pclus,root;

	// if percolated cluster, label as white
	pclus = this->get_pclus();
	for (i=0; i<N; i++){

		if ((this->get_lattice_site(i))==1){					
			if (pclus > 0){
				root = this->findroot(i);
				if (root == pclus)
					sp = 'X';
				else
					sp = 'Y';
			}

			// print lattice info 
			obj << setw(w) << sp;
			for (d=0; d<NDIM; d++){
				// get location of site	
				sposd = g*floor((i % (int)pow(L,d+1))/(pow(L,d)));
				obj << setw(w) << setprecision(p) << sposd;
			}
			obj << setw(w) << setprecision(p) << la;
			obj << endl;
		}		
	}	
}





/* 
==================================

		CELL LIST METHODS		 

================================== 
*/


void voidcluster::add_pcell(int c, int p){
	pcell[c].push_back(p);
}

void voidcluster::add_scell(int c, int s){
	scell[c].push_back(s);
}

void voidcluster::clear_cells(){
	int i;
	for (i=0; i<NCELLS; i++){
		pcell[i].clear();
		scell[i].clear();
	}
}

void voidcluster::get_cell_positions(){
	int c,d;
	int ind;

	int NDIM;
	NDIM = this->get_NDIM();	

	// loop over cells, get pos of each one
	for (c=0; c<NCELLS; c++){
		for (d=0; d<NDIM; d++){
			// get indices xi,yi,zi,wi,etc
			ind = floor((c % (int)pow(NCL,d+1))/(pow(NCL,d)));

			// get position in given dimension of cell
			cellpos[c][d] = sp[d]*(ind+0.5);
		}
	}
}

void voidcluster::get_cell_neighbors(){
	int i,zi;
	int NDIM;
	NDIM = this->get_NDIM();

	if (NDIM == 3){
		for (i=0; i<NCELLS; i++){
			zi = i/(NCL*NCL);

			// faces
			cellneighbors[i][0] = (i+NCELLS-1) % NCELLS; 			// left neighbor (i-1)
			cellneighbors[i][1] = i-NCL;						// bottom neighbor (j-1)
			cellneighbors[i][2] = (i+1) % NCELLS; 				// right neighbor (i+1)
			cellneighbors[i][3] = i+NCL;						// top neighbor (j+1)
			cellneighbors[i][4] = (i+NCELLS-NCL*NCL) % NCELLS;			// backward neighbor (k-1)
			cellneighbors[i][5] = (i+NCL*NCL) % NCELLS;				// forward neighbor (k+1)		

			// y-direction bc
			if (i-zi*NCL*NCL < NCL)
				cellneighbors[i][1] = i-NCL+NCL*NCL;	
			if ((i+NCL)/(NCL*NCL) > zi)
				cellneighbors[i][3] = i-NCL*NCL+NCL;	

			// edges

			// * xy plane
			cellneighbors[i][6] = cellneighbors[i][1] - 1;	// bottom-left neightbor (i-1,j-1)
			cellneighbors[i][7] = cellneighbors[i][1] + 1;	// bottom-right neighbor (i+1,j-1)
			cellneighbors[i][8] = cellneighbors[i][3] - 1; 	// top-left neighbor (i-1,j+1)
			cellneighbors[i][9] = cellneighbors[i][3] + 1;	// top-right neighbor (i+1.j+1)

			// * xz plane
			cellneighbors[i][10] = (cellneighbors[i][4] - 1) % NCELLS;	// back-left neighbor (i-1,k-1)
			cellneighbors[i][11] = (cellneighbors[i][4] + 1) % NCELLS;	// back-right neighbor (i+1,k-1)
			cellneighbors[i][12] = (cellneighbors[i][5] - 1) % NCELLS;	// front-left neighbor (i-1,k+1)		
			cellneighbors[i][13] = (cellneighbors[i][5] + 1) % NCELLS;	// back-front neighbor (i+1,k+1)

			// * yz plane
			cellneighbors[i][14] = cellneighbors[i][4] - NCL;	// back-bottom neighbor (j-1,k-1)
			cellneighbors[i][15] = cellneighbors[i][4] + NCL;	// back-top neighbor (j+1,k-1)
			cellneighbors[i][16] = cellneighbors[i][5] - NCL;	// front-bottom neighbor (j-1,k+1)
			cellneighbors[i][17] = cellneighbors[i][5] + NCL; 	// front-top neighbor (j+1,k+1)

			// y-direction bc
			if (i-zi*NCL*NCL < NCL){
				cellneighbors[i][14] = cellneighbors[i][4]-NCL+NCL*NCL;
				cellneighbors[i][16] = cellneighbors[i][5]-NCL+NCL*NCL;	
			}
			if ((i+NCL)/(NCL*NCL) > zi){
				cellneighbors[i][15] = cellneighbors[i][4]-NCL*NCL+NCL;
				cellneighbors[i][17] = cellneighbors[i][5]-NCL*NCL+NCL;			
			}
			


			// cubic vertices
			cellneighbors[i][18] = (cellneighbors[i][14] - 1) % NCELLS; // back-bottom-left neighbor (i-1,j-1,k-1)
			cellneighbors[i][19] = (cellneighbors[i][14] + 1) % NCELLS; // back-bottom-right neighbor (i+1,j-1,k-1)
			cellneighbors[i][20] = (cellneighbors[i][15] - 1) % NCELLS; // back-top-left neighbor (i-1,j+1,k-1)
			cellneighbors[i][21] = (cellneighbors[i][15] + 1) % NCELLS; // back-top-right neighbor (i+1,j+1,k-1)
			cellneighbors[i][22] = (cellneighbors[i][16] - 1) % NCELLS; // front-bottom-left neighbor (i-1,j-1,k+1)
			cellneighbors[i][23] = (cellneighbors[i][16] + 1) % NCELLS; // front-bottom-right neighbor (i+1,j-1,k+1)
			cellneighbors[i][24] = (cellneighbors[i][17] - 1) % NCELLS; // front-top-left neighbor (i-1,j+1,k+1)
			cellneighbors[i][25] = (cellneighbors[i][17] + 1) % NCELLS; // front-top-right neighbor (i+1,j+1,k+1)



			if (i % NCL == 0){
				// left BC
				cellneighbors[i][0] = i+NCL-1;
				cellneighbors[i][6] = cellneighbors[i][1]+NCL-1;
				cellneighbors[i][8] = cellneighbors[i][3]+NCL-1;
				cellneighbors[i][10] = cellneighbors[i][4]+NCL-1;
				cellneighbors[i][12] = cellneighbors[i][5]+NCL-1;	
				cellneighbors[i][18] = cellneighbors[i][14]+NCL-1;
				cellneighbors[i][20] = cellneighbors[i][15]+NCL-1;
				cellneighbors[i][22] = cellneighbors[i][16]+NCL-1;
				cellneighbors[i][24] = cellneighbors[i][17]+NCL-1;
			}
			if ((i+1) % NCL == 0){
				// right BC
				cellneighbors[i][2] = i-NCL+1;
				cellneighbors[i][7] = cellneighbors[i][1]-NCL+1;
				cellneighbors[i][9] = cellneighbors[i][3]-NCL+1;
				cellneighbors[i][11] = cellneighbors[i][4]-NCL+1;
				cellneighbors[i][13] = cellneighbors[i][5]-NCL+1;
				cellneighbors[i][19] = cellneighbors[i][14]-NCL+1;
				cellneighbors[i][21] = cellneighbors[i][15]-NCL+1;
				cellneighbors[i][23] = cellneighbors[i][16]-NCL+1;
				cellneighbors[i][25] = cellneighbors[i][17]-NCL+1;
			}
		}
	}	
	else if(NDIM == 2){
		cout << "ERROR: NDIM=2 not supported when using neighborlist..." << endl;
		throw "NDIM not supported!\n";
	}
	else{
		cout << "ERROR: NDIM not supported when using neighborlist..." << endl;
		throw "NDIM not supported!\n";
	}
}

int voidcluster::get_pcell(int p){
	int c,d,minm;
	double cposi;
	double dr,h,mindist;
	int NDIM;
	NDIM = this->get_NDIM();

	mindist = 3*B[0];
	minm = -1;

	// loop over cells, get pos of each one
	for (c=0; c<NCELLS; c++){

		h = 0;
		for (d=0; d<NDIM; d++){
			// get cell pos
			cposi = cellpos[c][d];

			// get distance to particle i
			dr = pos[p][d] - cposi;
			dr = dr - B[d]*round(dr/B[d]);

			h += dr*dr;					
		}
		h = sqrt(h);

		// check if min dist
		if (h < mindist){			
			mindist = h;
			minm = c;
		}
	}	

	// if no min set, error
	if (minm == -1){
		cout << "ERROR: No native cell found, error in get_new_cell()...";
		throw "NO NATIVE CELL FOUND!!\n";
	}
	else
		return minm;
}

int voidcluster::get_scell(int s){
	int c,d,minm;
	double cposd,sposd;
	double dr,h,mindist;
	int L,NDIM;
	L = this->get_L();
	NDIM = this->get_NDIM();

	mindist = 3*B[0];
	minm = -1;

	// loop over cells, get pos of each one
	for (c=0; c<NCELLS; c++){

		h = 0;
		for (d=0; d<NDIM; d++){
			// get cell pos
			cposd = cellpos[c][d];

			// get pos of site s
			sposd = g*floor((s % (int)pow(L,d+1))/(pow(L,d)));

			// get distance to particle i
			dr = sposd - cposd;
			dr = dr - B[d]*round(dr/B[d]);
			h += dr*dr;					
		}

		// check if min dist
		if (h < mindist){					
			mindist = h;
			minm = c;
		}
	}
	// if no min set, error
	if (minm == -1){
		cout << "ERROR: No native cell found, error in get_new_cell()...";
		throw "NO NATIVE CELL FOUND!!\n";
	}
	else
		return minm;
}

void voidcluster::label_cells(){
	int p,d;	
	int cind;
	long long int s,N;
	N = this->get_N();

	this->clear_cells();

	// for all atoms in system, check cell location
	for (p=0; p<NP; p++){
		// check which pcell particle i is in		
		cind = this->get_pcell(p);

		// append particle to new cell
		this->add_pcell(cind,p);
		this->add_pcell_neighbors(cind,p);
		pclabel[p] = cind;
	}

	// for all lattice sites in system, check cell location
	for (s=0; s<N; s++){
		// check which pcell site i is in
		cind = this->get_scell(s);

		// append site to new cell
		this->add_scell(cind,s);
		sclabel[s] = cind;
	}
}

void voidcluster::add_pcell_neighbors(int c, int p){
	long long int i;
	int cn;
	for (i=0; i<NCN; i++){
		cn = cellneighbors[c][i];
		this->add_pcell(cn,p);
	}
}


void voidcluster::print_pcell(){
	cout << "Printing cell info";
	cout << ": NCELLS = " << NCELLS;
	cout << "; NCL = " << NCL;
	cout << "; NCN = " << NCN;
	cout << endl;

	int i,n,m;
	for (i=0; i<NCELLS; i++){
		cout << "cell " << i << ": ";
		n = pcell[i].size();
		for (m=0; m<n; m++)
			cout << setw(5) << pcell[i].at(m);
		cout << endl;
	}
}

void voidcluster::print_cellpos(){
	cout << "Printing cell pos info";
	cout << endl;

	int NDIM = this->get_NDIM();

	int i,n,d;
	for (i=0; i<NCELLS; i++){
		cout << "cell pos " << i << ": ";
		n = pcell[i].size();
		for (d=0; d<NDIM; d++)
			cout << setw(10) << cellpos[i][d];
		cout << endl;
	}
}



void voidcluster::print_cell_xyz(ofstream &obj){
	cout << "Printing both particle and void to .xyz file, read in Ovito..." << endl;

	// get relevant variables from base class
	int NDIM = this->get_NDIM();
	int L = this->get_L();

	// local variables
	int i,d,dd,Ntot,w,p;

	// output column width and precision
	w = 15;
	p = 5;

	// get total number of void points and particles
	Ntot = NP+NCELLS;

	// print .xyz header
	obj << Ntot << endl;
	obj << "Lattice=\"";
	for (d=0; d<NDIM; d++){
		for(dd=0; dd<NDIM; dd++){
			if (dd==d)
				obj << B[d];
			else
				obj << " 0.0 ";
		}
	}
	obj << "\" ";
	obj << '\t';
	obj << "Properties=species:S:1:pos:R:" <<  NDIM << ":radius:R:1" << endl;

	// print particle positions
	char sp;
	int c,n0,nn,f0;
	int midc = 1;
	n0 = this->pcell[midc].size();
	for (i=0; i<NP; i++){
		f0 = 0;
		for (nn=0; nn<n0; nn++){
			if (pcell[midc].at(nn) == i){
				f0 = 1;
				obj << setw(w) << 'C';
				break;
			}
		}

		if (f0 == 0)
			obj << setw(w) << 'N';

		for (d=0; d<NDIM; d++)
			obj << setw(w) << setprecision(p) << pos[i][d];
		obj << setw(w) << setprecision(p) << rad[i];
		obj << endl;
	}	

	sp = 'O';
	for (c=0; c<NCELLS; c++){
		obj << setw(w) << sp;
		for (d=0; d<NDIM; d++){
			obj << setw(w) << setprecision(p) << cellpos[c][d];
		}
		obj << setw(w) << setprecision(p) << 0.05*rad[0];
		obj << endl;
	}


}


void voidcluster::load_zhe_pos(string &zstr){
	int i; 
	double px,py,pz,Ltmp;

	// load zhe's data
	ifstream zobj(zstr.c_str());

	cout << "loading zhe data..." << endl;
	for (i=0; i<NP; i++){
		zobj >> px >> py >> pz;
		pos[i][0] = px;
		pos[i][1] = py;
		pos[i][2] = pz;
		cout << "pos at i = " << i;
		cout << setw(8) << pos[i][0];
		cout << setw(8) << pos[i][1];
		cout << setw(8) << pos[i][2];
		cout << endl;
	}

	zobj.close();
}


















