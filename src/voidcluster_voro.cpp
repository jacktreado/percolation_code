// Voronoi edge cluster percolation implementation

#include "clustertree.h"
#include "voidcluster.h"
#include "voro++.hh"
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
using namespace voro;
const double PI = 3.1415926;

// check if tmp point is contained within vectors of real points
double merge_vertices(double B[], double xtmp, double ytmp, double ztmp, vector<double>& gvx, vector<double>& gvy, vector<double>& gvz){
	double dx,dy,dz,dr;
	int k,sz;
	sz = gvx.size();

	for (k=0; k<sz; k++){
		dx = xtmp - gvx.at(k);
		dx = dx - B[0]*round(dx/B[0]);
		dy = ytmp - gvy.at(k);
		dy = dy - B[1]*round(dy/B[1]);
		dz = ztmp - gvz.at(k);
		dz = dz - B[2]*round(dz/B[2]);

		dr = sqrt(dx*dx + dy*dy + dz*dz);

		// if vertices are very close, then merge
		if (dr < 2e-2)
			return k;
	}	

	// if you got this far, then new vertex found! so add to vectors
	k = sz;
	gvx.push_back(xtmp);
	gvy.push_back(ytmp);
	gvz.push_back(ztmp);
	return k;
}

void voidcluster::find_voro_perc(double epsilon, double seed, double aH, double aL){
	// local variables
	int i,j,n,m;
	int NDIM = 3;	

	// generate voronoi from particle positions

	// Create a container_poly with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// 2 atoms within each computational block.
	n = 2;
	int init = 5;
	bool pbc = true;
	container conp(0,B[0],0,B[1],0,B[2],n,n,n,
			pbc,pbc,pbc,init);

	for (i=0; i<NP; i++)
		conp.put(i,pos[i][0],pos[i][1],pos[i][2]);

	// output postion, edge, vertex info


	string vorostr = "voro_N" + to_string(NP) + "_seed" + to_string(seed) + ".dat";
	// REAL ONE: conp.print_custom("%i %w %g %o %s %t %A %P",vorostr.c_str());
	conp.print_custom("%i %w %g %o %s %a %t %P",vorostr.c_str());

	// read in particle, vertex, edge info
	ifstream vobj(vorostr.c_str());
	vector<double>* vvertx;
	vector<double>* vverty;
	vector<double>* vvertz;
	vector<double>* vvertc;
	vector<int>** vvertn;	
	vector<int>* vfacec;
	int* Nf;
	int** Nv;

	// allocate array of vectors
	vvertx = new vector<double>[NP];
	vverty = new vector<double>[NP];
	vvertz = new vector<double>[NP];
	vvertc = new vector<double>[NP];
	vvertn = new vector<int>*[NP];

	// allocate arrays for number of faces
	Nf = new int[NP];
	Nv = new int*[NP];

	int ind,vn,vnn,vntmp1,vntmp2,invec,en,vc,fn,fe,f;
	double val;
	int pos1,pos2;
	string linestr,numstr;
	stringstream ss("test");

	for (i=0; i<NP; i++){
		// get cell index, number of cell vertices, number of cell edges
		vobj >> ind;
		vobj >> vn;
		vobj >> en;

		// cout << "ind = " << ind << "; ";
		// cout << "vn = " << vn << "; ";
		// cout << "en = " << en << "; ";
		// cout << endl << endl;

		// allocate space for vertex neighbor vectors
		vvertn[ind] = new vector<int>[vn];

		for (j=0; j<vn; j++){
			vobj >> vc;
			vvertc[ind].push_back(vc);
		}
		vobj >> fn;
		// save number of faces per particles
		Nf[ind] = fn;
		Nv[ind] = new int[fn];

		for (j=0; j<fn; j++){
			vobj >> fe;
			Nv[ind][j] = fe;
		}
		vfacec = new vector<int>[fn];

		// save rest of file line into string
		getline(vobj,linestr);

		// parse vertex positions
		pos1 = 0;
		pos2 = 0;
		val = 0;
		// cout << "linestr size = " << linestr.size() << "; pos2 = " << pos2 << " and " << (pos2 < linestr.size()) << endl;		
		
		// cout << "linestr = " << linestr << endl;		
		// loop over faces, vertices per face
		for (f=0; f<fn; f++){
			// get locations in string of beginning of list
			pos1 = linestr.find_first_of("(",pos2+1);
			pos2 = linestr.find_first_of(",",pos1+1);

			// loop over face list, add vertices to vector
			for (j=0; j<Nv[ind][f]; j++){
				numstr = linestr.substr(pos1+1,pos2-pos1-1);
				ss << numstr;
				ss >> vnn;
				ss.clear();

				// add vertex to list of vertices
				vfacec[f].push_back(vnn);

				// get new string positions
				pos1 = pos2;
				if (j<Nv[ind][f]-2)
					pos2 = linestr.find_first_of(",",pos1+1);
				else if (j==Nv[ind][f]-2)
					pos2 = linestr.find_first_of(")",pos1+1);

				// cout << "** particle ind = " << ind << "; face f = " << f << "; vnn = " << vnn << "; numstr = " << numstr << endl;
			}
		}

		// loop over faces, get vertex neighbors (check for redundancy)
		for (f=0; f<fn; f++){
			fe = Nv[ind][f];
			// cout << "ind = " << ind << ", face = " << f << "..." << endl;
			for (j=0; j<fe; j++){
				// get next neighbor, adjacent vertex
				vntmp1 = vfacec[f].at(j);
				vntmp2 = vfacec[f].at((j+1) % fe);

				// check for redundancy
				invec = 0;
				// cout << "invec = " << invec << "; vntmp1 = " << vntmp1 << "; vntmp2 = " << vntmp2 << "; "; 
				while (invec < vvertn[ind][vntmp1].size()){
					// cout << "vvertn[" << ind << "][" << vntmp1 << "] = " << vvertn[ind][vntmp1].at(invec) << endl;
					if (vntmp2 == vvertn[ind][vntmp1].at(invec))
						invec = -1;
					else
						invec++;
				}

				// if no redundant vertex found, add to list
				if (invec >= 0)
					vvertn[ind][vntmp1].push_back(vntmp2);
			}
		}

		while (pos2 < linestr.size()-1){
			// get locations of first dimension	
			pos1 = linestr.find_first_of("(",pos2+1);
			pos2 = linestr.find_first_of(",",pos1+1);

			// get string between pos1 and pos2 (x)
			numstr = linestr.substr(pos1+1,pos2-pos1-1);
			ss << numstr;
			ss >> val;
			ss.clear();
			vvertx[ind].push_back(val);

			// cout << "x str = " << numstr << "; ";

			// get location of second dimension
			pos1 = pos2;
			pos2 = linestr.find_first_of(",",pos1+1);

			// get string between pos1 and pos2 (y)
			numstr = linestr.substr(pos1+1,pos2-pos1-1);
			ss << numstr;
			ss >> val;
			ss.clear();
			vverty[ind].push_back(val);

			// cout << "y str = " << numstr << "; ";

			// get location of second dimension
			pos1 = pos2;
			pos2 = linestr.find_first_of(")",pos1+1);

			// get string between pos1 and pos2 (y)
			numstr = linestr.substr(pos1+1,pos2-pos1-1);
			ss << numstr;
			ss >> val;
			ss.clear();
			vvertz[ind].push_back(val);

			// cout << "z str = " << numstr << "; ";
			// cout << "pos2 = " << pos2 << ", linestr.size() = " << linestr.size() << endl;
		}
		// cout << endl;		
		// cout << linestr << endl;

		for (f=0; f<fn; f++)
			vfacec[f].clear();
		delete [] vfacec;
	}

	// get total number of vertices
	int Nvtot = 0;
	for (i=0; i<NP; i++){
		Nvtot += vvertc[i].size();
	}

	cout << "Finished reading data, plotting xyz..." << endl;

	// GET LOCATIONS OF VERTICES, UNITE ONES AT SAME POSITION, GET MASTER LIST

	// calculate vertex index k for all vertices
	vector<double> gvx;
	vector<double> gvy;
	vector<double> gvz;	
	double xtmp,ytmp,ztmp;
	int k,a;
	int* vmap;
	vmap = new int[Nvtot];

	a = 0;
	cout << "Getting vertices to merge together" << endl;
	for (i=0; i<NP; i++){		
		for (m=0; m<vvertc[i].size(); m++){
			xtmp = vvertx[i].at(m);
			ytmp = vverty[i].at(m);
			ztmp = vvertz[i].at(m);
			if (i==0 && m==0){
				gvx.push_back(xtmp);
				gvy.push_back(ytmp);
				gvz.push_back(ztmp);
				k = 0;
			}
			else
				k = merge_vertices(B,xtmp,ytmp,ztmp,gvx,gvy,gvz);

			vmap[a+m] = k;
			// cout << "p = " << i << "; site " << a+m << "; vmap[site] = " << k << endl;
		}
		a += vvertc[i].size();
	}

	// declare number of merged vertices
	int Nvmerge = gvx.size();
	vector<int>* gvn;
	gvn = new vector<int>[Nvmerge];	

	// get particles assigned to each vertex
	vector<int>* gvp;
	gvp = new vector<int>[Nvmerge];

	// loop over vertices, add unique neighbors to each vertex	
	for (k=0; k<Nvmerge; k++){
		// loop over all non-merged vertices, find k
		a = 0;
		for (i=0; i<NP; i++){		
			for (m=0; m<vvertc[i].size(); m++){
				if (vmap[a+m]==k){
					// add list of particles to vertex k
					gvp[k].push_back(i);

					// push back neighbors of (i,m) to vertex k	
					for (j=0; j<vvertn[i][m].size(); j++)
						gvn[k].push_back(vmap[a+vvertn[i][m].at(j)]);

				}
			}
			a += vvertc[i].size();
		}
	}

	// remove all redundant entries in gvn
	for (k=0; k<Nvmerge; k++){
		sort(gvn[k].begin(),gvn[k].end());
		gvn[k].erase(unique(gvn[k].begin(),gvn[k].end()),gvn[k].end());
	}

	// translate vertices to center of mass frame (cvx, cvy, cvz)
	cout << "Getting vertices center of mass..." << endl;
	double xsum, ysum, zsum, comx, comy, comz;
	xsum = 0;
	ysum = 0;
	zsum = 0;
	for (k=0; k<Nvmerge; k++){
		xsum += gvx.at(k);
		ysum += gvy.at(k);
		zsum += gvz.at(k);
	}
	comx = xsum/Nvmerge;
	comy = ysum/Nvmerge;
	comz = zsum/Nvmerge;

	// print xyz output to check that everything is correct
	
	int d,dd;
	ofstream xyzobj("vertex_edge.xyz");
	xyzobj << Nvmerge << endl;
	xyzobj << "Lattice=\"";
	for (d=0; d<NDIM; d++){
		for(dd=0; dd<NDIM; dd++){
			if (dd==d)
				xyzobj << B[d];
			else
				xyzobj << " 0.0 ";
		}
	}
	xyzobj << "\" ";
	xyzobj << '\t';
	xyzobj << "Properties=species:S:1:pos:R:" << NDIM << ":radius:R:1" << endl;
	int kmatch = 1192;
	for (k=0; k<Nvmerge; k++){
			j = 0;
			if (k == kmatch)
				xyzobj << "H";
			else{
				for (i=0; i<gvn[k].size(); i++){
					if (gvn[k].at(i) == kmatch){
						xyzobj << "N";
						j = 1;
					}
				}
				if (j == 0)
					xyzobj << "C";				
			}
			xyzobj << setw(15) << gvx.at(k);
			xyzobj << setw(15) << gvy.at(k);
			xyzobj << setw(15) << gvz.at(k);
			xyzobj << setw(15) << 0.01;
			xyzobj << endl;
	}
	

	// label all edges with reset lattice
	vector<int*> vedges;
	int l,Nedges;
	int* ptrtmp;
	ptrtmp = nullptr;
	Nedges = 0;
	for (k=0; k<Nvmerge; k++){
		for (l=0; l<gvn[k].size(); l++){
			if (gvn[k].at(l) > k){
				vedges.push_back(ptrtmp);
				vedges.at(Nedges) = new int[2];
				vedges.at(Nedges)[0] = k;
				vedges.at(Nedges)[1] = gvn[k].at(l);
				// cout << "en = " << Nedges << "; vedges.at(" << Nedges << ")[0] = " << vedges.at(Nedges)[0] << "; vedges.at(" << Nedges << ")[1] = " << vedges.at(Nedges)[1] << endl;
				Nedges++;
			}
		}
	}
	cout << endl << endl;

	// reset lattice, get ready to add nearest neighbors
	cout << "Resetting lattice with edges" << endl;
	this->reset_lattice(Nedges);

	// populate nearest neighbor info
	cout << "Getting nearest neighbor info" << endl;
	int n0,n1,nnn,e0,e1;
	vector<int> vnntmp;
	vector<int> ennn;
	for (e0=0; e0<Nedges; e0++){
		// add edges of vertex 1 to nn for e
		n0 = vedges.at(e0)[0];
		n1 = vedges.at(e0)[1];
		for (e1=0; e1<Nedges; e1++){
			if (e1 != e0){
				// if edges are connected by either n0 or n1, add to nn
				if (vedges.at(e1)[0] == n0 || vedges.at(e1)[0] == n1 || vedges.at(e1)[1] == n0 || vedges.at(e1)[1] == n1)
					vnntmp.push_back(e1);
			}
		}
		nnn = vnntmp.size();
		ennn.push_back(nnn);
		this->initialize_nn(e0,nnn);
		for (e1=0; e1<nnn; e1++)
			this->set_nn(e0,e1,vnntmp.at(e1));
		vnntmp.clear();
	}

	// get locations of points of closest approaches for all vertices

	// populate lattice
	double ex,ey,ez;	
	vector<int>* pl;
	pl = new vector<int>[Nedges];

	// vectors for plane geometry
	double v1x,v1y,v1z,v0x,v0y,v0z,dr,dx,dy,dz;
	int checkx, checky, checkz;
	double p0[NDIM];
	double p1[NDIM];
	double p2[NDIM];
	double pq[NDIM];
	double pr[NDIM];
	double nrm[NDIM];
	double lp[NDIM];
	double ln[NDIM];
	double pint[NDIM];
	double cpa[NDIM];
	double cpax[Nedges];
	double cpay[Nedges];
	double cpaz[Nedges];
	double nmag,lnmag,d1,d2,dval;


	for (e0=0; e0<Nedges; e0++){
		// cout << "e0 = " << e0 << ": ";

		// get neighbor vertices for e0
		n0 = vedges.at(e0)[0];
		n1 = vedges.at(e0)[1];

		// get vertex position
		v0x = gvx.at(n0);
		v0y = gvy.at(n0);
		v0z = gvz.at(n0);
		v1x = gvx.at(n1);
		v1y = gvy.at(n1);
		v1z = gvz.at(n1);

		// check distances, make sure edges and particles are not across boundary
		checkx = 0;
		checky = 0;
		checkz = 0;
		if (abs(v1x-v0x)>0.5*B[0]){
			if (v1x > v0x)
				v1x -= B[0];
			else
				v1x += B[0];
			checkx = 1;
		}

		if (abs(v1y-v0y)>0.5*B[1]){
			if (v1y > v0y)
				v1y -= B[1];
			else
				v1y += B[1];
			checky = 1;
		}

		if (abs(v1z-v0z)>0.5*B[2]){
			if (v1z > v0z)
				v1z -= B[2];
			else
				v1z += B[2];
			checkz = 1;
		}

		// get edge positions (with translated edges)
		ex = 0.5*(v0x+v1x);
		ey = 0.5*(v0y+v1y);
		ez = 0.5*(v0z+v1z);

		

		// get plane of 3 adjacent particle centers

		// get intersection of list of particles adjacent to vertices n0 & n1
		// l0 = new int[gvp[n0].size()];
		// for (i=0; i<gvp[n0].size(); i++)
		// 	l0[i] = gvp[n0].at(i);

		// l1 = new int[gvp[n1].size()];
		// for (i=0; i<gvp[n1].size(); i++)
		// 	l1[i] = gvp[n1].at(i);

		// get intersection of particles in gvp[n0] and gvp[n1]
		for (i=0; i<gvp[n0].size(); i++){
			for (j=0; j<gvp[n1].size(); j++){
				// if l0[i] = l1[j], then this is an intersection, so add to pl		
				if (gvp[n0].at(i)==gvp[n1].at(j)){
					pl[e0].push_back(gvp[n0].at(i));
					break;
				}
			}
		}

		// remove all redundant entries in pl
		if (pl[e0].size() > 3){
			sort(pl[e0].begin(),pl[e0].end());
			pl[e0].erase(unique(pl[e0].begin(),pl[e0].end()),pl[e0].end());
		}

		// cout << "e0 = " << e0 << " n0 = " << n0 << " n1 = " << n1 << "; plsize = " << pl[e0].size() << "; pl = ";
		// for (i=0; i<pl[e0].size(); i++)
		// 	cout << setw(10) << pl[e0].at(i);
		// cout << endl;

		// if (pl.size() > 3){
		// 	cout << "particle list > 3, ending..." << endl;
		// 	pl.clear();
		// 	throw "pl > 3";
		// }							
		

		// get particle vectors
		for (int d=0; d<NDIM; d++){
			p0[d] = pos[pl[e0].at(0)][d];				
			p1[d] = pos[pl[e0].at(1)][d];
			p2[d] = pos[pl[e0].at(2)][d];

			if (d==0){
				if (abs(p0[d]-v0x)>0.5*B[d]){
					if (p0[d]>v0x)
						p0[d]-=B[d];
					else
						p0[d]+=B[d];
				}

				if (abs(p1[d]-v0x)>0.5*B[d]){
					if (p1[d]>v0x)
						p1[d]-=B[d];
					else
						p1[d]+=B[d];
				}

				if (abs(p2[d]-v0x)>0.5*B[d]){
					if (p2[d]>v0x)
						p2[d]-=B[d];
					else
						p2[d]+=B[d];
				}
			}

			if (d==1){
				if (abs(p0[d]-v0y)>0.5*B[d]){
					if (p0[d]>v0y)
						p0[d]-=B[d];
					else
						p0[d]+=B[d];
				}

				if (abs(p1[d]-v0y)>0.5*B[d]){
					if (p1[d]>v0y)
						p1[d]-=B[d];
					else
						p1[d]+=B[d];
				}

				if (abs(p2[d]-v0y)>0.5*B[d]){
					if (p2[d]>v0y)
						p2[d]-=B[d];
					else
						p2[d]+=B[d];
				}
			}

			if (d==2){
				if (abs(p0[d]-v0z)>0.5*B[d]){
					if (p0[d]>v0z)
						p0[d]-=B[d];
					else
						p0[d]+=B[d];
				}

				if (abs(p1[d]-v0z)>0.5*B[d]){
					if (p1[d]>v0z)
						p1[d]-=B[d];
					else
						p1[d]+=B[d];
				}

				if (abs(p2[d]-v0z)>0.5*B[d]){
					if (p2[d]>v0z)
						p2[d]-=B[d];
					else
						p2[d]+=B[d];
				}
			}

			pq[d] = p1[d]-p0[d];
			pr[d] = p2[d]-p0[d];
		}

		// get plane normal vector (scale to unit length)
		nrm[0] = pq[1]*pr[2]-pq[2]*pr[1];
		nrm[1] = -pq[0]*pr[2]+pq[2]*pr[0];
		nrm[2] = pq[0]*pr[1]-pq[1]*pr[0];
		nmag = sqrt(nrm[0]*nrm[0]+nrm[1]*nrm[1]+nrm[2]*nrm[2]);
		
		nrm[0] /= nmag;
		nrm[1] /= nmag;
		nrm[2] /= nmag;

		
		// cout << "e0 = " << e0 << "; nmag = " << nmag << ", nrm[0] = " << nrm[0] << ", nrm[1] = " << nrm[1] << ", nrm[2] = " << nrm[2] << endl;
		// cout << "e0 = " << e0 << "; pq[0] = " << pq[0] << ", pq[1] = " << pq[1] << ", pq[2] = " << pq[2] << endl;
		// cout << "e0 = " << e0 << "; pr[0] = " << pr[0] << ", pr[1] = " << pr[1] << ", pr[2] = " << pr[2] << endl;

		// get midpoint of edge e0
		lp[0] = ex;
		lp[1] = ey;
		lp[2] = ez;

		// get unit vector that points along line segment
		ln[0] = gvx.at(n0) - lp[0];
		ln[1] = gvy.at(n0) - lp[1];
		ln[2] = gvz.at(n0) - lp[2];
		lnmag = sqrt(ln[0]*ln[0]+ln[1]*ln[1]+ln[2]*ln[2]);

		ln[0] /= lnmag;
		ln[1] /= lnmag;
		ln[2] /= lnmag;

		// get d, or scale of line intersecting with plane
		d2 = ln[0]*nrm[0] + ln[1]*nrm[1] + ln[2]*nrm[2];
		d1 = (p0[0]-lp[0])*nrm[0] + (p0[1]-lp[1])*nrm[1] + (p0[2]-lp[2])*nrm[2];
		dval = d1/d2;
		
		// cout << "e0 = " << e0 << "; lnmag = " << lnmag << ", d1 = " << d1 << ", d2 = " << d2 << ", dval = " << dval << endl;			
		/*
		if (1-abs(d2) > 1e-2){
			cout << "edge not perpendicular to plane, throwing error..." << endl;
			for (i=0; i<NP; i++){
				vvertx[i].clear();
				vverty[i].clear();
				vvertz[i].clear();
				vvertc[i].clear();		
				for (j=0; j<vvertc[i].size(); j++)
					vvertn[i][f].clear();
				delete [] vvertn[i];	
				delete [] Nv[i];		
			}
			for (k=0; k<Nvmerge; k++){
				gvn[k].clear();
				gvp[k].clear();
			}
			for (e0=0; e0<Nedges; e0++)
				pl[e0].clear();

			delete [] pl;
			delete [] vvertx;
			delete [] vverty;
			delete [] vvertz;
			delete [] vvertc;
			delete [] vvertn;
			delete [] Nf;
			delete [] Nv;
			delete [] vmap;
			delete [] gvn;
			delete [] gvp;

			vobj.close();
			remove(vorostr.c_str());
			throw "non-perpendicular edge-plane found!";
		}		
		*/

		// get p: intersection (or closest approach) of plane with particles in l to edge e0
		pint[0] = dval*ln[0] + lp[0];
		pint[1] = dval*ln[1] + lp[1];
		pint[2] = dval*ln[2] + lp[2];

		// if dval <= distance form midpoint to edge, then point is on edge
		if (dval <= lnmag){
			cpa[0] = pint[0];
			cpa[1] = pint[1];
			cpa[2] = pint[2];
		}
		// else, point is off edge, and no intersection but closest point of approach
		else{
			// measure distance to both vertices
			d1 = sqrt(pow(pint[0]-gvx.at(n0),2) + pow(pint[1]-gvy.at(n0),2) + pow(pint[2]-gvz.at(n0),2));
			d2 = sqrt(pow(pint[0]-gvx.at(n1),2) + pow(pint[1]-gvy.at(n1),2) + pow(pint[2]-gvz.at(n1),2));
			if (d1 > d2){
				cpa[0] = gvx.at(n1);
				cpa[1] = gvy.at(n1);
				cpa[2] = gvz.at(n1);
			}
			else{
				cpa[0] = gvx.at(n0);
				cpa[1] = gvy.at(n0);
				cpa[2] = gvz.at(n0);
			}				
		}
		cpax[e0] = cpa[0];
		cpay[e0] = cpa[1];
		cpaz[e0] = cpa[2];
	}







	int edgeoverlap = 0;

	// bisection variables
	double check,check_new,r;
	int kmax;
	k = 1;
	int kk = 1;
	kmax = 1e4;
	check = 10*epsilon;
	check_new = check;
	r = 0.5*(aH+aL);

	// perc check variables
	int cf = 1;
	double subdiv = 10.0;
	double loc,dloc;
	int perc = 0;
	string percdir = "NA";

	// output information to console
	int w = 15;
	int prc = 5;
	double poro; 
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

	cout << "Looping over a values, looking for percolation" << endl;
	while ((check > epsilon || perc == 0) && k < kmax){
	// while (check > epsilon && k < kmax){
		check = r;

		this->set_radius(r);
		this->reset_sys();
		this->reset_ptr();
		percdir = "NA";


		// CHECK INTERSECTION OF EVERY PARTICLE WITH POINT OF CLOSEST APPROACH

		// populate lattice
		// cout << endl << "Printing particle adjacencies" << endl;
		for (e0=0; e0<Nedges; e0++){
			edgeoverlap = 0;

			// // JUST CHECK ALL POSSIBLE OVERLAPS WITH ALL PARTICLES
			// dr = 0;
			// edgeoverlap = 0;
			// for (i=0; i<NP; i++){
			// 	dx = ex - pos[i][0];
			// 	dx = dx - round(dx/B[0])*B[0];

			// 	dy = ey - pos[i][1];
			// 	dy = dy - round(dy/B[1])*B[1];

			// 	dz = ez - pos[i][2];
			// 	dz = dz - round(dz/B[2])*B[2];
				
			// 	dr = sqrt(dx*dx + dy*dy + dz*dz);
			// 	if (dr < r){
			// 		edgeoverlap = 1;
			// 		break;
			// 	}
			// }

			for (kk=0; kk<pl[e0].size(); kk++){
				dx = cpax[e0] - pos[pl[e0].at(kk)][0];
				dx = dx - round(dx/B[0])*B[0];

				dy = cpay[e0] - pos[pl[e0].at(kk)][1];
				dy = dy - round(dy/B[1])*B[1];

				dz = cpaz[e0] - pos[pl[e0].at(kk)][2];
				dz = dz - round(dz/B[2])*B[2];

				dr = sqrt(dx*dx + dy*dy + dz*dz);
				if (dr < r){
					edgeoverlap = 1;
					break;
				}
			}

			// if intersects, then p is not void (lattice[e0] = 0)
			// else, p is in void (lattice[e0] = 1)
			if (edgeoverlap == 1)
				this->set_lattice(e0,0);
			else
				this->set_lattice(e0,1);
		}

		// merge clusters
		this->merge_clusters(ennn);
		this->post_process_voro();

		// check percolation in X direction (post-process)
		perc = 0;
		cf = 1;
		loc = 0;
		dloc = B[0]/subdiv;
		while(cf == 1 && loc < B[0]-dloc){
			cf = 0;

			// check for all points between loc & loc+dloc
			for (e0=0; e0<Nedges; e0++){
				// get neighbor edges
				n0 = vedges.at(e0)[0];
				n1 = vedges.at(e0)[1];

				// get edge position
				v0x = gvx.at(n0);
				v1x = gvx.at(n1);

				if (abs(v1x-v0x)>0.5*B[0]){
					if (v1x > v0x)
						v1x -= B[0];
					else
						v1x += B[0];
				}

				// get edge positions (with translated edges)
				ex = 0.5*(v0x+v1x)-comx+0.5;
				if (ex < 0)
					ex += B[0];
				else if (ex > B[0])
					ex -= B[0];

				// check for ex in x bin
				if (ex > loc && ex <= loc+dloc){
					// check if edge is occupied with biggest cluster
					if (this->findroot(e0)==this->get_pclus()){					
						cf = 1;
						break;
					}	
				}
			}

			// if max cluster found in slice, advance
			// cout << "loc = " << loc << ", cf = " << cf << endl;
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

				// cout << "loc = " << loc << endl;				

				// check for all points between loc & loc+dloc
				for (e0=0; e0<Nedges; e0++){
					// get neighbor edges
					n0 = vedges.at(e0)[0];
					n1 = vedges.at(e0)[1];

					// get edge position
					v0y = gvy.at(n0);
					v1y = gvy.at(n1);

					if (abs(v1y-v0y)>0.5*B[1]){
						if (v1y > v0y)
							v1y -= B[1];
						else
							v1y += B[1];
					}

					// get edge positions (with translated edges)
					ey = 0.5*(v0y+v1y)-comy+0.5;
					if (ey < 0)
						ey += B[1];
					else if (ey > B[1])
						ey -= B[1];

					// check for ex in x bin
					if (ey > loc && ey <= loc+dloc){
						// cout << "edge e0 = " << e0 << " in bin: ey = " << ey << ", root = " << this->findroot(e0) << ", pclus = " << this->get_pclus() << endl;
						// check if edge is occupied with biggest cluster
						if (this->findroot(e0)==this->get_pclus()){
							// cout << "so breaking..." << endl;	
							cf = 1;
							break;
						}	
					}
				}
				if (cf == 1)
					loc += dloc;
				// cout << ", cf = " << cf << "..." << endl;				
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
				for (e0=0; e0<Nedges; e0++){
					// get neighbor edges
					n0 = vedges.at(e0)[0];
					n1 = vedges.at(e0)[1];

					// get edge position
					v0z = gvz.at(n0);
					v1z = gvz.at(n1);

					if (abs(v1z-v0z)>0.5*B[2]){
						if (v1z > v0z)
							v1z -= B[2];
						else
							v1z += B[2];
					}

					// get edge positions (with translated edges)
					ez = 0.5*(v0z+v1z)-comz+0.5;
					if (ez < 0)
						ez += B[2];
					else if (ez > B[2])
						ez -= B[2];

					// check for ex in x bin
					if (ez > loc && ez <= loc+dloc){
						// check if edge is occupied with biggest cluster
						if (this->findroot(e0)==this->get_pclus()){					
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

		if (perc == 1)
			aL = r;
		else
			aH = r;

		r = 0.5*(aH+aL);
		check = abs(check-r)/r;
		check_new = check;
		k++;
	}
	// if converged, end
	if (k < kmax && check < epsilon){
		cout << "percolation found!";
		ac = r;
		vctot = poro;
		cout << "; ac = " << ac << ", poro = " << vctot << endl;
	}
	else{
		cout << "percolation not found in k < kmax :(" << endl;
		ac = -1;
		vctot = -1;
	}


	// print xyz output to check that everything is correct
	// xyzobj << endl;
	
	xyzobj << Nedges+Nvmerge << endl;
	xyzobj << "Lattice=\"";
	for (d=0; d<NDIM; d++){
		for(dd=0; dd<NDIM; dd++){
			if (dd==d)
				xyzobj << B[d];
			else
				xyzobj << " 0.0 ";
		}
	}
	xyzobj << "\" ";
	xyzobj << '\t';
	xyzobj << "Properties=species:S:1:pos:R:" << NDIM << ":radius:R:1" << endl;
	for (k=0; k<Nedges; k++){
			// get neighbor edges
			n0 = vedges.at(k)[0];
			n1 = vedges.at(k)[1];

			if (this->get_lattice_site(k)==1)
				xyzobj << setw(15) << 'O';
			else
				xyzobj << setw(15) << 'S';
			xyzobj << setw(15) << cpax[k];
			xyzobj << setw(15) << cpay[k];
			xyzobj << setw(15) << cpaz[k];
			xyzobj << setw(15) << 0.01;
			xyzobj << endl;
	}
	for (k=0; k<Nvmerge; k++){
			j = 0;
			if (k == kmatch)
				xyzobj << "H";
			else{
				for (i=0; i<gvn[k].size(); i++){
					if (gvn[k].at(i) == kmatch){
						xyzobj << "N";
						j = 1;
					}
				}
				if (j == 0)
					xyzobj << "X";				
			}
			xyzobj << setw(15) << gvx.at(k);
			xyzobj << setw(15) << gvy.at(k);
			xyzobj << setw(15) << gvz.at(k);
			xyzobj << setw(15) << 0.005;
			xyzobj << endl;
	}
	/*
	for (k=0; k<NP; k++){
			if (k==2 || k==0 || k==6)
				xyzobj << "C";
			else
				xyzobj << "Y";
			xyzobj << setw(15) << pos[k][0];
			xyzobj << setw(15) << pos[k][1];
			xyzobj << setw(15) << pos[k][2];
			xyzobj << setw(15) << ac;
			xyzobj << endl;
	}
	*/


	container_poly conp2(0,B[0],0,B[1],0,B[2],n,n,n,
			pbc,pbc,pbc,init);

	for (i=0; i<NP; i++)
		conp2.put(i,pos[i][0],pos[i][1],pos[i][2],0.001);

	// draw cells and particles
	conp2.draw_cells_pov("ov_sphere_c.pov");
	conp2.draw_particles_pov("ov_sphere_p.pov");
	

	/*

	// print vertex positions
	cout << "vertex positions: " << endl;
	for (i=0; i<NP; i++){
		cout << "vvertx[" << i << "] = ";
		for (j=0; j<vvertx[i].size(); j++)
			cout << setw(12) << vvertx[i].at(j);
		cout << endl;

		cout << "vverty[" << i << "] = ";
		for (j=0; j<vverty[i].size(); j++)
			cout << setw(12) << vverty[i].at(j);
		cout << endl;

		cout << "vvertz[" << i << "] = ";
		for (j=0; j<vvertz[i].size(); j++)
			cout << setw(12) << vvertz[i].at(j);
		cout << endl;

		cout << "vvertc[" << i << "] = ";
		for (j=0; j<vvertc[i].size(); j++)
			cout << setw(12) << vvertc[i].at(j);
		cout << endl;

		cout << "vvertn[" << i << "] = ";
		for (j=0; j<vvertc[i].size(); j++){
			cout << "( " << j << ": ";
			for (f=0; f<vvertn[i][j].size(); f++)
				cout << vvertn[i][j].at(f) << " ";
			cout << ") ";
		}
		cout << endl << endl;

	}

	cout << "Plotting vertex points: " << endl;
	for (k=0; k<Nvmerge; k++){		
		cout << "k = " << k << ": pos = " << setw(12) << gvx.at(k) << setw(12) << gvy.at(k) << setw(12) << gvz.at(k) << endl;
		cout << "	** neighbors: " << endl;
		for (j=0; j<gvn[k].size(); j++)
			cout << setw(12) << gvn[k].at(j);
		cout << endl << "	** particles: " << endl;
		for (j=0; j<gvp[k].size(); j++)
			cout << setw(12) << gvp[k].at(j);
		cout << endl;
	}

	*/

	for (i=0; i<NP; i++){
		vvertx[i].clear();
		vverty[i].clear();
		vvertz[i].clear();
		vvertc[i].clear();		
		for (j=0; j<vvertc[i].size(); j++)
			vvertn[i][f].clear();
		delete [] vvertn[i];
		delete [] Nv[i];		
	}
	for (k=0; k<Nvmerge; k++){
		gvn[k].clear();
		gvp[k].clear();
	}
	for (e0=0; e0<Nedges; e0++)
		pl[e0].clear();

	delete [] pl;
	delete [] vvertx;
	delete [] vverty;
	delete [] vvertz;
	delete [] vvertc;
	delete [] vvertn;
	delete [] Nf;
	delete [] Nv;
	delete [] vmap;
	delete [] gvn;
	delete [] gvp;

	vobj.close();
	remove(vorostr.c_str());
	//xyzobj.close();
}