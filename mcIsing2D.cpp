#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>
#include <omp.h>
#include <chrono>
#include <algorithm>

#include "util.h"
#include "mcIsing2D.h"

using namespace std;

extern int seed;
extern mt19937 rng;
extern uniform_real_distribution<double> unif;
extern int CACHE_LINE_SIZE;

int **addedI;
vector<int> toCheckI;
int sitesToCheckI;
int clusterIterI;

//Basic utilites
// x = 0,1,..., p.S1-1  and
// t = 0,1,..., p.Lt-1
// i = x + p.S1*t so
// x = i % p.S1 and
// t = i/p.S1 
inline int xp(int i, Param p) {
  //if( (i+1)%p.S1 == 0 ) return (i/p.S1)*p.S1;
  //else return i+1;
  return (i + 1) % p.S1 + p.S1 * (i / p.S1);
}

inline int xm(int i, Param p) {
  //if( i%p.S1 == 0 ) return (i/p.S1)*p.S1 + (p.S1-1);
  //else return i-1;  
  return (i - 1 + p.S1) % p.S1 + p.S1 * (i / p.S1);
}

inline int tp(int i, Param p) {
  //if(i/p.S1 == p.Lt-1) return (i - (p.Lt-1)*p.S1);
  //else return i+p.S1;
  return i % p.S1 + p.S1 * ((i / p.S1 + 1) % p.Lt);
}

inline int ttm(int i, Param p){
  //if(i/p.S1 == 0) return (i + (p.Lt-1)*p.S1);
  //else return i-p.S1;
  return i % p.S1 + p.S1 * ((i / p.S1 - 1 + p.Lt) % p.Lt);
}

//----------------------//
// Short Range Routines //
//----------------------//

// declare variables to implement Cluster algorithms
//Declare variables to implement Cluster algorithms
int ising_wc_ave = 0;
int ising_wc_size = 0;
int ising_wc_calls = 0;
int ising_wc_poss = 0;
int ising_wc_t_size = 0;
int ising_wc_s_size = 0;

int ising_sw_ave = 0;
int ising_sw_size = 0;
int ising_sw_calls = 0;
int ising_sw_t_size = 0;
int ising_sw_s_size = 0;

int ising_accept = 0;
int ising_tries  = 0;

void metropolisUpdateISR(int *s, Param p, int iter) {

  double DeltaE = 0.0;
  double extern_h = p.h;
  double J = p.J;
  
  for (int i = 0; i < p.surfaceVol; i++) {

    //External field
    DeltaE += -extern_h*s[i];
    
    //Interaction
    DeltaE += -J*(s[xp(i, p)] + s[xm(i, p)] +
		  s[tp(i, p)] + s[ttm(i, p)]);
    
    ising_tries++;
    
    if(DeltaE < 0.0) {
      //  cout<< " Acepted  " << endl;
      ising_accept += 1;
      s[i] *= -1;
    }
    else if ( unif(rng)  < exp(-DeltaE)) {
      //  cout<< " Acepted  " << endl;
      ising_accept += 1;
      s[i] *= -1;
    }     
  }// end loop over lattice volume

  if( iter < p.n_therm && iter%p.n_skip == 0 ) {
    cout<<"At iter "<<iter<<" the Metro acceptance rate is "<<(double)ising_accept/(double)ising_tries<<endl;
  }
  
}

void swendsenWangUpdateISR(int *s, Param p, int iter) {
  

}

void swendsenWangClusterAddISR(int i, int *s, int cSpin, int clusterNum, 
			       int *clusterDef, Param p) {
  
}

void wolffUpdateISR(int *s, Param p, int iter) {

  ising_wc_calls++;
  
  bool *cluster = new bool[p.surfaceVol];
  for (int i = 0; i < p.surfaceVol; i++)
    cluster[i] = false;

  // choose a random spin and grow a cluster
  int i = int(unif(rng) * p.surfaceVol);
  int cSpin = s[i];
  
  //This function is recursive and will call itself
  //until all four attempts in the lattice directions
  // (+x, -x, +t, -t) have failed to ncrease the cluster.
  ising_wc_size = 1;
  wolffClusterAddISR(i, s, cSpin, cluster, p);

  ising_wc_ave += ising_wc_size;

  if( iter%p.n_skip == 0) {
    setprecision(4);
    cout<<"Ave. cluster size at iter "<<iter<<" = "<<ising_wc_ave<<"/"<<ising_wc_calls<<" = "<<1.0*ising_wc_ave/ising_wc_calls<<endl;
  }
  delete[] cluster;
  
}

void wolffClusterAddISR(int i, int *s, int cSpin,
			bool *cluster, Param p) {

  double J = p.J;
  
  // The site belongs to the cluster, so flip it.
  cluster[i] = true;
  s[i] *= -1;

  // - If the (aligned) neighbour spin does not already belong to the
  // cluster, then try to add it to the cluster.
  // - If the site has already been added, then we may skip the test.

  //Forward in T
  if(!cluster[ tp(i,p) ] && s[tp(i,p)] == cSpin) {
    if(unif(rng) < 1 - exp(-2*J)) {
      ising_wc_size++;
      ising_wc_t_size++;
      //cout<<"->tp";
      wolffClusterAddISR(tp(i,p), s, cSpin, cluster, p);
    }
  }

  //Forward in X
  if(!cluster[ xp(i,p) ] && s[xp(i,p)] == cSpin) {
    if(unif(rng) < 1 - exp(-2*J)) {
      ising_wc_size++;
      ising_wc_s_size++;
      //cout<<"->xp";
      wolffClusterAddISR(xp(i,p), s, cSpin, cluster, p);
    }
  }

  
  //Backard in T 
  if(!cluster[ ttm(i,p) ] && s[ttm(i,p)] == cSpin) {  
    if(unif(rng) < 1 - exp(-2*J)) {
      ising_wc_size++;
      ising_wc_t_size++;
      //cout<<"->tm";
      wolffClusterAddISR(ttm(i,p), s, cSpin, cluster, p);
    }
  }

  //Backward in X
  if(!cluster[ xm(i,p) ] && s[xm(i,p)] == cSpin) {
    if (unif(rng) < 1 - exp(-2*J)) {
      ising_wc_size++;
      ising_wc_s_size++;
      //cout<<"->xm";
      wolffClusterAddISR(xm(i,p), s, cSpin, cluster, p);
    }
  }   
}


//---------------------//
// Long Range Routines //
//---------------------//

void init_connectivityI(Param p) {
  
  addedI = (int**)malloc(p.surfaceVol*sizeof(int*));  
  for(int i=0; i<p.surfaceVol; i++) {
    addedI[i] = (int*)malloc((p.surfaceVol+1)*sizeof(int));
#ifdef USE_OMP
#pragma omp parallel for
#endif
    for(int j=0; j<p.surfaceVol+1; j++) {
      addedI[i][j] = -1;
    }
  }
}



void Ising2D::metropolisUpdateILR(Param p, int iter) {
  
  
}

void Ising2D::swendsenWangUpdateILR(Param p, int iter) {
  
}

void swendsenWangClusterAddILR(int i, int *s, int cSpin, int clusterNum, 
			       int *clusterDef, vector<int> Rcluster, 
			       double *LR_couplings, Param p) {
  
}

void clusterPossibleILR(int i, int *s, int cSpin,
		       bool *Pcluster, vector<int> &Rcluster,
		       Param p) {


}

void wolffUpdateILR(int *s, Param p, double *LR_couplings, int iter) {

  ising_wc_calls++;

  if(ising_wc_calls == p.n_therm) {
    ising_wc_calls = 1;
    ising_wc_ave = 0;
  }
  
  //Choose a random spin.
  int i = int(unif(rng) * p.surfaceVol);
  int cSpin = s[i];
  
  // The site belongs to the cluster, so flip it.
  ising_wc_size = 1;
  s[i] *= -1;
  
  //This function is recursive and will call itself
  //until all attempts to increase the cluster size
  //have failed.
  wolffClusterAddILR(i, s, cSpin, LR_couplings, p);
  
  ising_wc_ave += ising_wc_size;

  if(iter%p.n_skip == 0) {
    setprecision(4);
    cout<<"Average (CPU) cluster size at iter "<<iter<<" = "<<ising_wc_ave<<"/"<<ising_wc_calls<<" = "<<1.0*ising_wc_ave/ising_wc_calls<<" = "<<100.0*ising_wc_ave/(ising_wc_calls*p.surfaceVol)<<"%"<<endl;
  }
}

void wolffClusterAddILR(int i, int *s, int cSpin, double *LR_couplings, Param p) {

  int S1 = p.S1;
  int x_len = S1/2 + 1;
  double J = p.J;
  
  int t1,x1,t2,x2,dt,dx;
  t1 = i / p.S1;
  x1 = i % p.S1;
  
  double prob = 0.0;
  double rand = 0.0;
  //We now loop over the possible lattice sites, adding sites
  //(creating bonds) with the specified LR probablity.
  for(int j=0; j<p.surfaceVol; j++) {
    if(s[j] == cSpin && j != i) {
      
      //Index divided by circumference, using the int floor feature/bug,
      //gives the timeslice index.
      t2 = j / p.S1;
      dt = abs(t2-t1) > p.Lt/2 ? p.Lt - abs(t2-t1) : abs(t2-t1);
      
      //The index modulo the circumference gives the spatial index.
      x2 = j % p.S1;            
      dx = abs(x2-x1) > p.S1/2 ? p.S1 - abs(x2-x1) : abs(x2-x1);      
      
      prob = 1 - exp(-2*J*LR_couplings[dx + dt*x_len]);
      rand = unif(rng);
      if(rand < prob) {
	ising_wc_size++;
	// The site belongs to the cluster, so flip it.
	s[j] *= -1;
	wolffClusterAddILR(j, s, cSpin, LR_couplings, p);
      }
    }
  }

  
}

double energyISR(int *s, Param p, double &KE) {  

  KE = 0.0;
  
  int S1 = p.S1;
  int Lt = p.Lt;
  int vol = S1*Lt;  
  int s_lc = 1;
  
  for (int i = 0; i < vol; i++) {    
    s_lc = s[i];    
    KE += s_lc*(s[xp(i,p)] + s[xm(i,p)] + s[tp(i,p)] + s[ttm(i,p)]);  
  }
  return -0.5*p.J*KE;
}

double energyILR(int *s, Param p, double *LR_couplings, double &KE) {  

  KE = 0.0;
  int S1 = p.S1;
  int Lt = p.Lt;
  int vol = S1*Lt;  
  int x_len = S1/2 + 1;
  int s_lc = 1;
  
  for (int i = 0; i < vol; i++) {    
    s_lc = s[i];

    int t1,x1;
    t1 = i / S1;
    x1 = i % S1;

    //KE terms
#ifdef USE_OMP
    
    int chunk = CACHE_LINE_SIZE/sizeof(double);
    
#pragma omp parallel 
    {
      double local = 0.0;
      double val = 0.0;
      int t2,dt,x2,dx;
#pragma omp for nowait 
      for(int j=0; j<vol; j+=chunk) {
	for(int k=0; k<chunk; k++) { 

	  if( (j+k) != i ) {
	    //Index divided by circumference, using the int floor feature/bug,
	    //gives the timeslice index.
	    t2 = (j+k) / S1;
	    dt = abs(t2-t1) > Lt/2 ? Lt - abs(t2-t1) : abs(t2-t1);
	    
	    //The index modulo the circumference gives the spatial index.
	    x2 = (j+k) % S1;            
	    dx = abs(x2-x1) > S1/2 ? S1-abs(x2-x1) : abs(x2-x1);      
	    
	    val = s_lc*s[j]*LR_couplings[dx+dt*x_len];	    
	    local += val;
	  }
	}
      }
#pragma omp atomic 
      KE += local;
    }    
#else

    int t2,dt,x2,dx;
    for(int j=0; j<p.surfaceVol; j++) {
      
      if( j != i ) {
	//Index divided by circumference, using the int floor feature/bug,
	//gives the timeslice index.
	t2 = j / p.S1;
	dt = abs(t2-t1) > p.Lt/2 ? p.Lt - abs(t2-t1) : abs(t2-t1);
	
	//The index modulo the circumference gives the spatial index.
	x2 = j % p.S1;            
	dx = abs(x2-x1) > p.S1/2 ? p.S1-abs(x2-x1) : abs(x2-x1);      

	KE += s_lc*s[j]*LR_couplings[dx+dt*x_len];
	
      }
    }
#endif
  }
  
  return -0.5*p.J*KE;
}











#if 0

  //This implementation attempts to parallelises only over those sites which 
  //are candidate sites. This has the advantage that the number of computations 
  //reduce in number as the cluster is filled.

  int newSites = 0;
  int lc_sze = toCheckI.size()/omp_get_num_threads();
  int LC_arr_len = p.surfaceVol/2+1;

  int t1,x1;
  t1 = i / p.S1;
  x1 = i % p.S1;

  //We now loop over the possible lattice sites, adding sites
  //(creating bonds) with the specified LR probablity.
#pragma omp parallel 
  {
    int sitesRemoved_local = 0;
    int newSites_local = 0;
    double prob = 0.0;
    double rand = 0.0;
    int t2,dt,x2,dx;
    int added_local[lc_sze][2];
    for(int k=0; k<lc_sze; k++) {
      added_local[k][0] = -1;
      added_local[k][1] = 0;
    }
    
#pragma omp for nowait 
    for(int j=0; j<toCheckI.size(); j+=lc_sze) {
      for(int k=0; k<lc_sze; k++) { 
	int idx = toCheckI[j+k];

	//cout<<"Checking: j="<<j<<" k="<<k<<" size="<<toCheckI.size()<<" toCheckI="<<toCheckI[j+k]<<endl;
	
	//Index divided by circumference, using the int floor feature/bug,
	//gives the timeslice index.
	t2 = (idx) / p.S1;
	dt = abs(t2-t1) > p.Lt/2 ? p.Lt - abs(t2-t1) : abs(t2-t1);
	
	//The index modulo the circumference gives the spatial index.
	x2 = (idx) % p.S1;            
	dx = abs(x2-x1) > p.S1/2 ? p.S1 - abs(x2-x1) : abs(x2-x1);      
	
	prob = 1 - exp(2*phi_lc*phi[idx]*LR_couplings[dt+dx*LC_arr_len]);
	rand = unif(rng);
	//cout<<rand<<" "<<prob<<endl;
	if(rand < prob) {
	  sqnl_wc_size++;
	  cout<<"Hit: j="<<j<<" k="<<k<<" size="<<toCheckI.size()<<" toCheckI="<<toCheckI[j+k]<<endl;
	  //toSeed.push_back(idx);
	  //++sitesAdded_local;
	  added_local[newSites_local][0] = idx;
	  added_local[newSites_local][1] = j+k;
	  ++newSites_local;
	}
      }
    }

#pragma omp critical
    newSites += newSites_local;
#pragma omp for
    for(int k=0; k<lc_sze; k++) {
      if(added_local[k][0] > 0) {
	cout<<"Adding site "<<added_local[k][0]<<endl;
	addedI[i][added_local[k][0]] = 1;
	s[added_local[k][0]] *= -1;
	phi[added_local[k][0]] *= -1;
	cout<<"Added site "<<added_local[k][0]<<endl;
       	
      }
    }
  }
  
  sitesToCheckI += newSites;
  
  if(newSites > 0) {
    cout<<"SITES TO ADD!"<<endl;
    //'added' now contains all of the spins on the wave front. Each
    //of these new sites must be explored, but no new exploration 
    //need test these sites. First, sort the added array so that all
    //the hits are at the beginning, order is unimportant.
    int pos = 0;
    int l = 0;
    for(l=0; l<p.surfaceVol; l++) {
      if(addedI[i][l] > 0) {
	cout<<"ADDING "<<l<<endl;
	addedI[i][pos] = l;
	pos++;

	cout<<"REDUCE THE SEARCH SPACE"<<endl;
	cout<<"Removing element "<<added_local[k][1]<<" (idx="<<added_local[k][0]<<")"<<endl;
	toCheckI.erase (toCheckI.begin() + added_local[k][1]);
	cout<<"SEARCH SPACE REDUCED"<<endl;

      }
    }
    cout<<"ADDED!"<<endl;
  }

  addedI[i][p.surfaceVol] = newSites;
  
  cout<<"clusterIterI: "<<clusterIterI<<" SITES TO CHECK: "<<sitesToCheckI;
  cout<<" Checked site: "<<i<<" New Sites: "<<newSites;
  
  for(int k=0; k<addedI[i][p.surfaceVol]; k++) {
    cout<<" Going to "<<addedI[i][k]<<" k="<<k<<" loop max="<<addedI[i][p.surfaceVol]<<endl;
    sitesToCheckI -= 1;
    wolffClusterAddLR(addedI[i][k], s, cSpin, LR_couplings, phi, p);
  }

#ifdef USE_OMP
  //update the toCheckI array with all the sites of the same sign.
  toCheckI.push_back(i);
  for(int j=0; j<p.surfaceVol; j++) {
    if(phi[j]*cSpin > 0 && j != i) {
      toCheckI.push_back(j);
    }
  }
  //cout<<endl<<"START: "<<toCheckI.size()<<endl;
  for(int j=0; j<p.surfaceVol; j++) {
    //if(phi[j]*cSpin > 0) cout<<j<<":"<<phi[j]<<" ";
  }  
#endif

#endif
