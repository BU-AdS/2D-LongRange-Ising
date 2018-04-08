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

#define CACHE_LINE_SIZE sysconf(_SC_LEVEL1_DCACHE_LINESIZE)

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
int sql_wc_ave = 0;
int sql_wc_size = 0;
int sql_wc_calls = 0;
int sql_wc_t_size = 0;
int sql_wc_s_size = 0;

int sql_sw_ave = 0;
int sql_sw_size = 0;
int sql_sw_calls = 0;
int sql_sw_t_size = 0;
int sql_sw_s_size = 0;

int sql_accept = 0;
int sql_tries  = 0;

void metropolisUpdateISR(int *s, Param p, int iter) {

  int s_old = 0;  
  double DeltaE = 0.0;
  
}

void swendsenWangUpdateISR(double *phi_arr, int *s, Param p, int iter) {
  

}

void swendsenWangClusterAddISR(int i, int *s, int cSpin, int clusterNum, 
			       int *clusterDef, double *phi_arr, Param p) {
  
}

void wolffUpdateISR(double *phi_arr, int *s, Param p, int iter) {
  
}

void wolffClusterAddISR(int i, int *s, int cSpin,
			bool *cluster, double *phi_arr, Param p) {
  
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

//Declare variables to implement Cluster algorithms
int sqnl_wc_ave = 0;
int sqnl_wc_size = 0;
int sqnl_wc_calls = 0;
int sqnl_wc_poss = 0;
int sqnl_wc_t_size = 0;
int sqnl_wc_s_size = 0;

int sqnl_sw_ave = 0;
int sqnl_sw_size = 0;
int sqnl_sw_calls = 0;
int sqnl_sw_t_size = 0;
int sqnl_sw_s_size = 0;

int sqnl_accept = 0;
int sqnl_tries  = 0;

void Ising2D::metropolisUpdateILR(Param p, int iter) {
  

}

void Ising2D::swendsenWangUpdateILR(Param p, int iter) {
  
}

void swendsenWangClusterAddILR(int i, int *s, int cSpin, int clusterNum, 
				int *clusterDef, vector<int> Rcluster, 
				double *LR_couplings,
				double *phi_arr, Param p) {
  
}

void clusterPossibleILR(int i, int *s, int cSpin,
		       bool *Pcluster, vector<int> &Rcluster,
		       Param p) {


}

void wolffUpdateILR(double *phi, int *s, Param p,
		   double *LR_couplings, int iter, int n) {
  

}

void wolffClusterAddILR(int i, int *s, int cSpin, double *LR_couplings,
		       double *phi, Param p) {
  
}

double actionILR(double *phi_arr, int *s, Param p,
		double *LR_couplings, 
		double &KE, double &PE) {  
  
  return PE + KE;
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
