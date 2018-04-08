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
#include "data_proc.h"
#include "mc_util.h"
#include "mc_update_sr.h"

using namespace std;

extern int seed;
extern mt19937_64 rng;
extern uniform_real_distribution<double> unif;

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

void metropolisUpdateSR(double *phi_arr, int *s, Param p, int iter) {

  int s_old = 0;  
  double phi_new = 0.0;
  double phi_new_sq = 0.0;
  double phi = 0.0;
  double phi_sq = 0.0;
  double lambda_p = 0.25*p.lambda;
  double musqr_p  = 0.50*p.musqr;
  
  double DeltaE = 0.0;
  
  for (int i = 0; i < p.surfaceVol; i++) {

    //Set some values we use a lot
    phi = phi_arr[i];
    DeltaE = 0.0;
    phi_new = phi + p.delta_phi * (2.0*unif(rng) - 1.0);
    phi_new_sq = phi_new*phi_new;
    phi_sq = phi*phi;
    
    if (s[i] * phi < 0) {
      printf("ERROR s and phi NOT aligned! (MUP)\n");
      exit(0);
    }

    //PE
    DeltaE += lambda_p*(phi_new_sq*phi_new_sq - phi_sq*phi_sq);
    DeltaE += musqr_p *(phi_new_sq            - phi_sq);
    
    //KE
    DeltaE += 2.0 * (phi_new_sq-phi_sq);
    DeltaE +=       (phi-phi_new)*(phi_arr[xp(i, p)] + phi_arr[xm(i, p)] +
				   phi_arr[tp(i, p)] + phi_arr[ttm(i, p)]);
    
    sql_tries++;
    
    if(DeltaE < 0.0) {
      //  cout<< " Acepted  " << endl;
      s_old = s[i];
      phi_arr[i] = phi_new;
      sql_accept += 1;
      s[i] = (phi_new > 0) ? 1 : -1;
    }
    else if ( unif(rng)  < exp(-DeltaE)) {
      //  cout<< " Acepted  " << endl;
      s_old = s[i];
      phi_arr[i] = phi_new;
      sql_accept += 1;
      s[i] = (phi_new > 0) ? 1 : -1;
    }     
  }// end loop over lattice volume 

  /*
  // TUNING ACCEPTANCE 
  if (iter < p.n_therm/2 && (iter+1) % p.n_skip/10 == 0) {
    if ((double) sql_accept / (double) sql_tries < 0.5) {
      p.delta_phi -= 0.001;
    } else {
      p.delta_phi += 0.001;
    }
    if(p.n_cluster*1.0*sql_wc_ave/sql_wc_calls < p.surfaceVol && iter > p.n_skip) {
      p.n_cluster++;
    } else {
      p.n_cluster--;
      if(p.n_cluster < 3) p.n_cluster++;
    }
  }
  */

  if( iter < p.n_therm && iter%p.n_skip == 0 ) {
    cout<<"At iter "<<iter<<" the Metro acceptance rate is "<<(double)sql_accept/(double)sql_tries<<endl;
  }
}

double actionSR(double *phi_arr, int *s, Param p,
		double & KE, double & PE) {  
  
  KE = 0.0;
  PE = 0.0;
  double phi_sq;
  double phi;
  double lambda_p = 0.25*p.lambda;
  double musqr_p  = 0.50*p.musqr;

  
  for (int i = 0; i < p.surfaceVol; i++)
    if (s[i] * phi_arr[i] < 0)
      printf("ERROR s and phi NOT aligned (actionPhi Square) ! \n");
  
  //PE terms
#ifdef USE_OMP
#pragma omp parallel for
#endif
  for (int i = 0; i < p.surfaceVol; i++) {    
    phi = phi_arr[i];
    phi_sq = phi*phi;
    
    PE += lambda_p * phi_sq*phi_sq;
    PE += musqr_p  * phi_sq;
    
    KE += 0.5 * (phi - phi_arr[xp(i,p)]) * (phi - phi_arr[xp(i,p)]);
    KE += 0.5 * (phi - phi_arr[tp(i,p)]) * (phi - phi_arr[tp(i,p)]);
  }
  
  return PE + KE;
}

void swendsenWangUpdateSR(double *phi_arr, int *s, Param p, int iter) {
  
  sql_sw_calls++;  
  int clusterNum = 0;
  
  //Integer array holding the cluster number each site
  //belongs to. If zero, the site is unchecked.
  int *clusterDef = new int[p.surfaceVol];
  for (int i = 0; i < p.surfaceVol; i++)
    clusterDef[i] = 0;
  
  //Integer array holding the spin value of each cluster,
  int *clusterSpin = new int[p.surfaceVol];
  for (int i = 0; i < p.surfaceVol; i++)
    clusterSpin[i] = 1;

  for (int i = 0; i < p.surfaceVol; i++) {
    if(clusterDef[i] == 0) {
      //This is the start of a new cluster.
      clusterNum++; 
      clusterDef[i] = clusterNum;
      s[i] < 0 ? clusterSpin[clusterNum] = -1 : clusterSpin[clusterNum] = 1;
      
      //This function will call itself recursively until it fails to 
      //add to the cluster
      swendsenWangClusterAddSR(i, s, clusterSpin[clusterNum], clusterNum, 
				clusterDef, phi_arr, p);
    }
  }
  
  //Loop over the defined clusters, flip with probabilty 0.5.
  for(int i=1; i<=clusterNum; i++) {
    if(unif(rng) < 0.5) clusterSpin[i] = -1;
    else clusterSpin[i] = 1;
  }
  
  //Apply spin flip. If a site is not in a cluster, its cluster value
  //is zero and its spin is not flipped.
  for (int i = 0; i < p.surfaceVol; i++) {
    phi_arr[i] *= clusterSpin[clusterDef[i]];
    s[i]       *= clusterSpin[clusterDef[i]];
  }      

  delete[] clusterDef;
  delete[] clusterSpin;
}

void swendsenWangClusterAddSR(int i, int *s, int cSpin, int clusterNum, 
			       int *clusterDef, double *phi_arr, Param p) {
  
  //The site belongs to the (clusterNum)th cluster
  clusterDef[i] = clusterNum;

  //If the (aligned) neighbor spin does not already belong to the
  // cluster, then try to add it to the cluster.
  // - If the site has already been added, then we may skip the test.

  //Forward in T
  if(clusterDef[tp(i,p)] == 0 && s[tp(i,p)] == cSpin) {
    if(unif(rng) < 1 - exp(-2*phi_arr[i]*phi_arr[tp(i,p)])) {
      //cout<<"->tp";
      swendsenWangClusterAddSR(tp(i,p), s, cSpin, clusterNum, 
				clusterDef, phi_arr, p);
    }
  }

  //Forward in X
  if(clusterDef[xp(i,p)] == 0 && s[xp(i,p)] == cSpin) {
    if(unif(rng) < 1 - exp(-2*phi_arr[i]*phi_arr[xp(i,p)])) {
      //cout<<"->xp";
      swendsenWangClusterAddSR(xp(i,p), s, cSpin, clusterNum, 
				clusterDef, phi_arr, p);
    }
  }
  
  //Backard in T 
  if(clusterDef[ttm(i,p)] == 0 && s[ttm(i,p)] == cSpin) {  
    if(unif(rng) < 1 - exp(-2*phi_arr[i]*phi_arr[ttm(i,p)])) {
      //cout<<"->tm";
      swendsenWangClusterAddSR(ttm(i,p), s, cSpin, clusterNum, 
				clusterDef, phi_arr, p);
    }
  }
  
  //Backward in X
  if(clusterDef[xm(i,p)] == 0 && s[xm(i,p)] == cSpin) {
    if (unif(rng) < 1 - exp(-2*phi_arr[i]*phi_arr[xm(i,p)])) {
      //cout<<"->xm";
      swendsenWangClusterAddSR(xm(i,p), s, cSpin, clusterNum, 
				clusterDef, phi_arr, p);
    }
  } 
}

void wolffUpdateSR(double *phi_arr, int *s, Param p, int iter) {
  
  sql_wc_calls++;

  bool *cluster = new bool[p.surfaceVol];
  for (int i = 0; i < p.surfaceVol; i++)
    cluster[i] = false;

  // choose a random spin and grow a cluster
  int i = int(unif(rng) * p.surfaceVol);
  int cSpin = s[i];
  
  //This function is recursive and will call itself
  //until all four attempts in the lattice directions
  // (+x, -x, +t, -t) have failed to ncrease the cluster.
  sql_wc_size = 1;
  wolffClusterAddSR(i, s, cSpin, cluster, phi_arr, p);

  sql_wc_ave += sql_wc_size;

  if( iter%p.n_skip == 0) {
    setprecision(4);
    cout<<"Using "<<p.n_cluster<<" Wolff hits."<<endl; 
    cout<<"Ave. cluster size at iter "<<iter<<" = "<<sql_wc_ave<<"/"<<sql_wc_calls<<" = "<<1.0*sql_wc_ave/sql_wc_calls<<endl;
  }
  delete[] cluster;
}

void wolffClusterAddSR(int i, int *s, int cSpin,
			bool *cluster, double *phi_arr, Param p) {
  
  // The site belongs to the cluster, so flip it.
  cluster[i] = true;
  s[i] *= -1;
  phi_arr[i] *= -1;

  // - If the (aligned) neighbor spin does not already belong to the
  // cluster, then try to add it to the cluster.
  // - If the site has already been added, then we may skip the test.

  //Forward in T
  if(!cluster[ tp(i,p) ] && s[tp(i,p)] == cSpin) {
    if(unif(rng) < 1 - exp(2*phi_arr[i]*phi_arr[tp(i,p)])) {
      sql_wc_size++;
      sql_wc_t_size++;
      //cout<<"->tp";
      wolffClusterAddSR(tp(i,p), s, cSpin, cluster, phi_arr, p);
    }
  }

  //Forward in X
  if(!cluster[ xp(i,p) ] && s[xp(i,p)] == cSpin) {
    if(unif(rng) < 1 - exp(2*phi_arr[i]*phi_arr[xp(i,p)])) {
      sql_wc_size++;
      sql_wc_s_size++;
      //cout<<"->xp";
      wolffClusterAddSR(xp(i,p), s, cSpin, cluster, phi_arr, p);
    }
  }

  
  //Backard in T 
  if(!cluster[ ttm(i,p) ] && s[ttm(i,p)] == cSpin) {  
    if(unif(rng) < 1 - exp(2*phi_arr[i]*phi_arr[ttm(i,p)])) {
      sql_wc_size++;
      sql_wc_t_size++;
      //cout<<"->tm";
      wolffClusterAddSR(ttm(i,p), s, cSpin, cluster, phi_arr, p);
    }
  }

  //Backward in X
  if(!cluster[ xm(i,p) ] && s[xm(i,p)] == cSpin) {
    if (unif(rng) < 1 - exp(2*phi_arr[i]*phi_arr[xm(i,p)])) {
      sql_wc_size++;
      sql_wc_s_size++;
      //cout<<"->xm";
      wolffClusterAddSR(xm(i,p), s, cSpin, cluster, phi_arr, p);
    }
  } 
}

