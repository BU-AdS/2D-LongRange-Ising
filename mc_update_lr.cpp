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
#include "mc_update_lr.h"

using namespace std;

extern int seed;
extern mt19937 rng;
extern uniform_real_distribution<double> unif;

#define CACHE_LINE_SIZE sysconf(_SC_LEVEL1_DCACHE_LINESIZE)

int **added;

void init_connectivity(Param p) {

  added = (int**)malloc(p.surfaceVol*sizeof(int*));  
  for(int i=0; i<p.surfaceVol; i++) {
    added[i] = (int*)malloc(p.surfaceVol*sizeof(int));
#ifdef USE_OMP
#pragma omp parallel for
#endif
    for(int j=0; j<p.surfaceVol; j++) {
      added[i][j] = -1;
    }
  }
}
  

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

int metropolisUpdateLR(double *phi_arr, int *s, Param &p,
		       double *LR_couplings,
		       double &delta_mag_phi, int iter) {

  delta_mag_phi = 0.0;
  
  int s_old     = 0;
  int delta_mag = 0;
  
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
      printf("ERROR s and phi NOT aligned! (MetroUpdate LR)\n");
      cout<<"S["<<i<<"] = "<<s[i]<<" and phi["<<i<<"] = "<<phi<<endl;
      exit(0);
    }
    
    //PE
    DeltaE += lambda_p*(phi_new_sq*phi_new_sq - phi_sq*phi_sq);
    DeltaE += musqr_p *(phi_new_sq            - phi_sq);
    
    //KE
    double pmpn = phi-phi_new;
    double pnsmps = 0.5*phi_new_sq-phi_sq;
    
#ifdef USE_OMP
    int chunk = p.surfaceVol/omp_get_num_threads();//;CACHE_LINE_SIZE/sizeof(double);
#pragma omp parallel for reduction(+:DeltaE)
    for(int j=0; j<p.surfaceVol; j+=chunk) {
      for(int k=0; k<chunk; k++) {
	int idx = j+k;
	//cout<<omp_get_num_threads()<<" "<<omp_get_thread_num()<<" "<<idx<<endl;
	//double val = (pmpn*phi_arr[j*omp_nt + k] + pnsmps)*LR_couplings[i+(j*omp_nt + k)*p.surfaceVol]*LR_couplings[i+(j*omp_nt + k)*p.surfaceVol];
	//val = (pmpn*phi_arr[j+k] + pnsmps)*LR_couplings[j+k + i*p.surfaceVol]*LR_couplings[j+k + i*p.surfaceVol];
	DeltaE += (pmpn*phi_arr[idx] + pnsmps)*LR_couplings[idx + i*p.surfaceVol]*LR_couplings[idx + i*p.surfaceVol];
      }
    }
#elif
    for(int j=0; j<p.surfaceVol; j++) {
      DeltaE += (pmpn*phi_arr[j] + pnsmps)*LR_couplings[j+i*p.surfaceVol]*LR_couplings[j+i*p.surfaceVol];
    }
#endif

    sqnl_tries++;
    
    if(DeltaE < 0.0) {
      //  cout<< " Acepted  " << endl;
      s_old = s[i];
      delta_mag_phi += phi_new - phi_arr[i];
      phi_arr[i] = phi_new;
      sqnl_accept += 1;
      s[i] = (phi_new > 0) ? 1 : -1;
      delta_mag += s[i] - s_old;
    }
    else if ( unif(rng) < exp(-DeltaE)) {
      //  cout<< " Acepted  " << endl;
      s_old = s[i];
      delta_mag_phi += phi_new - phi_arr[i];
      phi_arr[i] = phi_new;
      sqnl_accept += 1;
      s[i] = (phi_new > 0) ? 1 : -1;
      delta_mag += s[i] - s_old;
    }     
  }// end loop over lattice volume 

  if( iter < p.n_therm && iter%p.n_skip == 0 ) {
    cout<<"At iter "<<iter<<" the Metro acceptance rate is "<<(double)sqnl_accept/(double)sqnl_tries<<endl;
    cout<<"and delta_phi is "<<p.delta_phi<<endl;
  }

  return delta_mag;
}

void swendsenWangUpdateLR(double *phi_arr, int *s, Param p,
			    double *LR_couplings, 
			    double &delta_mag_phi, int iter) {

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

  //Tracks which sites are potentially in the cluster.
  bool *Pcluster = new bool[p.surfaceVol];

  //Records which sites are potentially in the cluster.
  //This will have a modifiable size, hence the vector
  //style.
  vector<int> Rcluster;

  for (int i = 0; i < p.surfaceVol; i++) {
    if(clusterDef[i] == 0) {
      //This is the start of a new cluster.
      clusterNum++; 
      clusterDef[i] = clusterNum;
      s[i] < 0 ? clusterSpin[clusterNum] = -1 : clusterSpin[clusterNum] = 1;

      //First, we must identify the maximum possible cluster, then we 
      //can loop over only those sites 

      for (int i = 0; i < p.surfaceVol; i++) Pcluster[i] = false;
      
      sqnl_wc_poss = 0;
      clusterPossibleLR(i, s, clusterSpin[clusterNum], 
				      Pcluster, Rcluster, p);
      
      //This function will call itself recursively until it fails to 
      //add to the cluster
      swendsenWangClusterAddLR(i, s, clusterSpin[clusterNum], clusterNum, 
				 clusterDef, Rcluster, LR_couplings, 
				 phi_arr, p);

      Rcluster.clear();
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

  if( iter%p.n_skip == 0 && iter < p.n_therm) {
    setprecision(4);
    cout<<"Using "<<p.n_cluster<<" SW hits."<<endl; 
    cout<<"Ave. number of clusters at iter "<<iter<<" = "<<clusterNum<<endl;
  }


  delete[] clusterDef;
  delete[] clusterSpin;
  delete[] Pcluster;
}

void swendsenWangClusterAddLR(int i, int *s, int cSpin, int clusterNum, 
				int *clusterDef, vector<int> Rcluster, 
				double *LR_couplings,
				double *phi_arr, Param p) {
  
  double probsw = 0.0;
  double randsw = 0.0;
  int idxj = 0;
  int idxij = 0;

  //The site belongs to the (clusterNum)th cluster
  clusterDef[i] = clusterNum;

  //We now loop over the possible lattice sites, adding sites
  //(creating bonds) with the specified LR probablity.
  for(int j=0; j<sqnl_wc_poss; j++) {
    idxj = Rcluster[j];
    idxij= i + idxj*p.surfaceVol;
    if(LR_couplings[idxij] > 1e-7 && clusterDef[idxj] != clusterNum ){
      probsw = 1 - exp(-2*phi_arr[i]*phi_arr[idxj]*LR_couplings[idxij]);
      randsw = unif(rng);
      if(randsw < probsw) {
	sqnl_wc_size++;
	swendsenWangClusterAddLR(idxj, s, cSpin, clusterNum, clusterDef, 
				   Rcluster, LR_couplings, phi_arr, p);
	
      }    
    }
  }
}

void clusterPossibleLR(int i, int *s, int cSpin,
		       bool *Pcluster, vector<int> &Rcluster,
		       Param p) {

  //The site possibily belongs to the cluster...
  Pcluster[i] = true;
  // ...so record it
  Rcluster.push_back(i);
  sqnl_wc_poss++;

  //If the neighbor's spin matches the cluster spin (cSpin)
  //then it is a possible cluster candidate. This is all we
  //need to identify at this point. The routine checks that
  //the candidate site is both not already identified, and
  //that it has the correct spin.

  //Forward in T
  if(!Pcluster[tp(i,p)] && s[tp(i,p)] == cSpin) {
    //cout<<"->tp";
    clusterPossibleLR(tp(i,p), s, cSpin, Pcluster, Rcluster, p);
  }

  //Forward in X
  if(!Pcluster[xp(i,p)] && s[xp(i,p)] == cSpin) {
    //cout<<"->xp";
    clusterPossibleLR(xp(i,p), s, cSpin, Pcluster, Rcluster, p);
  }

  //Backward in T
  if(!Pcluster[ttm(i,p)] && s[ttm(i,p)] == cSpin) {
    //cout<<"->tm";
    clusterPossibleLR(ttm(i,p), s, cSpin, Pcluster, Rcluster, p);
  }

  //Backward in X
  if(!Pcluster[xm(i,p)] && s[xm(i,p)] == cSpin) {
    //cout<<"->xm";
    clusterPossibleLR(xm(i,p), s, cSpin, Pcluster, Rcluster, p);
  }
}

void wolffUpdateLR(double *phi, int *s, Param p,
		   double *LR_couplings,
		   double &delta_mag_phi, int iter) {
  
  init_connectivity(p);

  sqnl_wc_calls++;
  sqnl_wc_poss = 0;
  
  //Choose a random spin and identify the possible cluster
  //with a flood fill.
  int i = int(unif(rng) * p.surfaceVol);
  int cSpin = s[i];
  
  // The site belongs to the cluster, so flip it.
  sqnl_wc_size = 1;
  s[i] *= -1;
  phi[i] *= -1;

  //This function is recursive and will call itself
  //until all attempts `to increase the cluster size
  //have failed.
  wolffClusterAddLR(i, s, cSpin, LR_couplings, phi, p);

  sqnl_wc_ave += sqnl_wc_size;

  if( iter%p.n_skip == 0) {
    setprecision(4);
    cout<<"Using "<<p.n_cluster<<" Wolff hits."<<endl; 
    cout<<"Ave. cluster size at iter "<<iter<<" = "<<sqnl_wc_ave<<"/"<<sqnl_wc_calls<<" = "<<1.0*sqnl_wc_ave/sqnl_wc_calls<<endl;
  }
  for(int a=0; a<p.surfaceVol; a++) free(added[a]);
  free(added);
}

void wolffClusterAddLR(int i, int *s, int cSpin, 
		       double *LR_couplings, double *phi, Param p) {
   
  double phi_lc = phi[i];
  
#ifdef USE_OMP

  int newSites = 0;
  int chunk = CACHE_LINE_SIZE/sizeof(double);
  int lc_sze = p.surfaceVol/omp_get_num_threads();

  //We now loop over the possible lattice sites, adding sites
  //(creating bonds) with the specified LR probablity.
#pragma omp parallel 
  {
    int newSites_local = 0;
    int lc_pos = 0;
    double prob = 0.0;
    double rand = 0.0;
    //double added_local[lc_sze];
    //for(int k=0; k<lc_sze; k++) {
    //added_local[k] = added[i][lc_sze*omp_get_thread_num()+k];
    //}  
#pragma omp for nowait 
    for(int j=0; j<omp_get_num_threads(); j++) {
      for(int k=0; k<lc_sze; k++) {
	int idx = j*lc_sze+k;
	if(s[idx] == cSpin) {
	  prob = 1 - exp(2*phi_lc*phi[idx]*LR_couplings[idx + i*p.surfaceVol]);
	  rand = unif(rng);
	  if(rand < prob) {
	    sqnl_wc_size++;
	    added[i][idx] = 1;
	    ++newSites_local;
	  }
	}
      }
    }
#pragma omp atomic
    newSites += newSites_local;
  }
  
  
  //cout<<newSites<<endl;
  
  if(newSites > 0) {
    //'added' now contains all of the spins on the wave front. Each
    //of these new sites must be explored, but no new exploration 
    //need test these sites. First, sort the added array so that all
    //the hits are at the beginning, order is unimportant.
    int pos = 0;
    int l = 0;
    for(l=0; l<p.surfaceVol; l++) {
      if(added[i][l] > 0) {
	added[i][pos] = l;
	pos++;
      }
    }

    //for(l=0; l<pos; l++) cout<<added[i][l]<<" ";
    //cout<<endl;

    //sort(added.begin(), added.begin()+p.surfaceVol, greater<int>());

    // These sites belongs to the cluster, so flip them.
    for(int k=0; k<newSites; k++) {
      s[added[i][k]] *= -1;
      phi[added[i][k]] *= -1;
    }
    
    for(int k=0; k<newSites; k++)
      wolffClusterAddLR(added[i][k], s, cSpin, LR_couplings, phi, p);
  }
  
#elif

  double prob = 0.0;
  double rand = 0.0;
  
  //We now loop over the possible lattice sites, adding sites
  //(creating bonds) with the specified LR probablity.
  for(int j=0; j<p.surfaceVol; j++) {
    if(s[j] == cSpin) {
      idx = j;
      prob = 1 - exp(2*phi[i]*phi[idx]*LR_couplings[i + idx*p.surfaceVol]);
      rand = unif(rng);
      if(rand < prob) {
	sqnl_wc_size++;
	// The site belongs to the cluster, so flip it.
	s[idx] *= -1;
	phi[idx] *= -1;  
	wolffClusterAddLR(idx, s, cSpin, Rcluster, LR_couplings, phi, p);
      }
    }
  }
#endif
}

