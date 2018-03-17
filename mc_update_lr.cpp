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
vector<int> toCheck;
int sitesToCheck;
int clusterIter;

void init_connectivity(Param p) {
  
  added = (int**)malloc(p.surfaceVol*sizeof(int*));  
  for(int i=0; i<p.surfaceVol; i++) {
    added[i] = (int*)malloc((p.surfaceVol+1)*sizeof(int));
#ifdef USE_OMP
#pragma omp parallel for
#endif
    for(int j=0; j<p.surfaceVol+1; j++) {
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
		       double *LR_couplings, double *denom,
		       double &delta_mag_phi, int iter) {

  delta_mag_phi = 0.0;
  
  int s_old     = 0;
  int delta_mag = 0;
  int Lt =p.Lt;
  int S1= p.S1;
  int x_len = S1/2 + 1;
  int t_len = Lt/2 + 1;

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
      cout<<"ERROR s and phi NOT aligned! iter = "<<iter<<" (MetroUpdate LR)"<<endl;
      cout<<"S["<<i<<"] = "<<s[i]<<" and phi["<<i<<"] = "<<phi<<endl;
      for (int j = 0; j < p.surfaceVol; j++) {    
	cout<<"S["<<j<<"] = "<<s[j]<<" and phi["<<j<<"] = "<<phi_arr[j]<<endl;
      }
      exit(0);
    }
    
    //PE
    DeltaE += lambda_p*(phi_new_sq*phi_new_sq - phi_sq*phi_sq);
    DeltaE += musqr_p *(phi_new_sq            - phi_sq);
    
    //KE
    double pmpn = phi-phi_new;
    double pnsmps = 0.5*(phi_new_sq-phi_sq);

#ifdef USE_OMP

    int t1,x1;
    t1 = i / S1;
    x1 = i % S1;

    int chunk = CACHE_LINE_SIZE/sizeof(double);
#pragma omp parallel
    {
      double local = 0.0;
      double val = 0.0;
      int t2,dt,x2,dx;
#pragma omp for nowait
      for(int j=0; j<p.surfaceVol; j+=chunk) {
	for(int k=0; k<chunk; k++) {

	  if( (j+k) != i ) {
	    //Index divided by circumference, using the int floor 
	    //feature/bug, gives the timeslice index.
	    t2 = (j+k) / S1;
	    dt = abs(t2-t1) > Lt/2 ? Lt - abs(t2-t1) : abs(t2-t1);	    
	    
	    //The index modulo the circumference gives the spatial index.
	    x2 = (j+k) % S1;            
	    dx = abs(x2-x1) > S1/2 ? S1-abs(x2-x1) : abs(x2-x1);      

	    val = ((pmpn*phi_arr[j+k] + pnsmps)*
		   LR_couplings[dx + dt*x_len]*LR_couplings[dx + dt*x_len]);
	    
	    local += val;
	  }
	}
      }
#pragma omp atomic
      DeltaE += local;
    }
#else
    int t1,t2,dt,x1,x2,dx;

    t1 = i / S1;
    x1 = i % S1;
    
    for(int j=0; j<p.surfaceVol; j++) {
      
      if(i!=j) {
	//Index divided by circumference, using the int floor 
	//feature/bug, gives the timeslice index.
	t2 = (j) / S1;
	dt = abs(t2-t1) > Lt/2 ? Lt - abs(t2-t1) : abs(t2-t1);
	
	//The index modulo the circumference gives the spatial index.
	x2 = (j) % S1;            
	dx = abs(x2-x1) > S1/2 ? S1-abs(x2-x1) : abs(x2-x1);      
	
	DeltaE += (pmpn*phi_arr[j] + pnsmps)*
	  LR_couplings[dx + dt*x_len]*LR_couplings[dx + dt*x_len];
	
      }
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
    cout<<"At iter "<<iter<<" the Metro acceptance rate is "<<(double)sqnl_accept/(double)sqnl_tries<<endl<<endl;
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
    //cout<<"Using "<<p.n_cluster<<" SW hits."<<endl; 
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

#ifdef USE_OMP
  init_connectivity(p);
#endif

  sqnl_wc_calls++;
  sqnl_wc_poss = 0;
  
  //Choose a random spin.
  int i = int(unif(rng) * p.surfaceVol);
  int cSpin = s[i];
  
  // The site belongs to the cluster, so flip it.
  sqnl_wc_size = 1;
  s[i] *= -1;
  phi[i] *= -1;

  //This function is recursive and will call itself
  //until all attempts to increase the cluster size
  //have failed.
  wolffClusterAddLR(i, s, cSpin, LR_couplings, phi, p);
  
  sqnl_wc_ave += sqnl_wc_size;

  if( iter%p.n_skip == 0) {
    setprecision(4);
    //cout<<"Using "<<p.n_cluster<<" Wolff hits."<<endl; 
    cout<<"Ave. cluster size at iter "<<iter<<" = "<<sqnl_wc_ave<<"/"<<sqnl_wc_calls<<" = "<<1.0*sqnl_wc_ave/sqnl_wc_calls<<endl;
  }
#ifdef USE_OMP
  for(int a=0; a<p.surfaceVol; a++) free(added[a]);
  free(added);
  toCheck.clear();
#endif
}

void wolffClusterAddLR(int i, int *s, int cSpin, 
		       double *LR_couplings, double *phi, Param p) {
  
  double phi_lc = phi[i];
  int S1 = p.S1;
  int Lt = p.Lt;
  int x_len = S1/2 + 1;
  int t_len = Lt/2 + 1;
  int arr_len = x_len*t_len;

#ifdef USE_OMP
  
  //This implementation parallelises over the entire surface of the lattice.
  //It performs a boolean check to see if the candidate site has the
  //correct spin, then performs the probablisitic test to add the site.
  
  int newSites = 0;
  int omp_nt = omp_get_num_threads();
  int lc_sze = p.surfaceVol/omp_nt;
  int t1,x1;
  t1 = i / p.S1;
  x1 = i % p.S1;

  //We now loop over the possible lattice sites, adding sites
  //(creating bonds) with the specified LR probablity.
#pragma omp parallel 
  {
    int newSites_local = 0;
    double prob = 0.0;
    double rand = 0.0;
    double added_local[lc_sze];
    for(int k=0; k<lc_sze; k++) {
      added_local[k] = -1;
    }
    int t2,dt,x2,dx;

#pragma omp for nowait 
    for(int j=0; j<p.surfaceVol; j+=lc_sze) {
      for(int k=0; k<lc_sze; k++) {
	int idx = j+k;
	if(s[idx] == cSpin) {
	  
	  //Index divided by circumference, using the int floor feature/bug,
	  //gives the timeslice index.
	  t2 = (idx) / p.S1;
	  dt = abs(t2-t1) > p.Lt/2 ? p.Lt - abs(t2-t1) : abs(t2-t1);
	  
	  //The index modulo the circumference gives the spatial index.
	  x2 = (idx) % p.S1;            
	  dx = abs(x2-x1) > p.S1/2 ? p.S1 - abs(x2-x1) : abs(x2-x1);      
	  
	  prob = 1 - exp(2*phi_lc*phi[idx]*LR_couplings[dx+dt*x_len]);
	  rand = unif(rng);
	  if(rand < prob) {
	    sqnl_wc_size++;
	    added_local[k] = 1;
	    ++newSites_local;
	  }
	}
      }
    }
#pragma omp critical
    newSites += newSites_local;
#pragma omp for
    for(int j=0; j<p.surfaceVol; j+=lc_sze) {
      for(int k=0; k<lc_sze; k++) {
	int idx = j+k;
	if(added_local[k] > 0) added[i][idx] = added_local[k];
      }
    }
  }
    
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

    // These sites belongs to the cluster, so flip them.
    for(int k=0; k<newSites; k++) {
      s[added[i][k]] *= -1;
      phi[added[i][k]] *= -1;
    }    
    for(int k=0; k<newSites; k++)
      wolffClusterAddLR(added[i][k], s, cSpin, LR_couplings, phi, p);
  }

#else

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
      
      prob = 1 - exp(2*phi_lc*phi[j]*LR_couplings[dx + dt*x_len]);
      rand = unif(rng);
      if(rand < prob) {
	sqnl_wc_size++;
	// The site belongs to the cluster, so flip it.
	s[j] *= -1;
	phi[j] *= -1;  
	wolffClusterAddLR(j, s, cSpin, LR_couplings, phi, p);
      }
    }
  }
  
#endif
  
}

#if 0

  //This implementation attempts to parallelises only over those sites which 
  //are candidate sites. This has the advantage that the number of computations 
  //reduce in number as the cluster is filled.

  int newSites = 0;
  int lc_sze = toCheck.size()/omp_get_num_threads();
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
    for(int j=0; j<toCheck.size(); j+=lc_sze) {
      for(int k=0; k<lc_sze; k++) { 
	int idx = toCheck[j+k];

	//cout<<"Checking: j="<<j<<" k="<<k<<" size="<<toCheck.size()<<" toCheck="<<toCheck[j+k]<<endl;
	
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
	  cout<<"Hit: j="<<j<<" k="<<k<<" size="<<toCheck.size()<<" toCheck="<<toCheck[j+k]<<endl;
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
	added[i][added_local[k][0]] = 1;
	s[added_local[k][0]] *= -1;
	phi[added_local[k][0]] *= -1;
	cout<<"Added site "<<added_local[k][0]<<endl;
       	
      }
    }
  }
  
  sitesToCheck += newSites;
  
  if(newSites > 0) {
    cout<<"SITES TO ADD!"<<endl;
    //'added' now contains all of the spins on the wave front. Each
    //of these new sites must be explored, but no new exploration 
    //need test these sites. First, sort the added array so that all
    //the hits are at the beginning, order is unimportant.
    int pos = 0;
    int l = 0;
    for(l=0; l<p.surfaceVol; l++) {
      if(added[i][l] > 0) {
	cout<<"ADDING "<<l<<endl;
	added[i][pos] = l;
	pos++;

	cout<<"REDUCE THE SEARCH SPACE"<<endl;
	cout<<"Removing element "<<added_local[k][1]<<" (idx="<<added_local[k][0]<<")"<<endl;
	toCheck.erase (toCheck.begin() + added_local[k][1]);
	cout<<"SEARCH SPACE REDUCED"<<endl;

      }
    }
    cout<<"ADDED!"<<endl;
  }

  added[i][p.surfaceVol] = newSites;
  
  cout<<"clusterIter: "<<clusterIter<<" SITES TO CHECK: "<<sitesToCheck;
  cout<<" Checked site: "<<i<<" New Sites: "<<newSites;
  
  for(int k=0; k<added[i][p.surfaceVol]; k++) {
    cout<<" Going to "<<added[i][k]<<" k="<<k<<" loop max="<<added[i][p.surfaceVol]<<endl;
    sitesToCheck -= 1;
    wolffClusterAddLR(added[i][k], s, cSpin, LR_couplings, phi, p);
  }

#ifdef USE_OMP
  //update the toCheck array with all the sites of the same sign.
  toCheck.push_back(i);
  for(int j=0; j<p.surfaceVol; j++) {
    if(phi[j]*cSpin > 0 && j != i) {
      toCheck.push_back(j);
    }
  }
  //cout<<endl<<"START: "<<toCheck.size()<<endl;
  for(int j=0; j<p.surfaceVol; j++) {
    //if(phi[j]*cSpin > 0) cout<<j<<":"<<phi[j]<<" ";
  }
  
#endif




#endif
