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

int metropolisUpdateSR(double *phi_arr, int *s, Param &p,
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
      delta_mag_phi += phi_new - phi_arr[i];
      phi_arr[i] = phi_new;
      sql_accept += 1;
      s[i] = (phi_new > 0) ? 1 : -1;
      delta_mag += s[i] - s_old;
    }
    else if ( unif(rng)  < exp(-DeltaE)) {
      //  cout<< " Acepted  " << endl;
      s_old = s[i];
      delta_mag_phi += phi_new - phi_arr[i];
      phi_arr[i] = phi_new;
      sql_accept += 1;
      s[i] = (phi_new > 0) ? 1 : -1;
      delta_mag += s[i] - s_old;
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
    cout<<"and delta_phi is "<<p.delta_phi<<endl;
  }

  return delta_mag;
}

void swendsenWangUpdateSR(double *phi_arr, int *s, Param p,
			   double &delta_mag_phi, int iter) {
  
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

void wolffUpdateSR(double *phi_arr, int *s, Param p,
		    double &delta_mag_phi, int iter) {
  
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
    cout<<"S/T cluster growth ratio at iter "<<iter<<" = "<<1.0*sql_wc_s_size/sql_wc_ave<<":"<<1.0*sql_wc_t_size/sql_wc_ave<<endl;
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



/*
//------------- Monte Carlo Update  ----------------//
void runMonteCarloSR(Param p) {
  
  //Create structure to hold observable data.
  observables obs(p.n_meas);
  
  //Create the LR couplings 
  //(if needed)
  double *LR_couplings;
  createLRcouplings(LR_couplings, p);

  //Create structures to hold spin and phi values.
  int *s = (int*)malloc(p.surfaceVol*sizeof(int));
  double *phi = (double*)malloc(p.surfaceVol*sizeof(double));
  for(int i = 0;i < p.surfaceVol; i++) {
    phi[i] = 2.0*unif(rng) - 1.0;
    s[i] = (phi[i] > 0) ? 1:-1;
    obs.MagPhi += phi[i];
  }  
  thermaliseSR(phi, s, p, obs.delta_mag_phi);

  //Running correlation function arrays.
  double *run_corr_t = (double*)malloc((p.Lt/2+1)*sizeof(double));
  for(int i=0; i<p.Lt/2+1; i++) {
    run_corr_t[i] = 0.0;
  }
  double *run_corr_s = (double*)malloc((p.S1/2+1)*sizeof(double));
  for(int i=0; i<p.S1/2+1; i++) {
    run_corr_s[i] = 0.0;
  }

  //Individual correlation functions.
  double **ind_corr_t = (double**)malloc((p.n_meas)*sizeof(double*));
  for(int i=0; i<p.n_meas; i++) {
    ind_corr_t[i] = (double*)malloc((p.Lt/2+1)*sizeof(double));
    for(int j=0; j<p.Lt/2+1; j++) {
      ind_corr_t[i][j] = 0.0;
    }
  }
  double **ind_corr_s = (double**)malloc((p.n_meas)*sizeof(double*));
  for(int i=0; i<p.n_meas; i++) {
    ind_corr_s[i] = (double*)malloc((p.S1/2+1)*sizeof(double));
    for(int j=0; j<p.S1/2+1; j++) {
      ind_corr_s[i][j] = 0.0;
    }
  }

  //The correlation functions must be weighted by how frequently
  //the i,j separation occurs  
  //correlation norms array.
  int **corr_norms = (int**)malloc((p.S1/2+1)*sizeof(int*));
  for(int i=0; i<p.S1/2+1; i++) {
    corr_norms[i] = (int*)malloc((p.Lt/2+1)*sizeof(int));
  }
  for(int i=0; i<p.S1/2+1; i++)
    for(int j=0; j<p.Lt/2+1; j++) 
      corr_norms[i][j] = 0;
  
  int s_idx = 0;
  int t_idx = 0;
  int s_size = p.S1;
  int t_size = p.Lt;

  //loop over sink/source *theta*
  for(int is=0; is<s_size; is++)
    for(int js=0; js<s_size; js++) {      
      s_idx = abs(is-js);
      if(s_idx > s_size/2) s_idx = s_size - s_idx ;      
      //loop over sink/source *temporal*
      for(int il=0; il<t_size; il++) 	  
	for(int jl=0; jl<t_size; jl++) {
	  t_idx = abs(il-jl);
	  if(t_idx >= t_size/2) t_idx = t_size - t_idx;
	  ++corr_norms[s_idx][t_idx];
	}
    }
    
  int idx = 0;
  double norm;

  long long metro = 0.0;
  long long cluster = 0.0;
  for(int iter = p.n_therm; iter < p.n_therm + p.n_skip*p.n_meas; iter++) {

    auto start = std::chrono::high_resolution_clock::now();    
    for(int i=0; i<p.n_cluster; i++) {
      if(p.useWolff) 
	wolffUpdateSR(phi, s, p, obs.delta_mag_phi, iter);      
      else 
	swendsenWangUpdateSR(phi, s, p, obs.delta_mag_phi, iter);
    }
    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    cluster += std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    
    start = std::chrono::high_resolution_clock::now();
    metropolisUpdateSR(phi, s, p, obs.delta_mag_phi, iter);
    elapsed = std::chrono::high_resolution_clock::now() - start;
    metro += std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    
    //Take measurements.
    if((iter+1) % p.n_skip == 0) {

      //update running average, dump to stdout. 
      //'idx' is ++ incremented in this function.
      measure(obs, phi, s, LR_couplings, idx, p);
      
      //Visualisation tools
      //visualiserSqr(phi, avePhiAb*norm, p);

      norm = 1.0/(idx);
      
      correlators(ind_corr_t, idx-1, run_corr_t, 
		  true,  phi, obs.avePhi*norm, p);
      correlators(ind_corr_s, idx-1, run_corr_s, 
		  false, phi, obs.avePhi*norm, p);
      
      
      //correlatorsSRImp(ind_corr_t, idx-1, run_corr_t, 
      //true,  phi, obs.avePhi*norm, s, p);
      //correlatorsSRImp(ind_corr_s, idx-1, run_corr_s, 
      //false, phi, obs.avePhi*norm, s, p);
      
      //Jacknife and dump the data
      if(idx%10 == 0) {	
	writeObservables(ind_corr_t, run_corr_s, ind_corr_t, run_corr_s,
			 corr_norms, idx, obs, p);
      }      
    }
  }

  free(run_corr_t);
  free(run_corr_s);
  free(corr_norms);
  free(s);
  free(phi);

}



void thermalise(double *phi, int *s, 
		Param p, double &delta_mag_phi) {

  long long metro = 0.0;
  long long cluster = 0.0;

  cout<<"Begin Metropolis cool down."<<endl;
  //Perform n_therm/2 pure metropolis hits during the hot phase.
  auto start = std::chrono::high_resolution_clock::now();        
  for(int iter = 0; iter < p.n_therm/4; iter++) {
    metropolisUpdateSR(phi, s, p, delta_mag_phi, iter);  
    if((iter+1)%100 == 0) cout<<"Metro cool down iter "<<iter+1<<endl;
  }
  auto elapsed = std::chrono::high_resolution_clock::now() - start;
  metro = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
  cout<<"Metro cool down complete. "<<p.n_therm/4<<" iters in "<<metro/(1.e6)<<"s"<<endl;
  metro = 0.0;

  //Cluster/Metro steps
  for(int iter = 0; iter < p.n_therm; iter++) {
    start = std::chrono::high_resolution_clock::now();        
    //Cluster steps.
    for(int i=0; i<p.n_cluster; i++) {
      if(p.useWolff) 
	wolffUpdateSR(phi, s, p, delta_mag_phi, iter);      
      else 
	swendsenWangUpdateSR(phi, s, p, delta_mag_phi, iter);
    }
    elapsed = std::chrono::high_resolution_clock::now() - start;
    cluster += std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

    //Metro step.
    start = std::chrono::high_resolution_clock::now();
    metropolisUpdateSR(phi, s, p, delta_mag_phi, iter);  
    elapsed = std::chrono::high_resolution_clock::now() - start;
    metro += std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

    if((iter+1)%p.n_skip == 0) {
      cout<<"Therm sweep "<<iter+1<<endl;
      cout<<"Average cluster time = "<<cluster/((iter+1)*1.0e6)<<"s"<<endl;
      cout<<"Average Metrop. time = "<<metro/((iter+1)*1.0e6)<<"s"<<endl;
    }
  }
  cout<<endl<<"Thermalisaton complete."<<endl<<endl;  

}

*/
