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

#include "util.h"
#include "hyp_util.h"
#include "data_proc.h"
#include "data_io.h"
#include "monte_carlo_square_nonlocal.h"

using namespace std;

extern int seed;
extern mt19937 rng;
extern uniform_real_distribution<double> unif;

//ferromagnetic coupling
const double J = 1.0;

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

//Declare variables to implement Wolff algorithm
int sqnl_wc_ave = 0;
int sqnl_wc_size = 0;
int sqnl_wc_calls = 0;
int sqnl_wc_poss = 0;
int sqnl_wc_t_size = 0;
int sqnl_wc_s_size = 0;

double actionSqNL(double *phi_arr, int *s, Param p,
		  double *LR_couplings, 
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
  
  for (int i = 0; i < p.surfaceVol; i++) {    
    phi = phi_arr[i];
    phi_sq = phi*phi;

    //PE terms
    PE += lambda_p * phi_sq*phi_sq;
    PE += musqr_p  * phi_sq;

    //KE terms
    double val = 0.0;
#pragma omp parallel for private(val) reduction(+:KE)
    for(int j=0; j<p.surfaceVol; j++) {      
      //Here we are overcounting by a factor of two, hence the 0.25
      //coefficient.
      //FIXME
      val = 0.25*(phi - phi_arr[j])*(phi - phi_arr[j])*LR_couplings[i+j*p.surfaceVol]*LR_couplings[i+j*p.surfaceVol];
      KE += val;
    }
  }
  return PE + KE;
}

int sqnl_accept = 0;
int sqnl_tries  = 0;

int metropolisUpdateSqNL(double *phi_arr, int *s, Param &p,
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
      printf("ERROR s and phi NOT aligned! (MUP)\n");
      exit(0);
    }
    
    //PE
    DeltaE += lambda_p*(phi_new_sq*phi_new_sq - phi_sq*phi_sq);
    DeltaE += musqr_p *(phi_new_sq            - phi_sq);
    
    //KE
    double pmpn = phi-phi_new;
    double pnsmps = 0.5*phi_new_sq-phi_sq;
    double val = 0.0;
#pragma omp parallel for private(val) reduction(+:DeltaE)
    for(int j=0; j<p.surfaceVol; j++) {
      val = (pmpn*phi_arr[j] + pnsmps)*LR_couplings[i+j*p.surfaceVol]*LR_couplings[i+j*p.surfaceVol];
      DeltaE += val;
    }
    
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
    else if ( unif(rng)  < exp(-DeltaE)) {
      //  cout<< " Acepted  " << endl;
      s_old = s[i];
      delta_mag_phi += phi_new - phi_arr[i];
      phi_arr[i] = phi_new;
      sqnl_accept += 1;
      s[i] = (phi_new > 0) ? 1 : -1;
      delta_mag += s[i] - s_old;
    }     
  }// end loop over lattice volume 

  /*  
  // TUNING ACCEPTANCE 
  if (iter < p.n_therm/2 && (iter+1) % p.n_skip/10 == 0) {
    if ((double) sqnl_accept / (double) sqnl_tries < 0.5) {
      p.delta_phi -= 0.001;
    } else {
      p.delta_phi += 0.001;
    }
    if(p.n_cluster*1.0*sqnl_wc_ave/sqnl_wc_calls < p.surfaceVol && iter > p.n_skip) {
      p.n_cluster++;
    } else {
      p.n_cluster--;
      if(p.n_cluster < 3) p.n_cluster++;
    }
  }
  
  if( iter < p.n_therm ) {
    cout<<"At iter "<<iter<<" the Acceptance rate is "<<(double)sqnl_accept/(double)sqnl_tries<<endl;
    cout<<"and delta_phi is "<<p.delta_phi<<endl;
  }
  */
  return delta_mag;
}

void swendsenWangUpdateSqNL(double *phi_arr, int *s, Param p,
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
      clusterPossibleSqNL(i, s, clusterSpin[clusterNum], 
			  Pcluster, Rcluster, p);
      
      //This function will call itself recursively until it fails to 
      //add to the cluster
      SWclusterAddSqNL(i, s, clusterSpin[clusterNum], clusterNum, clusterDef, 
		       Rcluster, LR_couplings, phi_arr, p);

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

double probsw = 0.0;
double randsw = 0.0;
int idxj = 0;
int idxij = 0;

void SWclusterAddSqNL(int i, int *s, int cSpin, int clusterNum, 
		      int *clusterDef, vector<int> Rcluster, 
		      double *LR_couplings,
		      double *phi_arr, Param p) {

  //The site belongs to the (clusterNum)th cluster
  clusterDef[i] = clusterNum;

  //We now loop over the possible lattice sites, adding sites
  //(creating bonds) with the specified LR probablity.
  for(int j=0; j<sqnl_wc_poss; j++) {
    idxj = Rcluster[j];
    idxij= i + idxj*p.surfaceVol;
    if(LR_couplings[idxij] > 1e-7 && clusterDef[idxj] != clusterNum ){
      probsw = 1 - exp(-2*J*phi_arr[i]*phi_arr[idxj]*LR_couplings[idxij]);
      randsw = unif(rng);
      if(randsw < probsw) {
	sqnl_wc_size++;
	SWclusterAddSqNL(idxj, s, cSpin, clusterNum, clusterDef, 
			 Rcluster, LR_couplings, phi_arr, p);
	
      }    
    }
  }
}

void wolffUpdateSqNL(double *phi_arr, int *s, Param p,
		     double *LR_couplings,
		     double &delta_mag_phi, int iter) {
  
  sqnl_wc_calls++;
  sqnl_wc_poss = 0;
  
  //Tracks which sites are potentially in the cluster.
  bool *Pcluster = new bool[p.surfaceVol];
#pragma omp parallel for
  for (int i = 0; i < p.surfaceVol; i++)
    Pcluster[i] = false;

  //Records which sites are potentially in the cluster.
  //This will have a modifiable size, hence the vector
  //style.
  vector<int> Rcluster;
  
  //Choose a random spin and identify the possible cluster
  //with a flood fill.
  int i = int(unif(rng) * p.surfaceVol);
  int cSpin = s[i];
  clusterPossibleSqNL(i, s, cSpin, Pcluster, Rcluster, p);

  //This function is recursive and will call itself
  //until all attempts `to increase the cluster size
  //have failed.
  sqnl_wc_size = 1;
  clusterAddSqNL(i, s, cSpin, Rcluster, LR_couplings, phi_arr, p);

  sqnl_wc_ave += sqnl_wc_size;

  if( iter%p.n_skip == 0 && iter < p.n_therm) {
    setprecision(4);
    cout<<"Using "<<p.n_cluster<<" Wolff hits."<<endl; 
    cout<<"Ave. cluster size at iter "<<iter<<" = "<<sqnl_wc_ave<<"/"<<sqnl_wc_calls<<" = "<<1.0*sqnl_wc_ave/sqnl_wc_calls<<endl;

  }
  
  delete[] Pcluster;
}

void clusterPossibleSqNL(int i, int *s, int cSpin,
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
    clusterPossibleSqNL(tp(i,p), s, cSpin, Pcluster, Rcluster, p);
  }
  
  //Forward in X
  if(!Pcluster[xp(i,p)] && s[xp(i,p)] == cSpin) {
    //cout<<"->xp";
    clusterPossibleSqNL(xp(i,p), s, cSpin, Pcluster, Rcluster, p);
  }
  
  //Backward in T 
  if(!Pcluster[ttm(i,p)] && s[ttm(i,p)] == cSpin) {  
    //cout<<"->tm";
    clusterPossibleSqNL(ttm(i,p), s, cSpin, Pcluster, Rcluster, p);
  }
  
  //Backward in X
  if(!Pcluster[xm(i,p)] && s[xm(i,p)] == cSpin) {
    //cout<<"->xm";
    clusterPossibleSqNL(xm(i,p), s, cSpin, Pcluster, Rcluster, p);
  }
}

void clusterAddSqNL(int i, int *s, int cSpin, vector<int> &Rcluster, 
		    double *LR_couplings,
		    double *phi_arr, Param p) {
  
  // The site belongs to the cluster, so flip it.
  s[i] *= -1;
  phi_arr[i] *= -1;

  double prob = 0.0;
  double rand = 0.0;

  int idx = 0;
  
  //We now loop over the possible lattice sites, adding sites
  //(creating bonds) with the specified LR probablity.
  for(int j=0; j<sqnl_wc_poss; j++) {
    idx = Rcluster[j];
    if(LR_couplings[i + idx*p.surfaceVol] > 1e-7){
      prob = 1 - exp(2*phi_arr[i]*phi_arr[idx]*LR_couplings[i + idx*p.surfaceVol]);
      rand = unif(rng);
      if(rand < prob) {
	sqnl_wc_size++;
	clusterAddSqNL(idx, s, cSpin, Rcluster, LR_couplings, phi_arr, p);
      }    
    }
  }
}

//------------- Monte Carlo Update  ----------------//
void runMonteCarloSqNL(vector<Vertex> &NodeList, Param p) {
  
  double KE = 0.0, PE = 0.0;
  double mag_phi = 0.0;
  double delta_mag_phi = 0.0;
  
  //Create the 1/r^{(2+\sigma)} LR couplings;
  double sigma = p.sigma;
  double *LR_couplings = new double[p.surfaceVol*p.surfaceVol];
  double r_x = 0.0;
  double r_y = 0.0;

  for(int i=0; i<p.surfaceVol; i++) {
#pragma omp parallel for
    for(int j=0; j<p.surfaceVol; j++) {
      
      r_x = abs(i%p.S1 - j%p.S1);
      if(r_x > p.S1/2) r_x = p.S1 - r_x;
      
      r_y = abs(i/p.Lt - j/p.Lt);
      if(r_y > p.Lt/2) r_y = p.Lt - r_y;
      
      if(i != j) LR_couplings[i + j*p.surfaceVol] = pow(sqrt(r_x*r_x + r_y*r_y), -(2+sigma));
      //Quick symmetry check...
      //if( (i == 10 && j == 20) || (i==20 && j == 10)) cout<<i<<" "<<j<<" "<<LR_couplings[i + j*p.surfaceVol]<<endl;
    }
  }
  
  int *s = (int*)malloc(p.surfaceVol*sizeof(int));
  double *phi = (double*)malloc(p.surfaceVol*sizeof(double));
  for(int i = 0;i < p.surfaceVol; i++) {
    phi[i] = 2.0*unif(rng) - 1.0;
    s[i] = (phi[i] > 0) ? 1:-1;
    mag_phi += phi[i];
  }  
  thermaliseSqNL(phi, s, p, delta_mag_phi, LR_couplings);
  
  //Arrays holding measurements for error analysis
  double E_arr[p.n_meas];
  double E2_arr[p.n_meas];
  double PhiAb_arr[p.n_meas];
  double Phi_arr[p.n_meas];
  double Phi2_arr[p.n_meas];
  double Phi4_arr[p.n_meas];
  for(int i=0; i<p.n_meas; i++) {
    E_arr[i]     = 0.0;
    E2_arr[i]    = 0.0;
    PhiAb_arr[i] = 0.0;
    Phi_arr[i]   = 0.0;
    Phi2_arr[i]  = 0.0;
    Phi4_arr[i]  = 0.0;
  }

  //Running averages
  double tmpE     = 0.0;  
  double aveE     = 0.0;
  double aveKE    = 0.0;
  double avePE    = 0.0;
  double aveE2    = 0.0;
  double avePhiAb = 0.0;
  double avePhi   = 0.0;
  double avePhi2  = 0.0;
  double avePhi4  = 0.0;
  double MagPhi   = 0.0;

  //Density of lattice sites
  double rhoVol  = 1.0/(double)p.surfaceVol;
  double rhoVol2 = rhoVol*rhoVol;

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
  double norm = 0.0;
  double val = 0.0;  
  for(int iter = p.n_therm; iter < p.n_therm + p.n_skip*p.n_meas; iter++) {
    for(int i=0; i<p.n_cluster; i++) {
      if(p.useWolff) 
	wolffUpdateSqNL(phi, s, p, LR_couplings, delta_mag_phi, iter);      
      else 
	swendsenWangUpdateSqNL(phi, s, p, LR_couplings, delta_mag_phi, iter);
    }    
    metropolisUpdateSqNL(phi, s, p, LR_couplings, delta_mag_phi, iter);
    
    //Take measurements.
    //------------------
    if((iter+1) % p.n_skip == 0) {
      
      tmpE   = actionSqNL(phi, s, p, LR_couplings, KE, PE);
      aveKE += rhoVol*KE;
      avePE += rhoVol*PE;
      aveE  += rhoVol*tmpE;
      aveE2 += rhoVol2*tmpE*tmpE;

      MagPhi = 0.0;      
#pragma omp parallel for private(val) reduction(+:MagPhi)
      for(int i = 0;i < p.surfaceVol; i++) {
	val = phi[i];
	MagPhi += val;
      }
      MagPhi *= rhoVol;
      
      avePhiAb  += abs(MagPhi);
      avePhi    += MagPhi;
      avePhi2   += MagPhi*MagPhi;
      avePhi4   += MagPhi*MagPhi*MagPhi*MagPhi;
      
      E_arr[idx]     = rhoVol*tmpE;
      E2_arr[idx]    = rhoVol*tmpE*tmpE;
      PhiAb_arr[idx] = abs(MagPhi);
      Phi_arr[idx]   = MagPhi;
      Phi2_arr[idx]  = MagPhi*MagPhi;
      Phi4_arr[idx]  = MagPhi*MagPhi*MagPhi*MagPhi;
      
      idx++;
      
      cout<<setprecision(8);
      norm = 1.0/(idx);

      //Dump to stdout
      cout<<"Measurement "<<(iter+1 - p.n_therm)/p.n_skip;
      cout<<" Sweep "<<iter+1<<endl;

      cout<<"Ave Energy= "<<aveE*norm<<endl;
      cout<<"Ave KE    = "<<aveKE*norm<<endl;
      cout<<"Ave PE    = "<<avePE*norm<<endl;
      cout<<"Ave |phi| = "<<avePhiAb*norm<<endl;
      cout<<"Ave phi   = "<<avePhi*norm<<endl;
      cout<<"Ave phi^2 = "<<avePhi2*norm<<endl;
      cout<<"Ave phi^4 = "<<avePhi4*norm<<endl;
      cout<<"Suscep    = "<<(avePhi2*norm-pow(avePhiAb*norm,2))/rhoVol<<endl;
      cout<<"Spec Heat = "<<(aveE2*norm-pow(aveE*norm,2))/rhoVol<<endl;
      cout<<"Binder    = "<<1.0-avePhi4/(3.0*avePhi2*avePhi2*norm)<<" ("<<norm*avePhi2*avePhi2/avePhi4<<")"<<endl;
      
      //Visualisation tool
      //visualiserSqr(phi, avePhiAb*norm, p);

      //Calculate correlaton functions and update the average.
      correlators(ind_corr_t, idx-1, run_corr_t, true,  phi, avePhi*norm, p);
      correlators(ind_corr_s, idx-1, run_corr_s, false, phi, avePhi*norm, p);

      double *auto_corr = (double*)malloc(idx*sizeof(double));
      if(idx > 10) autocorrelation(PhiAb_arr, avePhiAb*norm, idx, auto_corr);
      free(auto_corr);

      //dump the data
      if(idx%10 == 0) {

	double *jk_err_t = (double*)malloc((p.Lt/2+1)*sizeof(double));
	jackknife(ind_corr_t, run_corr_t, jk_err_t, 5, idx, p.Lt/2+1, p); 
	ofstream filet("correlators_t.dat");
	for(int i=0; i<p.Lt/2+1; i++) {
	  filet<<i<<" "<<(run_corr_t[i]/corr_norms[0][i])*norm;
	  filet<<" "<<(jk_err_t[i]/corr_norms[0][i])<<endl;
	  cout<<(run_corr_t[i]/corr_norms[0][i])*norm<<" ";
	  cout<<(jk_err_t[i]/corr_norms[0][i])<<" ";
	}
	filet.close();
	free(jk_err_t);
	cout<<endl;
	
	double *jk_err_s = (double*)malloc((p.S1/2+1)*sizeof(double));
	jackknife(ind_corr_t, run_corr_t, jk_err_s, 2, idx, p.S1/2+1, p);
	ofstream files("correlators_s.dat");
	for(int i=0; i<p.S1/2+1; i++) {
	  files<<i<<" "<<(run_corr_s[i]/corr_norms[i][0])*norm;
	  files<<" "<<(jk_err_s[i]/corr_norms[i][0])<<endl;
	  cout<<(run_corr_s[i]/corr_norms[i][0])*norm<<" ";
	  cout<<(jk_err_s[i]/corr_norms[i][0])<<" ";
	}
	files.close();
	free(jk_err_s);
      }
    }
  }


  
  free(run_corr_t);
  free(run_corr_s);
  free(corr_norms);
  free(s);
  free(phi);
}

void thermaliseSqNL(double *phi, int *s, Param p,
		    double &delta_mag_phi, double *LR_couplings) {
  
  long long metro = 0.0;
  long long cluster = 0.0;
  for(int iter = 0; iter < p.n_therm; iter++) {

    auto start = std::chrono::high_resolution_clock::now();    
    for(int i=0; i<p.n_cluster; i++) {
      if(p.useWolff) 
	wolffUpdateSqNL(phi, s, p, LR_couplings, delta_mag_phi, iter);      
      else 
	swendsenWangUpdateSqNL(phi, s, p, LR_couplings, delta_mag_phi, iter);
    }
    auto elapsed = std::chrono::high_resolution_clock::now() - start;

    cluster += std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    
    start = std::chrono::high_resolution_clock::now();
    metropolisUpdateSqNL(phi, s, p, LR_couplings, delta_mag_phi, iter);
    elapsed = std::chrono::high_resolution_clock::now() - start;

    metro += std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

    if((iter+1)%p.n_skip == 0) {
      cout<<"Therm sweep "<<iter+1<<" Cluster time = "<<cluster/(1e6)<<"s Metropolis time = "<<metro/(1e6)<<endl;
      cluster = 0.0;
      metro = 0.0;
    }
  }
}
