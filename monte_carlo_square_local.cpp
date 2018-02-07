#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>

#include "util.h"
#include "data_proc.h"
#include "data_io.h"
#include "monte_carlo_square_local.h"

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

// declare variables to implement Wolff algorithm
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

double actionSqL(double *phi_arr, int *s, Param p,
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

int sql_accept = 0;
int sql_tries  = 0;

int metropolisUpdateSqL(double *phi_arr, int *s, Param &p,
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

  // TUNING ACCEPTANCE 
  if (iter < p.n_therm/2 && (iter+1) % p.n_skip/10 == 0) {
    if ((double) sql_accept / (double) sql_tries < 0.5) {
      p.delta_phi -= 0.001;
    } else {
      p.delta_phi += 0.001;
    }
    if(p.n_wolff*1.0*sql_wc_ave/sql_wc_calls < p.surfaceVol && iter > p.n_skip) {
      p.n_wolff++;
    } else {
      p.n_wolff--;
      if(p.n_wolff < 3) p.n_wolff++;
    }
  }
  
  return delta_mag;
}

void swendsenWangUpdateSqL(double *phi_arr, int *s, Param p,
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
      
      cout<<"Site "<<i<<endl;
      
      //This function will call itself recursively until it fails to 
      //add to the cluster
      SWclusterAddSqL(i, s, clusterSpin[clusterNum], clusterNum, clusterDef, phi_arr, p);
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

void SWclusterAddSqL(int i, int *s, int cSpin, int clusterNum, 
		     int *clusterDef, double *phi_arr, Param p) {
  
  //The site belongs to the (clusterNum)th cluster
  clusterDef[i] = clusterNum;

  //ferromagnetic coupling
  const double J = 1.0;
  
  //If the (aligned) neighbor spin does not already belong to the
  // cluster, then try to add it to the cluster.
  // - If the site has already been added, then we may skip the test.

  //Forward in T
  if(clusterDef[tp(i,p)] == 0 && s[tp(i,p)] == cSpin) {
    if(unif(rng) < 1 - exp(-2*J*phi_arr[i]*phi_arr[tp(i,p)])) {
      //cout<<"->tp";
      SWclusterAddSqL(tp(i,p), s, cSpin, clusterNum, clusterDef, phi_arr, p);
    }
  }

  //Forward in X
  if(clusterDef[xp(i,p)] == 0 && s[xp(i,p)] == cSpin) {
    if(unif(rng) < 1 - exp(-2*J*phi_arr[i]*phi_arr[xp(i,p)])) {
      //cout<<"->xp";
      SWclusterAddSqL(xp(i,p), s, cSpin, clusterNum, clusterDef, phi_arr, p);
    }
  }
  
  //Backard in T 
  if(clusterDef[ttm(i,p)] == 0 && s[ttm(i,p)] == cSpin) {  
    if(unif(rng) < 1 - exp(-2*J*phi_arr[i]*phi_arr[ttm(i,p)])) {
      //cout<<"->tm";
      SWclusterAddSqL(ttm(i,p), s, cSpin, clusterNum, clusterDef, phi_arr, p);
    }
  }
  
  //Backward in X
  if(clusterDef[xm(i,p)] == 0 && s[xm(i,p)] == cSpin) {
    if (unif(rng) < 1 - exp(-2*J*phi_arr[i]*phi_arr[xm(i,p)])) {
      //cout<<"->xm";
      SWclusterAddSqL(xm(i,p), s, cSpin, clusterNum, clusterDef, phi_arr, p);
    }
  } 
}


void wolffUpdateSqL(double *phi_arr, int *s, Param p,
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
  clusterAddSqL(i, s, cSpin, cluster, phi_arr, p);

  sql_wc_ave += sql_wc_size;

  if( iter%p.n_skip == 0 && iter < p.n_therm) {
    setprecision(4);
    cout<<"Using "<<p.n_wolff<<" Wolff hits."<<endl; 
    cout<<"Ave. cluster size at iter "<<iter<<" = "<<sql_wc_ave<<"/"<<sql_wc_calls<<" = "<<1.0*sql_wc_ave/sql_wc_calls<<endl;
    cout<<"S/T cluster growth ratio at iter "<<iter<<" = "<<1.0*sql_wc_s_size/sql_wc_ave<<":"<<1.0*sql_wc_t_size/sql_wc_ave<<endl;
  }
  delete[] cluster;
}

void clusterAddSqL(int i, int *s, int cSpin,
		   bool *cluster, double *phi_arr, Param p) {
  
  // The site belongs to the cluster, so flip it.
  cluster[i] = true;
  s[i] *= -1;
  phi_arr[i] *= -1;

  // ferromagnetic coupling
  const double J = 1.0;
  
  // - If the (aligned) neighbor spin does not already belong to the
  // cluster, then try to add it to the cluster.
  // - If the site has already been added, then we may skip the test.

  //Forward in T
  if(!cluster[ tp(i,p) ] && s[tp(i,p)] == cSpin) {
    if(unif(rng) < 1 - exp(2*J*phi_arr[i]*phi_arr[tp(i,p)])) {
      sql_wc_size++;
      sql_wc_t_size++;
      //cout<<"->tp";
      clusterAddSqL(tp(i,p), s, cSpin, cluster, phi_arr, p);
    }
  }

  //Forward in X
  if(!cluster[ xp(i,p) ] && s[xp(i,p)] == cSpin) {
    if(unif(rng) < 1 - exp(2*J*phi_arr[i]*phi_arr[xp(i,p)])) {
      sql_wc_size++;
      sql_wc_s_size++;
      //cout<<"->xp";
      clusterAddSqL(xp(i,p), s, cSpin, cluster, phi_arr, p);
    }
  }

  
  //Backard in T 
  if(!cluster[ ttm(i,p) ] && s[ttm(i,p)] == cSpin) {  
    if(unif(rng) < 1 - exp(2*J*phi_arr[i]*phi_arr[ttm(i,p)])) {
      sql_wc_size++;
      sql_wc_t_size++;
      //cout<<"->tm";
      clusterAddSqL(ttm(i,p), s, cSpin, cluster, phi_arr, p);
    }
  }

  //Backward in X
  if(!cluster[ xm(i,p) ] && s[xm(i,p)] == cSpin) {
    if (unif(rng) < 1 - exp(2*J*phi_arr[i]*phi_arr[xm(i,p)])) {
      sql_wc_size++;
      sql_wc_s_size++;
      //cout<<"->xm";
      clusterAddSqL(xm(i,p), s, cSpin, cluster, phi_arr, p);
    }
  } 
}

//------------- Monte Carlo Update  ----------------//
void runMonteCarloSqL(vector<Vertex> &NodeList, Param p) {
  
  double KE = 0.0, PE = 0.0;
  double mag_phi = 0.0;
  double delta_mag_phi = 0.0;

  int *s = (int*)malloc(p.surfaceVol*sizeof(int));
  double *phi;
  phi = (double*)malloc(p.surfaceVol*sizeof(double));
  for(int i = 0;i < p.surfaceVol; i++) {
    phi[i] = 2.0*unif(rng) - 1.0;
    s[i] = (phi[i] > 0) ? 1:-1;
    mag_phi += phi[i];
  }  
  thermaliseSqL(phi, s, p, delta_mag_phi);
  
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

  //correlation function arrays. One keeps the running average,
  //one is temporary to pass around single measurent values. 
  double **corr_tmp = (double**)malloc((p.S1/2)*sizeof(double*));
  double **corr_ave = (double**)malloc((p.S1/2)*sizeof(double*));
  for(int i=0; i<p.S1/2; i++) {
    corr_tmp[i] = (double*)malloc((p.Lt/2)*sizeof(double));
    corr_ave[i] = (double*)malloc((p.Lt/2)*sizeof(double));
  }
  for(int i=0; i<p.S1/2; i++)
    for(int j=0; j<p.Lt/2; j++) 
      corr_ave[i][j] = 0.0;
  
  double **phi_sq_arr = (double**)malloc(p.S1*sizeof(double*));
  for(int i=0; i<p.S1; i++) {
    phi_sq_arr[i] = (double*)malloc(p.Lt*sizeof(double));
    for(int j=0; j<p.Lt; j++) phi_sq_arr[i][j] = 0.0;
  }

  int idx = 0;
  double norm;
  
  for(int iter = p.n_therm; iter < p.n_therm + p.n_skip*p.n_meas; iter++) {
    
    for(int i=0; i<p.n_wolff; i++) {
      //wolffUpdateSqL(phi, s, p, delta_mag_phi, iter);
      swendsenWangUpdateSqL(phi, s, p, delta_mag_phi, iter);
    }
    metropolisUpdateSqL(phi, s, p, delta_mag_phi, iter);
    
    //Take measurements.
    if((iter+1) % p.n_skip == 0) {
      
      for(int i = 0;i < p.S1; i++) {
	for(int j=0; j<p.Lt; j++) 
	  phi_sq_arr[i][j] += pow(phi[i + p.S1*j],2);
      }
      
      tmpE   = actionSqL(phi, s, p, KE, PE);
      aveKE += rhoVol*KE;
      avePE += rhoVol*PE;
      aveE  += rhoVol*tmpE;
      aveE2 += rhoVol2*tmpE*tmpE;

      MagPhi = 0.0;
      for(int i = 0;i < p.surfaceVol; i++) MagPhi += phi[i];
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
      cout<<"Measurement "<<(iter+1 - p.n_therm)/p.n_skip<<" Sweep "<<iter+1<<endl;
      cout<<"Ave Energy= "<<aveE*norm<<endl;
      cout<<"Ave KE    = "<<aveKE*norm<<endl;
      cout<<"Ave PE    = "<<avePE*norm<<endl;
      cout<<"Ave |phi| = "<<avePhiAb*norm<<endl;
      cout<<"Ave phi   = "<<avePhi*norm<<endl;
      cout<<"Ave phi^2 = "<<avePhi2*norm<<endl;
      cout<<"Ave phi^4 = "<<avePhi4*norm<<endl;
      cout<<"Suscep    = "<<(avePhi2*norm-pow(avePhiAb*norm,2))/rhoVol<<endl;
      cout<<"Spec Heat = "<<(aveE2*norm-pow(aveE*norm,2))/rhoVol<<endl;
      cout<<"Binder    = "<<1.0-avePhi4/(3.0*avePhi2*avePhi2*norm)<<endl;
      
      //Visualisation tools
      visualiserSqr(phi, avePhiAb*norm, p);
      //visualiserPhi2(phi_sq_arr, p, idx);      

      //Calculate correlaton functions and update the average.
      correlators(corr_tmp, corr_ave, idx, phi, avePhi*norm, p);
      
      ofstream filet("correlators_t.dat");
      for(int i=0; i<p.Lt/2; i++) {
	filet << i << " " << corr_tmp[0][i] << endl;
      }
      filet.close();
      
      ofstream files("correlators_s.dat");
      for(int i=0; i<p.S1/2; i++) {
	files << i << " " << corr_tmp[i][0] << endl;
      }
      files.close();
      //if(idx%50 == 0) corr_eigs(corr_run, p);
      
    }
  }

  //autocorrelation(PhiAb_arr, avePhiAb, p.n_meas);
  
  //correlators(corr_run, corr_ave, idx, phi_cyl, avePhi*norm, p);
  //corr_eigs(corr_run, p);
  
  //free(corr_run);
  //free(corr_ave);

  free(s);
  free(phi);
}

void thermaliseSqL(double *phi, int *s, 
		   Param p, double &delta_mag_phi) {
  
  for(int iter = 0; iter < p.n_therm; iter++) {
    for(int i=0; i<p.n_wolff; i++) {
      swendsenWangUpdateSqL(phi, s, p, delta_mag_phi, iter);
      //wolffUpdateSqL(phi, s, p, delta_mag_phi, iter);
    }
    metropolisUpdateSqL(phi, s, p, delta_mag_phi, iter);
    if((iter+1)%p.n_skip == 0) cout<<"Therm sweep "<<iter+1<<endl;
  }
}
