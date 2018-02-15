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

#include "hyp_util.h"
#include "mc_util.h"
#include "data_proc.h"
#include "mc_update_lr.h"
#include "mc_update_sr.h"

using namespace std;

extern int seed;
extern mt19937 rng;
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

//---------------------------------------//
//Square 2D Ising Monte Carlo base class //
//---------------------------------------//
MonteCarlo2DIsing::MonteCarlo2DIsing(Param p) {

  //Long Range coupling array
  LR_couplings = (double*)malloc(p.surfaceVol*p.surfaceVol*sizeof(double));
  for(int i=0; i<p.surfaceVol; i++){
    for(int j=0; j<p.surfaceVol; j++){
      LR_couplings[i+j*p.surfaceVol] = 0.0;;
    }
  }

  //Create arrays to hold spin and phi values.  
  s = (int*)malloc(p.surfaceVol*sizeof(int));
  phi = (double*)malloc(p.surfaceVol*sizeof(double));
  for(int i = 0;i < p.surfaceVol; i++) {
    phi[i] = 2.0*unif(rng) - 1.0;
    s[i] = (phi[i] > 0) ? 1:-1;
  } 

  //Running correlation function arrays.
  run_corr_t = (double*)malloc((p.Lt/2+1)*sizeof(double));
  for(int i=0; i<p.Lt/2+1; i++) {
    run_corr_t[i] = 0.0;
  }
  run_corr_s = (double*)malloc((p.S1/2+1)*sizeof(double));
  for(int i=0; i<p.S1/2+1; i++) {
    run_corr_s[i] = 0.0;
  }

  //Individual correlation functions.
  ind_corr_t = (double**)malloc((p.n_meas)*sizeof(double*));
  for(int i=0; i<p.n_meas; i++) {
    ind_corr_t[i] = (double*)malloc((p.Lt/2+1)*sizeof(double));
    for(int j=0; j<p.Lt/2+1; j++) {
      ind_corr_t[i][j] = 0.0;
    }
  }
  ind_corr_s = (double**)malloc((p.n_meas)*sizeof(double*));
  for(int i=0; i<p.n_meas; i++) {
    ind_corr_s[i] = (double*)malloc((p.S1/2+1)*sizeof(double));
    for(int j=0; j<p.S1/2+1; j++) {
      ind_corr_s[i][j] = 0.0;
    }
  }

  //The correlation functions must be weighted by how frequently
  //the i,j separation occurs  
  corr_norms = (int**)malloc((p.S1/2+1)*sizeof(int*));
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
}

void MonteCarlo2DIsing::runSimulation(Param p) {

  observables obs(p.n_meas);

  createLRcouplings(LR_couplings, p);
  
  thermalise(phi, s, p, obs);
  
  for(int iter = p.n_therm; iter < p.n_therm + p.n_skip*p.n_meas; iter++) {

    updateIter(obs, p, iter);
    
    //Take measurements.
    if((iter+1) % p.n_skip == 0) {
      
      //update running average, dump to stdout. 
      //'meas' is ++ incremented in this function.
      measure(obs, phi, s, LR_couplings, meas, p);

      norm = 1.0/(meas);

      //Visualisation tool
      //visualiserSqr(phi, avePhiAb*norm, p);
      
      //Calculate correlaton functions and update the average.
      correlators(ind_corr_t, meas-1, run_corr_t, true,  
		  phi, obs.avePhi*norm, p);
      correlators(ind_corr_s, meas-1, run_corr_s, false, 
		  phi, obs.avePhi*norm, p);
      
      //Jacknife and dump the data
      if(meas%10 == 0) {	
	writeObservables(ind_corr_t, run_corr_t, ind_corr_s, run_corr_s,
			 corr_norms, meas, obs, p);
      }      
    }
  }
  
  free(run_corr_t);
  free(run_corr_s);
  free(corr_norms);
  free(s);
  free(phi);
}

void MonteCarlo2DIsing::updateIter(observables obs, Param p, int iter) {
  
    auto start = std::chrono::high_resolution_clock::now();    
    //Cluster steps.
    for(int i=0; i<p.n_cluster; i++) {
      if(p.useWolff) {
	wolffUpdate(phi, s, p, obs.delta_mag_phi, iter);  
      }
      else { 
	swendsenWangUpdate(phi, s, p, obs.delta_mag_phi, iter);
      }
    }
    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    cluster += std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

    //Metro step.
    start = std::chrono::high_resolution_clock::now();
    metropolisUpdate(phi, s, p, obs.delta_mag_phi, iter);
    elapsed = std::chrono::high_resolution_clock::now() - start;
    metro += std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

}

void MonteCarlo2DIsing::thermalise(double *phi, int *s, Param p, 
				   observables obs) {
  
  long long metro = 0.0;
  long long cluster = 0.0;

  cout<<"Begin Metropolis cool down."<<endl;
  //Perform n_therm pure metropolis hits during the hot phase.
  auto start = std::chrono::high_resolution_clock::now();        
  for(int iter = 0; iter < p.n_therm; iter++) {
    metropolisUpdate(phi, s, p, obs.delta_mag_phi, iter);  
    if((iter+1)%p.n_skip == 0) cout<<"Metro cool down iter "<<iter+1<<endl;
  }
  auto elapsed = std::chrono::high_resolution_clock::now() - start;
  metro = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
  cout<<"Metro cool down complete. "<<p.n_therm<<" iters in "<<metro/(1.e6)<<"s"<<endl;
  metro = 0.0;

  //Cluster/Metro steps
  for(int iter = 0; iter < p.n_therm; iter++) {
    updateIter(obs, p, iter);
  }
  cout<<endl<<"Thermalisaton complete."<<endl<<endl;  
}

//------------------------------------------------//
// Wrappers for the different LR and SR functions //
//------------------------------------------------//
void MonteCarlo2DIsing::wolffUpdate(double *phi, int *s, Param p, 
				    double &delta_mag_phi, int iter) {
  
  if(p.coupling_type == SR) 
    wolffUpdateSR(phi, s, p, delta_mag_phi, iter); 
  else
    wolffUpdateLR(phi, s, p, LR_couplings, delta_mag_phi, iter); 
}

void MonteCarlo2DIsing::swendsenWangUpdate(double *phi, int *s, Param p, 
					   double &delta_mag_phi, int iter){
  
  if(p.coupling_type == SR) 
    swendsenWangUpdateSR(phi, s, p, delta_mag_phi, iter); 
  else
    swendsenWangUpdateLR(phi, s, p, LR_couplings, delta_mag_phi, iter); 
}

void MonteCarlo2DIsing:: metropolisUpdate(double *phi, int *s, Param p, 
					  double &delta_mag_phi, int iter) {
  
  if(p.coupling_type == SR) 
    metropolisUpdateSR(phi, s, p, delta_mag_phi, iter); 
  else
    metropolisUpdateLR(phi, s, p, LR_couplings, delta_mag_phi, iter); 
}

//-----------------------------------------------------------//
// Container class for keeping hold of observable quantities //
//-----------------------------------------------------------//
observables::observables(int meas) {

  n_meas = meas;
  
  E_arr = new double[n_meas];
  E2_arr = new double[n_meas];
  PhiAb_arr = new double[n_meas];
  Phi_arr = new double[n_meas];
  Phi2_arr = new double[n_meas];
  Phi4_arr = new double[n_meas];
  Suscep = new double[n_meas];
  SpecHeat = new double[n_meas];
  Binder = new double[n_meas];
  for(int i=0; i<n_meas; i++) {
    E_arr[i]     = 0.0;
    E2_arr[i]    = 0.0;
    PhiAb_arr[i] = 0.0;
    Phi_arr[i]   = 0.0;
    Phi2_arr[i]  = 0.0;
    Phi4_arr[i]  = 0.0;
    Suscep[i]    = 0.0;
    SpecHeat[i]  = 0.0;
    Binder[i]    = 0.0;
  }

  tmpE     = 0.0;  
  aveE     = 0.0;
  aveKE    = 0.0;
  avePE    = 0.0;
  aveE2    = 0.0;
  avePhiAb = 0.0;
  avePhi   = 0.0;
  avePhi2  = 0.0;
  avePhi4  = 0.0;
  MagPhi   = 0.0;  
   
}

double actionLR(double *phi_arr, int *s, Param p,
		double *LR_couplings, 
		double &KE, double &PE) {  
  
  KE = 0.0;
  PE = 0.0;
  double phi_sq;
  double phi;
  double lambda_p = 0.25*p.lambda;
  double musqr_p  = 0.50*p.musqr;

  for (int i = 0; i < p.surfaceVol; i++)
    if (s[i] * phi_arr[i] < 0)
      printf("ERROR s and phi NOT aligned (actionPhi SqNL) ! \n");
  
  for (int i = 0; i < p.surfaceVol; i++) {    
    phi = phi_arr[i];
    phi_sq = phi*phi;

    //PE terms
    PE += lambda_p * phi_sq*phi_sq;
    PE += musqr_p  * phi_sq;

    //KE terms
    for(int j=0; j<p.surfaceVol; j++) {      
      //Here we are overcounting by a factor of two, hence the 0.25
      //coefficient.
      //FIXME
      KE += 0.25*(phi - phi_arr[j])*(phi - phi_arr[j])*LR_couplings[i+j*p.surfaceVol]*LR_couplings[i+j*p.surfaceVol];
    }
  }
  return PE + KE;
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


void measure(observables &obs, double *phi, int *s, 
	     double *LR_couplings, int &idx, Param p) {

  double rhoVol  = 1.0/p.surfaceVol;
  double rhoVol2 = rhoVol*rhoVol;

  if(p.coupling_type == SR) {
    obs.tmpE = actionSR(phi, s, p, obs.KE, obs.PE);
  } else {
    obs.tmpE = actionLR(phi, s, p, LR_couplings, obs.KE, obs.PE);
  }

  obs.aveKE += rhoVol*obs.KE;
  obs.avePE += rhoVol*obs.PE;
  obs.aveE  += rhoVol*obs.tmpE;
  obs.aveE2 += rhoVol2*obs.tmpE*obs.tmpE;
  
  obs.MagPhi = 0.0;
  for(int i = 0;i < p.surfaceVol; i++) obs.MagPhi += phi[i];
  obs.MagPhi *= rhoVol;
  
  obs.avePhiAb  += abs(obs.MagPhi);
  obs.avePhi    += obs.MagPhi;
  obs.avePhi2   += obs.MagPhi*obs.MagPhi;
  obs.avePhi4   += obs.MagPhi*obs.MagPhi*obs.MagPhi*obs.MagPhi;
  
  obs.E_arr[idx]     = rhoVol*obs.tmpE;
  obs.E2_arr[idx]    = rhoVol*obs.tmpE*obs.tmpE;
  obs.PhiAb_arr[idx] = abs(obs.MagPhi);
  obs.Phi_arr[idx]   = obs.MagPhi;
  obs.Phi2_arr[idx]  = obs.MagPhi*obs.MagPhi;
  obs.Phi4_arr[idx]  = obs.MagPhi*obs.MagPhi*obs.MagPhi*obs.MagPhi;
  
  idx++;
  
  cout<<setprecision(8);
  double norm = 1.0/(idx);

  obs.Suscep[idx-1]   = (obs.avePhi2*norm-pow(obs.avePhiAb*norm,2))/rhoVol;
  obs.SpecHeat[idx-1] = (obs.aveE2*norm-pow(obs.aveE*norm,2))/rhoVol;
  obs.Binder[idx-1]   = 1.0-obs.avePhi4/(3.0*obs.avePhi2*obs.avePhi2*norm);
  
  //Dump to stdout
  cout<<"Measurement "<<idx<<endl;
  cout<<"Ave Energy= "<<obs.aveE*norm<<endl;
  cout<<"Ave KE    = "<<obs.aveKE*norm<<endl;
  cout<<"Ave PE    = "<<obs.avePE*norm<<endl;
  cout<<"Ave |phi| = "<<obs.avePhiAb*norm<<endl;
  cout<<"Ave phi   = "<<obs.avePhi*norm<<endl;
  cout<<"Ave phi^2 = "<<obs.avePhi2*norm<<endl;
  cout<<"Ave phi^4 = "<<obs.avePhi4*norm<<endl;
  cout<<"Suscep    = "<<obs.Suscep[idx-1]<<endl;
  cout<<"Spec Heat = "<<obs.SpecHeat[idx-1]<<endl;
  cout<<"Binder    = "<<obs.Binder[idx-1]<<endl;
  
}

void writeObservables(double **ind_corr_t, double *run_corr_t, 
		      double **ind_corr_s, double *run_corr_s,
		      int **corr_norms, int idx, observables obs, Param p){

  double norm = 1.0/idx;
  
  //The observable (spec heat, suscep, etc) are used only as a guide
  //to discern criticality, The real quantities of interest are
  //The critical exponents, especally \eta, and its dependence
  //on \sigma, hence we MUST have proper error estimates of the
  //correlation function values. Hardcoded to dump result every 10th
  //measurement, with a jk block of 5.
  double *jk_err_t = (double*)malloc((p.Lt/2+1)*sizeof(double));
  jackknife(ind_corr_t, run_corr_t, jk_err_t, 5, idx, p.Lt/2+1, p); 
  ofstream filet("correlators_t.dat");
  for(int i=0; i<p.Lt/2+1; i++) {
    filet<<i<<" "<<(run_corr_t[i]/corr_norms[0][i])*norm;
    filet<<" "<<(jk_err_t[i]/corr_norms[0][i])<<endl;
  }
  filet.close();
  free(jk_err_t);
  
  double *jk_err_s = (double*)malloc((p.S1/2+1)*sizeof(double));
  jackknife(ind_corr_s, run_corr_s, jk_err_s, 5, idx, p.S1/2+1, p);
  ofstream files("correlators_s.dat");
  for(int i=0; i<p.S1/2+1; i++) {
    files<<i<<" "<<(run_corr_s[i]/corr_norms[i][0])*norm;
    files<<" "<<(jk_err_s[i]/corr_norms[i][0])<<endl;
  }
  files.close();
  free(jk_err_s);
  
  ofstream file_obs("observables.dat");
  for(int i=0; i<idx; i++) {
    file_obs<<i<<" "<<obs.Suscep[i]<<" "<<obs.SpecHeat[i]<<" "<<obs.Binder[i]<<endl;
  }
  file_obs.close();
  
}

void createLRcouplings(double *LR_couplings, Param p) {
  
  double sigma = p.sigma;
  double r_x = 0.0;
  double r_y = 0.0;

  if(p.coupling_type == SR) {
    //do nothing, no LR interactions
  } else {    
    if(p.usePowLaw) {
      //Square lattice, 1/r^{d+\sigma} coupling
      for(int i=0; i<p.surfaceVol; i++) {
	for(int j=0; j<p.surfaceVol; j++) {

	  //Index divided by circumference, using the int floor feature/bug,
	  //gives the timeslice index.
	  int t1 = i / p.S1;
	  int t2 = j / p.S1;
	  int dt = abs(t2-t1) > p.Lt/2 ? p.Lt - abs(t2-t1) : abs(t2-t1);
		  
	  //The index modulo the circumference gives the spatial index.
	  int x1 = i % p.S1;
	  int x2 = j % p.S1;            
	  double dx = abs(x2-x1) > p.S1/2 ? p.S1-abs(x2-x1) : abs(x2-x1);      
		  
	  // r_x = abs(i%p.S1 - j%p.S1);
	  // if(r_x > p.S1/2) r_x = p.S1 - r_x;
	  
	  // r_y = abs(i/p.Lt - j/p.Lt);
	  // if(r_y > p.Lt/2) r_y = p.Lt - r_y;
	  
	  if(i != j) {
	    LR_couplings[i+j*p.surfaceVol] = pow(sqrt(dx*dx + dt*dt),-(2+sigma));
	    
	    //Temporal
	    if(i==0 && j%p.S1==0) cout<<"time "<<j/p.S1<<" "<<LR_couplings[i+j*p.surfaceVol]<<endl;
	    //Spatial
	    if(i==0 && j<p.S1) cout<<"space "<<j<<" "<<LR_couplings[i+j*p.surfaceVol]<<endl;
	  }
	}
      }
    } else {
      //Square lattice, user defined definition of coupling.
      LRAdSCouplings(LR_couplings, p);
    }
  }
}

inline double Coupling(double dt, double dth, Param p) {
  
  //return 1.0/(cosh(dt) - cos(dth));
  return exp(-(cosh(dt) - cos(dth))*(1+sqrt(1+p.msqr)));
}

void LRAdSCouplings(double *LR_couplings, Param &p){

  //Set the t scale so that LR(0,dt) \approx LR(dth,0). Do this for
  //a separation of one lattice spacing.

  double couplingNorm = 1.0/Coupling(0, M_PI/p.S1, p);
  p.t_scale = M_PI/p.S1;
  
  for(int i=0; i<p.surfaceVol; i++){
    for(int j=0; j<p.surfaceVol; j++){
      //Index divided by circumference, using the int floor feature/bug,
      //gives the timeslice for each index.
      int t1 = i / p.S1;
      int t2 = j / p.S1;
      double dt = abs(t2-t1) > p.Lt/2 ? p.Lt - abs(t2-t1) : abs(t2-t1);
      dt *= p.t_scale;
      
      //The index modulo the circumference gives the 'angle'
      //for each index.
      int th1 = i % p.S1;
      int th2 = j % p.S1;            
      double dth = abs(th2-th1) > p.S1/2 ? p.S1-abs(th2-th1) : abs(th2-th1);      
      dth *= M_PI/p.S1;
      
      if(i!=j) {
	LR_couplings[i+j*p.surfaceVol] = pow(couplingNorm*Coupling(dt, dth, p), (2+p.sigma)/2);
	
	//Temporal
	if(i==0 && j%p.S1==0) cout<<0<<" "<<j/p.S1<<" "<<dt<<" "<<LR_couplings[i+j*p.surfaceVol]<<endl;
	//Spatial
	if(i==0 && j<p.S1) cout<<0<<" "<<j<<" "<<dth<<" "<<LR_couplings[i+j*p.surfaceVol]<<endl;
      }
    }
  }
}

