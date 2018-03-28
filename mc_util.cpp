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
#include "data_io.h"
#include "mc_update_lr.h"
#include "mc_update_sr.h"

#ifdef USE_GPU
#include "gpu_mc_update_lr.cuh"
#endif

using namespace std;

extern int seed;
extern mt19937 rng;
extern uniform_real_distribution<double> unif;
extern int CACHE_LINE_SIZE;
extern int sze;

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

  long long time = 0.0;
  auto start1 = std::chrono::high_resolution_clock::now();
  cout<<"CPU malloc and init start..."<<endl;
  
  int x_len = p.S1/2 + 1;
  int t_len = p.Lt/2 + 1;
  int arr_len = x_len*t_len;
  
  //Long Range coupling array
  LR_couplings = (double*)malloc(arr_len*sizeof(double));
  //Arrays to hold spin and phi values.  
  s = (int*)malloc(p.surfaceVol*sizeof(int));
  phi = (double*)malloc(p.surfaceVol*sizeof(double));
  //Running correlation function arrays.
  run_corr_t = (double*)malloc((p.Lt/2+1)*sizeof(double));
  run_corr_s = (double*)malloc((p.S1/2+1)*sizeof(double));
  //Individual correlation functions.
  ind_corr_t = (double**)malloc((p.n_meas)*sizeof(double*));
  ind_corr_s = (double**)malloc((p.n_meas)*sizeof(double*));

  //Initialise Phi and spin fields.
  for(int i = 0;i < p.surfaceVol; i++) {
    phi[i] = 2.0*unif(rng) - 1.0;
    s[i] = (phi[i] > 0) ? 1 : -1;
    //cout<<"phi["<<i<<"] = "<<phi[i]<<" s["<<i<<"] = "<<s[i]<<endl;
  } 
  
  for(int j=0; j<t_len; j++){
    for(int i=0; i<x_len; i++){
      LR_couplings[i+j*x_len] = 0.0;
    }
  }
  
  for(int i=0; i<p.Lt/2+1; i++) {
    run_corr_t[i] = 0.0;
  }
  for(int i=0; i<p.S1/2+1; i++) {
    run_corr_s[i] = 0.0;
  }
  
  for(int i=0; i<p.n_meas; i++) {
    ind_corr_t[i] = (double*)malloc((p.Lt/2+1)*sizeof(double));
    ind_corr_s[i] = (double*)malloc((p.S1/2+1)*sizeof(double));
    
    for(int j=0; j<p.Lt/2+1; j++) ind_corr_t[i][j] = 0.0;
    for(int j=0; j<p.S1/2+1; j++) ind_corr_s[i][j] = 0.0;
  }
  
  auto elapsed1 = std::chrono::high_resolution_clock::now() - start1;
  time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
  cpu_added = (bool*)malloc(p.surfaceVol*sizeof(bool));
  cout<<"CPU malloc and init time = "<<time/(1.0e6)<<endl;

#ifdef USE_GPU

  //Device memory allocations.
  cout<<"CUDA malloc start..."<<endl;
  
  start1 = std::chrono::high_resolution_clock::now();
  //Allocate space on the GPU for the couplings and distances.
  cudaMalloc((void**) &gpu_LR_couplings, arr_len*sizeof(double));  
  cudaMalloc((void**) &gpu_phi, p.surfaceVol*sizeof(double));
  cudaMalloc((void**) &gpu_s, p.surfaceVol*sizeof(int));
  cudaMalloc((void**) &gpu_added, p.surfaceVol*sizeof(bool));  
  cudaMalloc((void**) &states, p.surfaceVol*sizeof(curandState_t));  
  cudaMalloc((void**) &gpu_rands, p.surfaceVol*sizeof(double));
  cudaMalloc((void**) &gpu_rands_aux, p.surfaceVol*sizeof(double));  
  cudaMalloc((void**) &gpu_result, p.surfaceVol*sizeof(double));
  
  elapsed1 = std::chrono::high_resolution_clock::now() - start1;
  time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
  cout<<"CUDA malloc time = "<<time/(1.0e6)<<endl;

  start1 = std::chrono::high_resolution_clock::now();
  GPU_initRand(p, seed, states);
  elapsed1 = std::chrono::high_resolution_clock::now() - start1;
  time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
  cout<<"CUDA rand states init time = "<<time/(1.0e6)<<endl;
  
  debug_arr1 = (double*)malloc(p.surfaceVol*sizeof(double));
  debug_arr2 = (double*)malloc(p.surfaceVol*sizeof(double));
  for(int j=0; j<p.S1; j++){
    for(int i=0; i<p.Lt; i++){
      debug_arr1[i+j*p.Lt] = 0.0;
      debug_arr2[i+j*p.Lt] = 0.0;
    }
  }
  
#endif
  
}

void MonteCarlo2DIsing::runSimulation(Param p) {

  //Create object to hold observable data.
  observables obs(p.n_meas);

  if(p.coupling_type != SR) {
    //Create LR couplings and kinetic term denomnators
    long long time = 0.0;
    auto start1 = std::chrono::high_resolution_clock::now();
    createLRcouplings(p);
    auto elapsed1 = std::chrono::high_resolution_clock::now() - start1;
    time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
    cout<<"Coupling creation time = "<<time/(1.0e6)<<endl;
  }
  
  thermalise(p);
    
  //reset timing variables.
  metro = 0.0;
  cluster = 0.0;
  for(int iter = p.n_therm + p.n_metro_cool; iter < p.n_therm + p.n_skip*p.n_meas + p.n_metro_cool; iter++) {
    
    updateIter(p, iter);
    
    //Take measurements.
    if((iter+1) % p.n_skip == 0) {

      //Timing
      int therm_iter = iter - (p.n_therm + p.n_metro_cool);
      cout<<"(Thermalised) Average time per Metro update   = "<<metro/(therm_iter*1.0e6)<<"s"<<endl;
      cout<<"(Thermalised) Average time per cluster update = "<<cluster/(therm_iter*p.n_cluster*1.0e6)<<"s"<<endl;

      //copy GPU arrays to the CPU FIXME: do all measurements on the GPU
      
      GPU_copyArraysToHost(p);      
      
      //update running average, dump to stdout. 
      //'meas' is ++ incremented in this function.
      measure(obs, phi, s, LR_couplings, meas, p);

      norm = 1.0/(meas);

      //Visualisation tool
      //visualiserSqr(phi, obs.avePhiAb*norm, p);
      
      //Calculate correlaton functions and update the average.
      /*
      correlatorsImp(ind_corr_t, meas-1, run_corr_t, true,  
		     phi, obs.avePhi*norm, s, p);
      correlatorsImp(ind_corr_s, meas-1, run_corr_s, false, 
		     phi, obs.avePhi*norm, s, p);
      */
      correlators(ind_corr_t, meas-1, run_corr_t, true,  
		  phi, obs.avePhi*norm, p);
      correlators(ind_corr_s, meas-1, run_corr_s, false, 
		  phi, obs.avePhi*norm, p);
      
      //Jacknife and dump the data
      if(meas%10 == 0) {	
	writeObservables(ind_corr_t, run_corr_t, ind_corr_s,
			 run_corr_s, meas, obs, p);
      }      
    }
  }
  
  free(run_corr_t);
  free(run_corr_s);
  free(s);
  free(phi);
  
}

void MonteCarlo2DIsing::updateIter(Param p, int iter) {
  
    auto start = std::chrono::high_resolution_clock::now();    
    //Cluster steps.
    if(p.useWolff)  wolffUpdate(p, iter);
    else swendsenWangUpdate(p, iter);
    
    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    cluster += std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

    //Metro step.
    start = std::chrono::high_resolution_clock::now();
    metropolisUpdate(p, iter);
    elapsed = std::chrono::high_resolution_clock::now() - start;
    metro += std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
}

void MonteCarlo2DIsing::thermalise(Param p) {
  
  cout<<"Begin Metropolis cool down."<<endl;
  //Perform n_therm pure metropolis hits during the hot phase.
  auto start = std::chrono::high_resolution_clock::now();        
  for(int iter = 0; iter < p.n_metro_cool; iter++) {
    metropolisUpdate(p, iter);  

    if((iter)%p.n_skip == 0 && iter > 0){
      auto elapsed = std::chrono::high_resolution_clock::now() - start;
      metro = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
      cout<<"Metro cool down iter "<<iter<<" average time per (hot phase) update = "<<metro/((iter+1)*1.0e6)<<"s"<<endl<<endl;
      //exit(0);
    }
  }

  auto elapsed = std::chrono::high_resolution_clock::now() - start;
  metro = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
  cout<<"Metro cool down complete. "<<p.n_metro_cool<<" iters at "<<metro/((p.n_metro_cool)*1.0e6)<<"s per (hot phase) iter."<<endl<<endl;

  //Reset timing variables
  metro = 0.0;
  cluster = 0.0;

  //Cluster/Metro steps
  for(int iter = p.n_metro_cool; iter < p.n_therm + p.n_metro_cool; iter++) {
    updateIter(p, iter);
    if( (iter+1)%p.n_skip == 0 ) {
      cout<<"Average time per Metro update   = "<<metro/((iter - p.n_metro_cool + 1)*1.0e6)<<"s"<<endl;  
      cout<<"Average time per cluster update = "<<cluster/((iter - p.n_metro_cool + 1)*p.n_cluster*1.0e6)<<"s"<<endl;  
    }      
  }
  cout<<endl<<"Thermalisaton complete."<<endl<<endl;  
}

//------------------------------------------------//
// Wrappers for the different LR and SR functions //
//------------------------------------------------//
void MonteCarlo2DIsing::wolffUpdate(Param p, int iter) {

  for(int i=0; i<p.n_cluster; i++) {
    if(p.coupling_type == SR) {
      wolffUpdateSR(phi, s, p, iter);
    }
    else {
      if(p.useGPUCluster) {
	if(!p.useGPUMetro) {
	  //We must copy from the host to the device
	  GPU_copyArraysToDevice(p);
	}
	
	int rand_site = int(unif(rng) * p.surfaceVol);	
	GPU_wolffUpdateLR(p, rand_site, iter, i);
	
	if(!p.useGPUMetro) {
	  //We must copy from the device to the host
	  GPU_copyArraysToHost(p);
	}	
      } else {
	wolffUpdateLR(phi, s, p, LR_couplings, iter, i);
      }
    }
  }
}

void MonteCarlo2DIsing::swendsenWangUpdate(Param p, int iter){
  
  if(p.coupling_type == SR) 
    swendsenWangUpdateSR(phi, s, p, iter); 
  else
    swendsenWangUpdateLR(phi, s, p, LR_couplings, iter); 
}

void MonteCarlo2DIsing::metropolisUpdate(Param p, int iter) {
  
  if(p.coupling_type == SR) {
    metropolisUpdateSR(phi, s, p, iter);
  }
  else {
    if(p.useGPUMetro) {
      if(!p.useGPUCluster) {
	//We must copy from the host to the device
	GPU_copyArraysToDevice(p);
      }
      
      GPU_metropolisUpdateLR(p, iter);
      if(p.doMetroCheck) {
	//Use the same fields, check against the CPU.
	metropolisUpdateLR(p, iter);
      }
      
      if(!p.useGPUCluster) {
	//We must copy from the device to the host
	GPU_copyArraysToHost(p);
      }
    } else {
      metropolisUpdateLR(p, iter);
    }
  }
}

//-------------------------------------------//
// Container class for observable quantities //
//-------------------------------------------//
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
  int S1 = p.S1;
  int Lt = p.Lt;
  int vol = S1*Lt;
  
  int x_len = S1/2 + 1;
    
  double phi_sq;
  double phi;
  double lambda_p = 0.25*p.lambda;
  double musqr_p  = 0.50*p.musqr;
  
  for (int i = 0; i < vol; i++)
    if (s[i] * phi_arr[i] < 0)
      printf("ERROR s and phi NOT aligned (actionPhi SqNL) ! \n");
  
  for (int i = 0; i < vol; i++) {    
    phi = phi_arr[i];
    phi_sq = phi*phi;
    
    //PE terms
    PE += lambda_p * phi_sq*phi_sq;
    PE += musqr_p  * phi_sq;

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
	    
	    //FIXME: Here we are overcounting by a factor of two, hence the 0.25 coefficient.
	    val = 0.25*((phi - phi_arr[j+k])*(phi - phi_arr[j+k])*
			LR_couplings[dx+dt*x_len]*LR_couplings[dx+dt*x_len]);
	    
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
	
	KE += 0.25*((phi - phi_arr[j])*(phi - phi_arr[j])*
		    LR_couplings[dx+dt*x_len]*LR_couplings[dx+dt*x_len]);
      }
    }
#endif
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
  for(int i = 0;i < p.surfaceVol; i++) {
    obs.MagPhi += phi[i];
    //cout<<phi[i]<<endl;
  }
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
		      int idx, observables obs, Param p){

  double norm = 1.0/idx;
  double inv_vol = 1.0/(p.Lt*p.S1);
  
  //The observable (spec heat, suscep, etc) are used only as a guide
  //to discern criticality, The real quantities of interest are
  //The critical exponents, especally \eta, and its dependence
  //on \sigma, hence we MUST have proper error estimates of the
  //correlation function values. Hardcoded to dump result every 10th
  //measurement, with a jk block of 5.
  double *jk_err_t = (double*)malloc((p.Lt/2+1)*sizeof(double));
  jackknife(ind_corr_t, run_corr_t, jk_err_t, 5, idx, p.Lt/2+1, p); 
  ofstream filet("correlators_t.dat");
  for(int i=1; i<p.Lt/2; i++) {
    filet<<i<<" "<<run_corr_t[i]*inv_vol*norm;
    filet<<" "<<jk_err_t[i]*inv_vol<<endl;
  }
  filet.close();
  free(jk_err_t);
  
  double *jk_err_s = (double*)malloc((p.S1/2+1)*sizeof(double));
  jackknife(ind_corr_s, run_corr_s, jk_err_s, 5, idx, p.S1/2+1, p);
  ofstream files("correlators_s.dat");
  for(int i=1; i<p.S1/2; i++) {
    files<<i<<" "<<run_corr_s[i]*inv_vol*norm;
    files<<" "<<jk_err_s[i]*inv_vol<<endl;
  }
  files.close();
  free(jk_err_s);
  
  ofstream file_obs("observables.dat");
  for(int i=0; i<idx; i++) {
    file_obs<<i<<" "<<obs.Suscep[i]<<" "<<obs.SpecHeat[i]<<" "<<obs.Binder[i]<<endl;
  }
  file_obs.close();
  
}

inline double Coupling(double dt, double dth, Param p) {
 
  dt  *= M_PI/p.S1;
  dth *= M_PI/p.S1;  

  return 1.0/pow(cosh(dt) - cos(dth), 1+p.sigma/2);

}

void MonteCarlo2DIsing::createLRcouplings(Param p) {
  
  double sigma = p.sigma;
  int S1 = p.S1;
  int Lt = p.Lt;
  int x_len = S1/2 + 1;
  int t_len = Lt/2 + 1;
  int arr_len = x_len*t_len;

  double couplingNormTheta = 1.0/Coupling(0, 1, p);
  double couplingNormTime = 1.0/Coupling(1, 0, p);

  for(int j=0; j<t_len; j++){
    for(int i=0; i<x_len; i++){
      
      //Associate dx with i, dt with j. dx runs fastest
      int idx = i + j*x_len;      
      
      if(p.usePowLaw) LR_couplings[idx] = pow(sqrt(i*i + j*j),-(2+sigma));
      else  LR_couplings[idx] = couplingNormTheta*Coupling((double)j, (double)i, p);
      
    }
  }
  
  LR_couplings[0] = couplingNormTheta;  
  
#ifdef USE_GPU
  cudaMemcpy(gpu_LR_couplings, LR_couplings,
	     arr_len*sizeof(double), cudaMemcpyHostToDevice);
#endif
  
}

