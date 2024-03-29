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

#include "phiFourth2D.h"
#include "data_proc.h"
#include "data_io.h"
#include "mcPhiFourth2D.h"
#include "gpuMC.cuh"

using namespace std;

extern int seed;
extern mt19937 rng;
extern uniform_real_distribution<double> unif;
extern int CACHE_LINE_SIZE;
extern int sze;

//---------------------//
// 2D Phi Fourth class //
//---------------------//
PhiFourth2D::PhiFourth2D(Param p) : Ising2D(p) {
  
  long long time = 0.0;
  auto start1 = std::chrono::high_resolution_clock::now();
  cout<<"CPU malloc and init start..."<<endl;

  int vol = p.surfaceVol;
  int arr_len = (p.S1/2 + 1)*(p.Lt/2 + 1);
  
  //Arrays to hold phi values.  
  phi = (double*)malloc(vol*sizeof(double));
  phi_cpy = (double*)malloc(vol*sizeof(double));

  //Initialise Phi and spin fields.
  for(int i = 0;i < vol; i++) {
    phi[i] = 2.0*unif(rng) - 1.0;
    s[i] = (phi[i] > 0) ? 1 : -1;
    //cout<<"phi["<<i<<"] = "<<phi[i]<<" s["<<i<<"] = "<<s[i]<<endl;
  }
  cpu_added = (bool*)malloc(vol*sizeof(bool));

  //FT arrays.
  ind_ft_corr_phi_phi3 = (double**)malloc(((p.Lt/2+1)*3)*sizeof(double*));
  run_ft_corr_phi_phi3 = (double*)malloc(3*(p.Lt/2+1)*sizeof(double));

  ind_ft_corr_phi3_phi3 = (double**)malloc(((p.Lt/2+1)*3)*sizeof(double*));
  run_ft_corr_phi3_phi3 = (double*)malloc(3*(p.Lt/2+1)*sizeof(double));

  //ft correlation function
  for(int i=0; i<(p.Lt/2 +1)*3; i++) {
    ind_ft_corr_phi_phi3[i] = (double*)malloc(p.n_meas*sizeof(double));
    ind_ft_corr_phi3_phi3[i] = (double*)malloc(p.n_meas*sizeof(double));
    run_ft_corr_phi_phi3[i] = 0.0;
    run_ft_corr_phi3_phi3[i] = 0.0;
    for(int j=0; j<p.n_meas; j++) {
      ind_ft_corr_phi_phi3[i][j] = 0.0;
      ind_ft_corr_phi3_phi3[i][j] = 0.0;     
    }
  }
  
  //Running phiphi^3 and phi3phi3 correlation function arrays.
  run_corr_phi_phi3 = (double*)malloc(arr_len*sizeof(double));
  run_corr_phi3_phi3 = (double*)malloc(arr_len*sizeof(double));
  for(int i=0; i<arr_len; i++) {
    run_corr_phi_phi3[i] = 0.0;
    run_corr_phi3_phi3[i] = 0.0;
  }
  
  //Individual phiphi^3 and phi3phi3 correlation functions.
  ind_corr_phi_phi3 = (double**)malloc(p.n_meas*sizeof(double*));
  ind_corr_phi3_phi3 = (double**)malloc(p.n_meas*sizeof(double*));
  for(int i=0; i<p.n_meas; i++) {
    ind_corr_phi_phi3[i] = (double*)malloc(arr_len*sizeof(double));
    ind_corr_phi3_phi3[i] = (double*)malloc(arr_len*sizeof(double));
    for(int j=0; j<arr_len; j++) {
      ind_corr_phi_phi3[i][j] = 0.0;
      ind_corr_phi3_phi3[i][j] = 0.0;
    }
  }
    
  auto elapsed1 = std::chrono::high_resolution_clock::now() - start1;
  time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
  cout<<"CPU malloc and init time = "<<time/(1.0e6)<<endl;
  
#ifdef USE_GPU

  //Device memory allocations.
  cout<<"CUDA malloc start..."<<endl;
  
  start1 = std::chrono::high_resolution_clock::now();
  //Allocate space on the GPU.
  cudaMalloc((void**) &gpu_phi, vol*sizeof(double));
  cudaMalloc((void**) &gpu_phi_cpy, vol*sizeof(double));
  
  elapsed1 = std::chrono::high_resolution_clock::now() - start1;
  time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
  cout<<"CUDA malloc time = "<<time/(1.0e6)<<endl;

#endif  
}


void PhiFourth2D::runSimulation(Param p) {

  //Create object to hold observable data.
  observables obs(p.n_meas);

  //Counts number of measurements.
  int meas = 0;

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
#ifdef USE_GPU      
      //GPU_copyArraysToHost(p);      
#endif
      //update running average, dump to stdout. 
      //'meas' is ++ incremented in this function.
      measure(obs, meas, p);
      
      norm = 1.0/(meas);

      //Visualisation tool
      //visualiserSqr(phi, obs.avePhiAb*norm, p);
      
      //Calculate correlaton functions and update the average.
      //int rand_site = int(unif(rng) * p.surfaceVol);
#ifdef USE_GPU
      //GPU_correlatorsImpWolff(rand_site, meas-1, obs.avePhi*norm, p);
#endif
     
      long long time = 0.0;
      auto start1 = std::chrono::high_resolution_clock::now();
      correlatorsImpWolff(meas-1, obs.avePhi*norm, p);
      
      auto elapsed1 = std::chrono::high_resolution_clock::now() - start1;
      time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
      cout<<" PhiPhi (Improved) Correlation Function Calculation time = "<<time/(1.0e6)<<endl;
      
      start1 = std::chrono::high_resolution_clock::now();
      correlatorsPhi3(meas-1, obs.avePhi*norm, obs.phi3Ave*norm, p);
      
      elapsed1 = std::chrono::high_resolution_clock::now() - start1;
      time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
      cout<<"PhiPhi3 and Phi3Phi3 Correlation Function Calculation time = "<<time/(1.0e6)<<endl;
      
      time = 0.0;
      start1 = std::chrono::high_resolution_clock::now();

      //FT the correlation function values.
      FTcorrelation(ind_ft_corr, run_ft_corr, ind_corr, norm_corr, meas-1, p);
      FTcorrelation(ind_ft_corr_phi_phi3, run_ft_corr_phi_phi3, ind_corr_phi_phi3, norm_corr, meas-1, p);
      FTcorrelation(ind_ft_corr_phi3_phi3, run_ft_corr_phi3_phi3, ind_corr_phi3_phi3, norm_corr, meas-1, p);
      
      elapsed1 = std::chrono::high_resolution_clock::now() - start1;
      time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
      cout<<"Correlation Function FT time = "<<time/(1.0e6)<<endl;
      
      if(meas%p.n_write == 0) {	
	time = 0.0;
	start1 = std::chrono::high_resolution_clock::now();
	//Jacknife and dump the data
	writeObservables(ind_corr, run_corr, norm_corr,
			 ind_ft_corr, run_ft_corr, meas, obs, p);
	
	writePhi3(ind_corr_phi_phi3, run_corr_phi_phi3,
		  ind_corr_phi3_phi3, run_corr_phi3_phi3,
		  ind_ft_corr_phi_phi3, run_ft_corr_phi_phi3,
		  ind_ft_corr_phi3_phi3, run_ft_corr_phi3_phi3,
		  norm_corr, meas, obs, p);
	
	elapsed1 = std::chrono::high_resolution_clock::now() - start1;
	time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
	cout<<"Correlation, observables, and FT Jackknife time = "<<time/(1.0e6)<<endl;
	
	time = 0.0;
	start1 = std::chrono::high_resolution_clock::now();
	//Calculate the autocorrelation of |phi|
	autocorrelation(obs.PhiAb_arr, obs.avePhiAb*norm, meas, auto_corr);
	
	elapsed1 = std::chrono::high_resolution_clock::now() - start1;
	time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
	cout<<"Autocorrelation calculation time = "<<time/(1.0e6)<<endl;      
      }
    }
  }
}

void PhiFourth2D::updateIter(Param p, int iter) {
  
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

void PhiFourth2D::thermalise(Param p) {
  
  cout<<"Begin Metropolis cool down."<<endl;
  //Perform n_therm pure metropolis hits during the hot phase.
  auto start = std::chrono::high_resolution_clock::now();        
  for(int iter = 0; iter < p.n_metro_cool; iter++) {
    metropolisUpdate(p, iter);  
    
    if((iter)%p.n_skip == 0 && iter > 0){
      auto elapsed = std::chrono::high_resolution_clock::now() - start;
      metro = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
      cout<<"Metro cool down iter "<<iter<<" average time per (hot phase) update = "<<metro/((iter+1)*1.0e6)<<"s"<<endl<<endl;
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
void PhiFourth2D::wolffUpdate(Param p, int iter) {

  if(p.coupling_type == SR) {
    for(int i=0; i<p.n_cluster; i++) wolffUpdateSR(phi, s, p, iter);
  }
  else {
    if(p.useGPUCluster) {
#ifdef USE_GPU
      if(!p.useGPUMetro) GPU_copyArraysToDevice(p);
      
      for(int i=0; i<p.n_cluster; i++) {
	int rand_site = int(unif(rng) * p.surfaceVol);	
	GPU_wolffUpdateLR(p, rand_site, iter, i);
      }
      
      if(!p.useGPUMetro) GPU_copyArraysToHost(p);
#else
      cout<<"GPU Wolf not built"<<endl;
#endif
    } else {
      for(int i=0; i<p.n_cluster; i++) wolffUpdateLR(phi, s, p, LR_couplings, iter, i);
    }
  }
}

void PhiFourth2D::swendsenWangUpdate(Param p, int iter){
  
  if(p.coupling_type == SR) 
    swendsenWangUpdateSR(phi, s, p, iter); 
  else
    swendsenWangUpdateLR(p, iter); 
}

void PhiFourth2D::metropolisUpdate(Param p, int iter) {
  
  if(p.coupling_type == SR) {
    metropolisUpdateSR(phi, s, p, iter);
  }
  else {
    if(p.useGPUMetro) {
#ifdef USE_GPU      
      if(!p.useGPUCluster) GPU_copyArraysToDevice(p);
      GPU_metropolisUpdateLR(p, iter);
      if(p.doMetroCheck) {
	//Use the same fields, check against the CPU.
	metropolisUpdateLR(p, iter);
      }      
      if(!p.useGPUCluster) GPU_copyArraysToHost(p);
#else
      cout<<"GPU Metro not built"<<endl;
#endif
    } else {
      metropolisUpdateLR(p, iter);
    }
  }
}


void PhiFourth2D::measure(observables &obs, int &idx, Param p) {
  
  double rhoVol  = 1.0/p.surfaceVol;
  double tmpPhi3 = 0.0;
  
  if(p.coupling_type == SR) {
    obs.tmpE = actionSR(phi, s, p, obs.KE, obs.PE);
  } else {
    obs.tmpE = actionLR(phi, s, p, LR_couplings, obs.KE, obs.PE);
  }

  obs.aveKE += obs.KE;
  obs.avePE += obs.PE;
  obs.aveE  += obs.tmpE;
  obs.aveE2 += obs.tmpE*obs.tmpE;
  
  obs.MagPhi = 0.0;
  for(int i = 0;i < p.surfaceVol; i++) {
    obs.MagPhi += phi[i];
    tmpPhi3 += phi[i]*phi[i]*phi[i];
  }
  obs.MagPhi  *= rhoVol;
  tmpPhi3     *= rhoVol;
  
  obs.avePhiAb += abs(obs.MagPhi);
  obs.avePhi   += obs.MagPhi;
  obs.avePhi2  += obs.MagPhi*obs.MagPhi;
  obs.phi3Ave  += tmpPhi3;  
  obs.avePhi4  += obs.MagPhi*obs.MagPhi*obs.MagPhi*obs.MagPhi;
  
  obs.E_arr[idx]     = obs.tmpE;
  obs.E2_arr[idx]    = obs.tmpE*obs.tmpE;
  obs.PhiAb_arr[idx] = abs(obs.MagPhi);
  obs.Phi_arr[idx]   = obs.MagPhi;
  obs.Phi2_arr[idx]  = obs.MagPhi*obs.MagPhi;
  obs.Phi3_arr[idx]  = obs.MagPhi*obs.MagPhi;
  obs.Phi4_arr[idx]  = obs.MagPhi*obs.MagPhi*obs.MagPhi*obs.MagPhi;
  
  idx++;
  
  cout<<setprecision(8);
  double norm = 1.0/(idx);

  obs.Suscep[idx-1]   = (obs.avePhi2*norm - pow(obs.avePhiAb*norm,2))/rhoVol;
  obs.SpecHeat[idx-1] = (obs.aveE2*norm - pow(obs.aveE*norm, 2));
  obs.Binder[idx-1]   = 1.0-obs.avePhi4/(3.0*obs.avePhi2*obs.avePhi2*norm);
  
  //Dump to stdout
  cout<<"Measurement "<<idx<<endl;
  cout<<"Ave Energy= "<<obs.aveE*norm<<endl;
  cout<<"Ave KE    = "<<obs.aveKE*norm<<endl;
  cout<<"Ave PE    = "<<obs.avePE*norm<<endl;
  cout<<"Ave |phi| = "<<obs.avePhiAb*norm<<endl;
  cout<<"Ave phi   = "<<obs.avePhi*norm<<endl;
  cout<<"Ave phi^2 = "<<obs.avePhi2*norm<<endl;
  cout<<"Ave phi^3 = "<<obs.phi3Ave*norm<<endl;
  cout<<"Ave phi^4 = "<<obs.avePhi4*norm<<endl;
  cout<<"Suscep    = "<<obs.Suscep[idx-1]<<endl;
  cout<<"Spec Heat = "<<obs.SpecHeat[idx-1]<<endl;
  cout<<"Binder    = "<<obs.Binder[idx-1]<<endl;
  
}
