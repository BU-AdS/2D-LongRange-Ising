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
#include "gpuMcPhiFourth2D.cuh"

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

  //Running phi^3 correlation function arrays.
  run_corr_phi3 = (double*)malloc(arr_len*sizeof(double));
  //Individual phi^3 correlation functions.
  ind_corr_phi3 = (double**)malloc(p.n_meas*sizeof(double*));
  
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
      correlatorsImpWolff(ind_corr, run_corr,
			  meas-1, obs.avePhi*norm, p);
      
      auto elapsed1 = std::chrono::high_resolution_clock::now() - start1;
      time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
      cout<<"Correlation Function Calculation time = "<<time/(1.0e6)<<endl;
      
      
      /*
	bool isTemporal;
	isTemporal = true;
	correlators(ind_corr_t, run_corr_t, isTemporal, meas-1, 
	obs.avePhi*norm, p);
	isTemporal = false;
	correlators(ind_corr_s, run_corr_s, isTemporal, meas-1,
	obs.avePhi*norm, p);
      */
      
      //Jacknife and dump the data
      if(meas%10 == 0) {	
	writeObservables(ind_corr, run_corr, norm_corr, meas, obs, p);
      }
      
      //Calculate the autocorrelation of |phi|
      autocorrelation(obs.PhiAb_arr, obs.avePhiAb*norm, meas, auto_corr);
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
void PhiFourth2D::wolffUpdate(Param p, int iter) {

  for(int i=0; i<p.n_cluster; i++) {
    if(p.coupling_type == SR) {
      wolffUpdateSR(phi, s, p, iter);
    }
    else {
      if(p.useGPUCluster) {
	if(!p.useGPUMetro) {
	  //We must copy from the host to the device
#ifdef USE_GPU
	  GPU_copyArraysToDevice(p);
#endif
	}
#ifdef USE_GPU	
	int rand_site = int(unif(rng) * p.surfaceVol);	
	GPU_wolffUpdateLR(p, rand_site, iter, i);
#endif
	if(!p.useGPUMetro) {
	  //We must copy from the device to the host
#ifdef USE_GPU
	  GPU_copyArraysToHost(p);
#endif
	}	
      } else {
	wolffUpdateLR(phi, s, p, LR_couplings, iter, i);
      }
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
      if(!p.useGPUCluster) {
	//We must copy from the host to the device
#ifdef USE_GPU
	GPU_copyArraysToDevice(p);
#endif
      }
#ifdef USE_GPU
      GPU_metropolisUpdateLR(p, iter);
#endif
      if(p.doMetroCheck) {
	//Use the same fields, check against the CPU.
	metropolisUpdateLR(p, iter);
      }
      
      if(!p.useGPUCluster) {
	//We must copy from the device to the host
#ifdef USE_GPU
	GPU_copyArraysToHost(p);
#endif
      }
    } else {
      metropolisUpdateLR(p, iter);
    }
  }
}


void PhiFourth2D::measure(observables &obs, int &idx, Param p) {
  
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
  }
  obs.MagPhi *= rhoVol;
  
  obs.avePhiAb += abs(obs.MagPhi);
  obs.avePhi   += obs.MagPhi;
  obs.avePhi2  += obs.MagPhi*obs.MagPhi;
  obs.avePhi4  += obs.MagPhi*obs.MagPhi*obs.MagPhi*obs.MagPhi;
  
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
