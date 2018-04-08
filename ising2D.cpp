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

#ifdef USE_GPU
#include "gpuMcPhiFourth2D.cuh"
#endif

using namespace std;

extern int seed;
extern mt19937 rng;
extern uniform_real_distribution<double> unif;
extern int CACHE_LINE_SIZE;
extern int sze;

//----------------//
// 2D Ising class //
//----------------//
Ising2D::Ising2D(Param p) {

  long long time = 0.0;
  auto start1 = std::chrono::high_resolution_clock::now();
  cout<<"CPU malloc and init start..."<<endl;

  int S1 = p.S1;
  int Lt = p.Lt;
  int vol = p.surfaceVol;
  int x_len = S1/2 + 1;
  int t_len = Lt/2 + 1;
  int arr_len = x_len*t_len;
  
  //Long Range coupling array
  LR_couplings = (double*)malloc(arr_len*sizeof(double));
  //Arrays to hold spin and phi values.  
  s = (int*)malloc(vol*sizeof(int));
  s_cpy = (int*)malloc(vol*sizeof(int));

  //Running correlation function arrays.
  run_corr_t = (double*)malloc((Lt/2+1)*sizeof(double));
  run_corr_s = (double*)malloc((S1/2+1)*sizeof(double));
  //Individual correlation functions.
  ind_corr_t = (double**)malloc((p.n_meas)*sizeof(double*));
  ind_corr_s = (double**)malloc((p.n_meas)*sizeof(double*));

  //Autocorrelation
  auto_corr = (double*)malloc((p.n_meas)*sizeof(double));
  
  //Initialise Phi and spin fields.
  for(int i = 0;i < vol; i++) {
    s[i] = (2.0*unif(rng) - 1.0 > 0) ? 1 : -1;
    //cout<<"s["<<i<<"] = "<<s[i]<<endl;
  } 

  //Long-range coupling
  for(int j=0; j<t_len; j++)
    for(int i=0; i<x_len; i++)
      LR_couplings[i+j*x_len] = 0.0;

  //Running correlation function arrays
  for(int i=0; i<Lt/2+1; i++) run_corr_t[i] = 0.0;
  for(int i=0; i<S1/2+1; i++) run_corr_s[i] = 0.0;

  //Individual correlation function arrays, autocorrelation.
  for(int i=0; i<p.n_meas; i++) {
    auto_corr[i] = 0.0;    
    ind_corr_t[i] = (double*)malloc((Lt/2+1)*sizeof(double));
    ind_corr_s[i] = (double*)malloc((S1/2+1)*sizeof(double));    
    for(int j=0; j<Lt/2+1; j++) ind_corr_t[i][j] = 0.0;
    for(int j=0; j<S1/2+1; j++) ind_corr_s[i][j] = 0.0;
  }
  
  auto elapsed1 = std::chrono::high_resolution_clock::now() - start1;
  time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
  cpu_added = (bool*)malloc(vol*sizeof(bool));
  cout<<"CPU malloc and init time = "<<time/(1.0e6)<<endl;
  
#ifdef USE_GPU

  //Device memory allocations.
  cout<<"CUDA malloc start..."<<endl;
  
  start1 = std::chrono::high_resolution_clock::now();
  //Allocate space on the GPU.
  cudaMalloc((void**) &states, vol*sizeof(curandState_t));  
  
  cudaMalloc((void**) &gpu_LR_couplings, arr_len*sizeof(double));  
  cudaMalloc((void**) &gpu_rands, vol*sizeof(double));
  cudaMalloc((void**) &gpu_rands_aux, vol*sizeof(double));  
  cudaMalloc((void**) &gpu_result, vol*sizeof(double));
  cudaMalloc((void**) &gpu_ind_corr_s, x_len*sizeof(double));
  cudaMalloc((void**) &gpu_ind_corr_t, t_len*sizeof(double));
  
  cudaMalloc((void**) &gpu_s, vol*sizeof(int));
  cudaMalloc((void**) &gpu_s_cpy, vol*sizeof(int));

  cudaMalloc((void**) &gpu_added, vol*sizeof(bool));
  cudaMalloc((void**) &gpu_cluster, vol*sizeof(bool));  
  
  elapsed1 = std::chrono::high_resolution_clock::now() - start1;
  time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
  cout<<"CUDA malloc time = "<<time/(1.0e6)<<endl;

  start1 = std::chrono::high_resolution_clock::now();
  GPU_initRand(p, seed, states);
  elapsed1 = std::chrono::high_resolution_clock::now() - start1;
  time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
  cout<<"CUDA rand states init time = "<<time/(1.0e6)<<endl;
  
  debug_arr1 = (double*)malloc(vol*sizeof(double));
  debug_arr2 = (double*)malloc(vol*sizeof(double));
  for(int j=0; j<S1; j++){
    for(int i=0; i<Lt; i++){
      debug_arr1[i+j*Lt] = 0.0;
      debug_arr2[i+j*Lt] = 0.0;
    }
  }
#endif  
}

void Ising2D::runSimulation(Param p) {

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
  
  //thermalise(p);
    
  //reset timing variables.
  metro = 0.0;
  cluster = 0.0;
  for(int iter = p.n_therm + p.n_metro_cool; iter < p.n_therm + p.n_skip*p.n_meas + p.n_metro_cool; iter++) {
    
    //updateIter(p, iter);
    
    //Take measurements.
    if((iter+1) % p.n_skip == 0) {

      //Timing
      int therm_iter = iter - (p.n_therm + p.n_metro_cool);
      cout<<"(Thermalised) Average time per Metro update   = "<<metro/(therm_iter*1.0e6)<<"s"<<endl;
      cout<<"(Thermalised) Average time per cluster update = "<<cluster/(therm_iter*p.n_cluster*1.0e6)<<"s"<<endl;

      //copy GPU arrays to the CPU FIXME: do all measurements on the GPU      
      //GPU_copyArraysToHost(p);      
      
      //update running average, dump to stdout. 
      //'meas' is ++ incremented in this function.
      //measure(obs, meas, p);
      
      norm = 1.0/(meas);

      //Visualisation tool
      //visualiserIsing(s, obs.avePhiAb*norm, p);
      
      //Calculate correlaton functions and update the average.
      //int rand_site = int(unif(rng) * p.surfaceVol);	      
      //GPU_correlatorsImpWolff(rand_site, meas-1, obs.avePhi*norm, p);

      /*
      long long time = 0.0;
      auto start1 = std::chrono::high_resolution_clock::now();
      correlatorsImpWolff(ind_corr_t, run_corr_t,
			  ind_corr_s, run_corr_s,
			  meas-1, obs.avePhi*norm, p);
      
      auto elapsed1 = std::chrono::high_resolution_clock::now() - start1;
      time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
      cout<<"Correlation Function Calculation time = "<<time/(1.0e6)<<endl;
      
      
      // bool isTemporal;
      // isTemporal = true;
      // correlators(ind_corr_t, run_corr_t, isTemporal, meas-1, 
      // obs.avePhi*norm, p);
      // isTemporal = false;
      // correlators(ind_corr_s, run_corr_s, isTemporal, meas-1,
      // obs.avePhi*norm, p);
      
      //Jacknife and dump the data
      if(meas%10 == 0) {	
	writeObservables(ind_corr_t, run_corr_t, ind_corr_s,
			 run_corr_s, meas, obs, p);
      }
      
      //Calculate the autocorrelation of |phi|
      autocorrelation(obs.PhiAb_arr, obs.avePhiAb*norm, meas, auto_corr);
      */
    }
  }
}
/*
void Ising2D::updateIter(Param p, int iter) {
  
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

void Ising2D::thermalise(Param p) {
  
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
void Ising2D::wolffUpdate(Param p, int iter) {

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

void Ising2D::swendsenWangUpdate(Param p, int iter){
  
  if(p.coupling_type == SR) 
    swendsenWangUpdateSR(phi, s, p, iter); 
  else
    swendsenWangUpdateLR(p, iter); 
}

void Ising2D::metropolisUpdate(Param p, int iter) {
  
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
*/

inline double Coupling(double dt, double dth, Param p) {

  
  dt  *= M_PI/p.S1;
  dth *= M_PI/p.S1;  

  return 1.0/pow( (cosh(dt) - cos(dth)) , 1+p.sigma/2);

}

void Ising2D::createLRcouplings(Param p) {
  
  double sigma = p.sigma;
  int S1 = p.S1;
  int Lt = p.Lt;
  int x_len = S1/2 + 1;
  int t_len = Lt/2 + 1;
  int arr_len = x_len*t_len;

  double couplingNormTheta = 1.0/Coupling(0, 1, p);
  //double couplingNormTime = 1.0/Coupling(1, 0, p);

  double sum_1 = 0.0;
  double sum_2 = 0.0;
  
  LR_couplings[0] = couplingNormTheta;  
  
  for(int j=0; j<t_len; j++){
    for(int i=0; i<x_len; i++){

      //Associate dx with i, dt with j. dx runs fastest
      int idx = i + j*x_len;      
      if(i!=0 || j!=0) {
	if(p.usePowLaw) LR_couplings[idx] = pow(sqrt(i*i + j*j),-(2+sigma));
	else  LR_couplings[idx] = couplingNormTheta*Coupling((double)j*p.t_scale, (double)i, p);
	sum_1 += LR_couplings[idx];
	sum_2 += LR_couplings[idx]*LR_couplings[idx];
      }
      printf("%.3e ", LR_couplings[idx]);
    }
    cout<<endl;
  }

  cout<<endl<<endl<<"Sum 1 = "<<sum_1<<" sum 2 = "<<sum_2<<endl<<endl;
  
#ifdef USE_GPU
  cudaMemcpy(gpu_LR_couplings, LR_couplings,
	     arr_len*sizeof(double), cudaMemcpyHostToDevice);
#endif
  
}

