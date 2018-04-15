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
#include "mcIsing2D.h"

#ifdef USE_GPU
//#include "gpuIsing2D.cuh"
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

  //Long Range coupling array
  isingProb = (double*)malloc(arr_len*sizeof(double));
  
  //Arrays to hold spin and phi values.  
  s = (int*)malloc(vol*sizeof(int));
  s_cpy = (int*)malloc(vol*sizeof(int));

  //Running correlation function arrays.
  run_corr = (double*)malloc(arr_len*sizeof(double));
  //Correlation function normalisation array.
  norm_corr = (int*)malloc(arr_len*sizeof(int));  
  //Individual correlation functions.
  ind_corr = (double**)malloc(p.n_meas*sizeof(double*));
  
  //Autocorrelation
  auto_corr = (double*)malloc(p.n_meas*sizeof(double));
  
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
  for(int i=0; i<arr_len; i++) {
    run_corr[i] = 0.0;
    norm_corr[i] = 0;
  }

  int t1,x1,t2,x2,dt,dx,idx;
  for(int i=0; i<p.surfaceVol; i++) {
    t1 = i/S1;
    x1 = i%S1;
    for(int j=0; j<p.surfaceVol; j++) {
      //Index divided by circumference, using the int floor feature/bug,
      //gives the timeslice index.
      t2 = j / S1;
      dt = abs(t2-t1) > p.Lt/2 ? p.Lt - abs(t2-t1) : abs(t2-t1);
      
      //The index modulo the circumference gives the spatial index.
      x2 = j % S1;            
      dx = abs(x2-x1) > p.S1/2 ? p.S1 - abs(x2-x1) : abs(x2-x1);      
      
      idx = dx + x_len*dt;
      norm_corr[idx]++;//FIXME
    }
  }

  //Individual correlation function arrays, autocorrelation.
  for(int i=0; i<p.n_meas; i++) {
    auto_corr[i] = 0.0;    
    ind_corr[i] = (double*)malloc(arr_len*sizeof(double));
    for(int j=0; j<arr_len; j++) ind_corr[i][j] = 0.0;
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

  cudaMalloc((void**) &gpu_isingProb, arr_len*sizeof(double));  
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
    //Create LR couplings
    long long time = 0.0;
    auto start1 = std::chrono::high_resolution_clock::now();
    createLRcouplings(p);
    exponentiateLRcouplings(p);
    auto elapsed1 = std::chrono::high_resolution_clock::now() - start1;
    time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
    cout<<"Coupling creation time = "<<time/(1.0e6)<<endl;
  }
  
  thermalise(p);
  
  //reset timing variables.
  cluster = 0.0;
  for(int iter = p.n_therm; iter < p.n_therm + p.n_skip*p.n_meas; iter++) {
    
    updateIter(p, iter);
    
    //Take measurements.
    if((iter+1) % p.n_skip == 0) {

      //Timing
      int therm_iter = iter - p.n_therm;
      cout<<"(Thermalised) Average time per cluster update = "<<cluster/(therm_iter*p.n_cluster*1.0e6)<<"s"<<endl;

      //If using the GPU, copy the spin array to the host.
      if(p.useGPUCluster) GPU_copyArraysToHostI(p);
      
      //update running average, dump to stdout. 
      //'meas' is ++ incremented in this function.
      measureI(obs, meas, p);
      
      norm = 1.0/(meas);

      //Visualisation tool
      visualiserIsing(s, p);
      
      //Calculate correlaton functions and update the average.
      //int rand_site = int(unif(rng) * p.surfaceVol);	      
      //GPU_correlatorsImpWolff(rand_site, meas-1, obs.avePhi*norm, p);
  
      long long time = 0.0;
      auto start1 = std::chrono::high_resolution_clock::now();
      correlatorsImpWolffI(ind_corr, run_corr,
			   meas-1, obs.avePhi*norm, p);
	
      auto elapsed1 = std::chrono::high_resolution_clock::now() - start1;
      time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed1).count();
      cout<<"Correlation Function Calculation time = "<<time/(1.0e6)<<endl;
      
      
      //Jacknife and dump the data
      if(meas%10 == 0) {	
	writeObservables(ind_corr, run_corr, norm_corr, meas, obs, p);
      }
      
      //Calculate the autocorrelation of |phi|
      autocorrelation(obs.PhiAb_arr, obs.avePhiAb*norm, meas, auto_corr);

    }
  }
}

void Ising2D::updateIter(Param p, int iter) {
  
  auto start = std::chrono::high_resolution_clock::now();    
  if(p.useWolff)  wolffUpdate(p, iter);
  else swendsenWangUpdate(p, iter);  
  auto elapsed = std::chrono::high_resolution_clock::now() - start;
  cluster += std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();  

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
  
  //Reset timing variables
  cluster = 0.0;
  
  //Cluster steps
  for(int iter = 0; iter < p.n_therm; iter++) {
    updateIter(p, iter);
    if( (iter+1)%p.n_skip == 0 ) {
      cout<<"Average time per cluster update = "<<cluster/((iter+1)*1.0e6)<<"s"<<endl; 
    }      
  }
  cout<<endl<<"Thermalisaton complete."<<endl<<endl;  
}

//------------------------------------------------//
// Wrappers for the different LR and SR functions //
//------------------------------------------------//
void Ising2D::wolffUpdate(Param p, int iter) {

  if(p.coupling_type == SR) {
    wolffUpdateISR(s, p, iter);
  }
  else {
    if(p.useGPUCluster) {
#ifdef USE_GPU
      int rand_site = int(unif(rng) * p.surfaceVol);
      GPU_wolffUpdateILR(p, rand_site, iter);
#else
      cout<<"ERROR: GPU support not built. Please revise your build options"<<endl;
      cout<<"[Generated by Ising2D::wolffUpdate(Param p, int iter)]"<<endl;
      exit(0);
#endif
    } else {
      wolffUpdateILR(s, p, isingProb, iter);
    }
  }
}


void Ising2D::swendsenWangUpdate(Param p, int iter){
  
  if(p.coupling_type == SR) 
    swendsenWangUpdateISR(s, p, iter); 
  else
    swendsenWangUpdateILR(p, iter); 
}

void Ising2D::metropolisUpdate(Param p, int iter) {
  
  if(p.coupling_type == SR) {
    metropolisUpdateISR(s, p, iter);
  } else {
    metropolisUpdateILR(p, iter);
  }
}



inline double coupling(double dt, double dth, Param p) {
  
  dt  *= M_PI/p.S1;
  dth *= M_PI/p.S1;  
  return pow( (cosh(dt) - cos(dth)) , -(1+p.sigma/2));
  
}

void Ising2D::createLRcouplings(Param &p) {
  
  double sigma = p.sigma;
  int S1 = p.S1;
  int Lt = p.Lt;
  int x_len = S1/2 + 1;
  int t_len = Lt/2 + 1;
  
  double couplingNormTheta = 1.0/coupling(0, 1, p);
  //double couplingNormTime = 1.0/coupling(1, 0, p);

  double sum_1 = 0.0;
  double sum_2 = 0.0;
  
  LR_couplings[0] = couplingNormTheta;  

  //Scales the temporal direction so that the coupling it unit valued
  //at one latice spacing.
  p.t_scale = acosh(1 + pow(couplingNormTheta,1.0/(1+p.sigma/2)))*S1/M_PI;

  for(int j=0; j<t_len; j++){
    for(int i=0; i<x_len; i++){

      //Associate dx with i, dt with j. dx runs fastest
      int idx = i + j*x_len;      
      if(i!=0 || j!=0) {
	if(p.usePowLaw) LR_couplings[idx] = pow(sqrt(i*i + j*j),-(2+sigma));
	else  LR_couplings[idx] = couplingNormTheta*coupling((double)j*p.t_scale, (double)i, p);
	sum_1 += LR_couplings[idx];
	sum_2 += LR_couplings[idx]*LR_couplings[idx];
      }
      printf("%.8e ", LR_couplings[idx]);
    }
    cout<<endl;
  }

  cout<<endl<<endl<<"Sum 1 = "<<sum_1<<" sum 2 = "<<sum_2<<endl<<endl;
  
#ifdef USE_GPU
  cudaMemcpy(gpu_LR_couplings, LR_couplings,
	     x_len*t_len*sizeof(double), cudaMemcpyHostToDevice);
#endif
  
}

//In the Ising cluster algorithm, the probability of flipping a site is
//wholly dependent on spatial separation and the J value. We can therefore
//pre-compute the value of 1-exp(-2*J*LRcoupling[idx]) used in Wolff.
void Ising2D::exponentiateLRcouplings(Param p) {
  
  int S1 = p.S1;
  int Lt = p.Lt;
  double J = p.J;
  int x_len = S1/2 + 1;
  int t_len = Lt/2 + 1;
    
  for(int j=0; j<t_len; j++){
    for(int i=0; i<x_len; i++){

      //Associate dx with i, dt with j. dx runs fastest
      int idx = i + j*x_len;      
      isingProb[idx] = 1 - exp(-2*J*LR_couplings[idx]);	
      
    }
  }
#ifdef USE_GPU
  cudaMemcpy(gpu_isingProb, isingProb,
	     x_len*t_len*sizeof(double), cudaMemcpyHostToDevice);
#endif

}

//NB we use the same variable names for both Ising and Phi4th,
//e.g. the phi4th field value phi is used for the ising spin value sigma.
void Ising2D::measureI(observables &obs, int &idx, Param p) {

  int vol = p.surfaceVol;
  double rhoVol  = 1.0/vol;
  double J=p.J;
  
  //Energy
  if(p.coupling_type == SR) {
    obs.tmpE = energyISR(s, p, obs.KE);
  } else {
    obs.tmpE = energyILR(s, p, LR_couplings, obs.KE);
  }

  obs.aveE  += rhoVol*obs.tmpE;
  obs.aveE2 += rhoVol*obs.tmpE*rhoVol*obs.tmpE;
  
  obs.MagPhi = 0.0;
  for(int i = 0;i < vol; i++) {
    obs.MagPhi += s[i];
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

  obs.Suscep[idx-1]   = (obs.avePhi2*norm - pow(obs.avePhiAb*norm, 2))*vol/J;
  obs.SpecHeat[idx-1] = (obs.aveE2*norm - pow(obs.aveE*norm, 2))*vol/(J*J);
  obs.Binder[idx-1] = 1.0-obs.avePhi4/(3.0*obs.avePhi2*obs.avePhi2*norm);
  
  //Dump to stdout
  cout<<"Measurement "<<idx<<endl;
  cout<<"Ave Energy= "<<obs.aveE*norm<<endl;
  cout<<"Ave |s|   = "<<obs.avePhiAb*norm<<endl;
  cout<<"Ave s     = "<<obs.avePhi*norm<<endl;
  cout<<"Ave s^2   = "<<obs.avePhi2*norm<<endl;
  cout<<"Ave s^4   = "<<obs.avePhi4*norm<<endl;
  cout<<"Suscep    = "<<obs.Suscep[idx-1]<<endl;
  cout<<"Spec Heat = "<<obs.SpecHeat[idx-1]<<endl;
  cout<<"Binder    = "<<obs.Binder[idx-1]<<endl;

}
