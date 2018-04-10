#ifndef ISING2D_H
#define ISING2D_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#ifdef USE_GPU
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#endif


#include "util.h"

//----------------//
// 2D Ising class //
//----------------//

class Ising2D {
  
 public:
  
  //Array for the LR couplings.
  double *LR_couplings;

  //Array to hold spin values.  
  int *s;
  
  //Running correlation function arrays.
  double *run_corr_t;
  double *run_corr_s;
  double **run_corr;
  
  //Individual correlation function arrays..
  double **ind_corr_t;
  double **ind_corr_s;
  double ***ind_corr;
  
  //Holds Autocorrelation data
  double *auto_corr;
  
  //1.0/meas
  double norm = 0.0;
  
  //Objects to hold MC update times
  long long metro = 0.0;
  long long cluster = 0.0;

  int *s_cpy;
  bool *cpu_added;

  double *debug_arr1;
  double *debug_arr2;

  #ifdef USE_GPU  
  // CUDA's random number library uses curandState_t to keep track of the seed value
  // we will store a random state for every thread  
  curandState_t* states;
  curandState_t* states_aux;
  void GPU_initRand(Param p, int seed, curandState_t *states);
#endif

  double *gpu_rands;
  double *gpu_rands_aux;
  double *gpu_result;
  double *gpu_LR_couplings;

  double *gpu_ind_corr_t;//FIXME
  double *gpu_ind_corr_s;//FIXME 
  
  bool *gpu_added;
  bool *gpu_cluster;
  
  int *gpu_s;
  int *gpu_s_cpy;

  double GPU_metropolisUpdateLR(Param p, int iter);
  void GPU_wolffUpdateLR(Param p, int rand_site, int iter, int i);
  void GPU_copyArraysToHost(Param p);
  void GPU_copyArraysToDevice(Param p);
  void GPU_correlatorsImpWolff(int i, int meas, double avePhi, Param p);
  
  //Constructor
  Ising2D(Param p);

  //Do n_cluster updates, and one metro update.
  void updateIter(Param p, int iter);

  //Do some initial Metro updates, then run updateIter
  //until thermalisation.
  void thermalise(Param p);
  
  void runSimulation(Param p);
  
  void measureI(observables &obs, int &idx, Param p);
  
  void createLRcouplings(Param p);

  void swendsenWangUpdateILR(Param p, int iter);  
  //Wrappers for the different LR and SR cluster functions.
  void wolffUpdate(Param p, int iter);  
  void swendsenWangUpdate(Param p, int iter);
  void metropolisUpdate(Param p, int iter);
  void metropolisUpdateILR(Param p, int iter);
  
  void correlatorsImpWolffI(double **ind_corr_s, double *run_corr_s,
			    double **ind_corr_t, double *run_corr_t,
			    int meas, double aveS, Param p);

};

#endif
