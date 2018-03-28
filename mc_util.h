#ifndef MC_UTIL_H
#define MC_UTIL_H

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

class observables {
  
 private:
  int n_meas = 1;
  
 public:
  //Arrays holding measurements for error analysis
  double *E_arr;
  double *E2_arr;
  double *PhiAb_arr;
  double *Phi_arr;
  double *Phi2_arr;
  double *Phi4_arr;
  double *Suscep;
  double *SpecHeat;
  double *Binder;
  
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

  //temps
  double KE = 0.0;
  double PE = 0.0;
  
  //constructor
  observables(int meas);
  
  //destructor
  //~observable();
  
};


//---------------------------------------//
//Square 2D Ising Monte Carlo base class //
//---------------------------------------//

class MonteCarlo2DIsing {
  
 public:
  
  //Array for the LR couplings.
  double *LR_couplings;

  //Arrays to hold spin and phi values.  
  int *s;
  double *phi;
  
  //Running correlation function arrays.
  double *run_corr_t;
  double *run_corr_s;

  //Individual correlation function arrays..
  double **ind_corr_t;
  double **ind_corr_s;

  //The current measurement number.
  int meas = 0;
  
  //1.0/meas
  double norm = 0.0;
  
  //Object to hold MC update times
  long long metro = 0.0;
  long long cluster = 0.0;

#ifdef USE_GPU
  //GPU objects
  // CUDA's random number library uses curandState_t to keep track of the seed value
  // we will store a random state for every thread  
  curandState_t* states;
  curandState_t* states_aux;

  double *gpu_rands;
  double *gpu_rands_aux;
  double *gpu_result;
  double *gpu_LR_couplings;
  double *gpu_phi;
  
  bool *cpu_added;
  bool *gpu_added;

  int *gpu_s;

  double *debug_arr1;
  double *debug_arr2;
  
  void GPU_initRand(Param p, int seed, curandState_t *states);

  double GPU_metropolisUpdateLR(Param p, int iter);
  void GPU_wolffUpdateLR(Param p, int rand_site, int iter, int i);
  void GPU_copyArraysToHost(Param p);
  void GPU_copyArraysToDevice(Param p);
  
#endif  
  
  //Constructor
  MonteCarlo2DIsing(Param p);

  //Do n_cluster updates, and one metro update.
  void updateIter(Param p, int iter);

  //Do some initial Metro updates, then run updateIter
  //until thermalisation.
  void thermalise(Param p);

  //Self explanatory...
  void runSimulation(Param p);
  void createLRcouplings(Param p);
  
  //--------------------------------------------------------//
  // Wrappers for the different LR and SR cluster functions //
  //--------------------------------------------------------//
  void wolffUpdate(Param p, int iter);
  
  void swendsenWangUpdate(Param p, int iter);
  
  void metropolisUpdate(Param p, int iter);

  int metropolisUpdateLR(Param &p, int iter);
    
};

//Extract measurements from simulation
void measure(observables &obs, double *phi, int *s, 
	     double *LR_c, int &idx, Param p);

//Calcuate the Jack-Knife error and write the correlation functions
//to file
void writeObservables(double **ind_corr_t, double *run_corr_t, 
		      double **ind_corr_s, double *run_corr_s,
		      int idx, observables obs, Param p);

//Square 2D Ising short range utils
double actionSR(double *phi_arr, int *s, 
		Param p, double &KE, double &PE);

//Square 2D Ising long range utils
double actionLR(double *phi_arr, int *s, Param p,
		double *LR_couplings,
		double &KE, double &PE);


#endif
