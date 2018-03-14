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
  double delta_mag_phi = 0.0;
  
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

  //Array for the LR kinetic denominators
  double *denom;
  
  //Arrays to hold spin and phi values.  
  int *s;
  double *phi;
  
  //Running correlation function arrays.
  double *run_corr_t;
  double *run_corr_s;

  //Individual correlation function arrays..
  double **ind_corr_t;
  double **ind_corr_s;

  //The correlation functions must be weighted by how frequently
  //the i,j separation occurs  
  int **corr_norms;

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

  //Array of random numbers on the GPU
  double* gpu_rands;
  
  double *gpu_LR_couplings;
  double *gpu_denom;
  double *gpu_phi;
  int *gpu_s;
  
  void GPU_wolffClusterAddLR(int i, Param p, int cSpin, double *gpu_rands, int *added);
  
  void GPU_wolffUpdateLR(Param p, int iter);
    
#endif  
  
  //Constructor
  MonteCarlo2DIsing(Param p);

  //Do n_cluster updates, and one metro update.
  void updateIter(observables obs, Param p, int iter);

  //Do some initial Metro updates, then run updateIter
  //until thermalisation.
  void thermalise(double *phi, int *s, Param p, observables obs);

  //Self explanatory...
  void runSimulation(Param p);
  void createDenom(Param p);
  void createLRcouplings(Param p);
  
  //--------------------------------------------------------//
  // Wrappers for the different LR and SR cluster functions //
  //--------------------------------------------------------//
  void wolffUpdate(double *phi, int *s, Param p, 
		   double &delta_mag_phi, int iter);
  
  void swendsenWangUpdate(double *phi, int *s, Param p,
			  double &delta_mag_phi, int iter);
  
  void metropolisUpdate(double *phi, int *s, Param p,
			double &delta_mag_phi, int iter);
  
};

//Extract measurements from simulation
void measure(observables &obs, double *phi, int *s, 
	     double *LR_c, int &idx, Param p);

//Calcuate the Jack-Knife error and write the correlation functions
//to file
void writeObservables(double **ind_corr_t, double *run_corr_t, 
		      double **ind_corr_s, double *run_corr_s,
		      int **corr_norms, int idx, observables obs, Param p);

//Square 2D Ising short range utils
double actionSR(double *phi_arr, int *s, 
		Param p, double &KE, double &PE);

//Square 2D Ising long range utils
double actionLR(double *phi_arr, int *s, Param p,
		double *LR_couplings,
		double &KE, double &PE);

//Coupling creation function
void LRAdSCouplings(double *LR_couplings, Param &p);


#endif
