#ifndef PHIFOURTH2D_H
#define PHIFOURTH2D_H

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

#include "ising2D.h"

class PhiFourth2D: public Ising2D {
  
 public: 
  double *phi;
  double *phi_cpy;

  double **ind_corr_phi_phi3;
  double *run_corr_phi_phi3;
  
  double **ind_corr_phi3_phi3;
  double *run_corr_phi3_phi3;
  
  //Constructor
  PhiFourth2D(Param p);
  
  //Do n_cluster updates, and one metro update.
  void updateIter(Param p, int iter);
  
  //Do some initial Metro updates, then run updateIter
  //until thermalisation.
  void thermalise(Param p);
  
  void runSimulation(Param p);

  //Extract measurements from simulation
  void measure(observables &obs, int &idx, Param p);
  
  void swendsenWangUpdateLR(Param p, int iter);
  //Wrappers for the different LR and SR cluster functions.
  void wolffUpdate(Param p, int iter);  
  void swendsenWangUpdate(Param p, int iter);
  void metropolisUpdate(Param p, int iter);
  void metropolisUpdateLR(Param p, int iter);
    
#ifdef USE_GPU
  
  double *gpu_phi;
  double *gpu_phi_cpy;
  double *gpu_ind_corr_s;

  double GPU_metropolisUpdateLR(Param p, int iter);
  void GPU_wolffUpdateLR(Param p, int rand_site, int iter, int i);
  void GPU_copyArraysToHost(Param p);
  void GPU_copyArraysToDevice(Param p);
  void GPU_correlatorsImpWolff(int i, int meas, double avePhi, Param p);

#endif
  
  //Correlation function calculator
  void correlators(int meas, double avePhi, Param p);
  void correlatorsImpSW(int meas, double avePhi, Param p); //FIXME
  void correlatorsImpWolff(int meas, double avePhi, Param p);  
  void correlatorsPhi3(int meas, double avePhi, double avePhi3, Param p);
  
  
};


#endif
