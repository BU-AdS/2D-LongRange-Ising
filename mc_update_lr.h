#ifndef MC_2D_ISING_LR_H
#define MC_2D_ISING_LR_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "mc_util.h"
#include "util.h"

//Long range cluster and Metropolis routines.
void wolffUpdateLR(double *phi_arr, int *s, Param p,
		   double *LR_couplings,
		   double &delta_mag_phi, int iter);

void wolffClusterAddLR(int i, int *s, int cSpin, 
		       std::vector<int> &Rcluster, 
		       double *LR_couplings, double *phi, Param p);

void swendsenWangUpdateLR(double *phi_arr, int *s, Param p,
			  double *LR_coupling,
			  double &delta_mag_phi, int iter);

void swendsenWangClusterAddLR(int i, int *s, int cSpin, int clusterNum, 
			      int *clusterDef, std::vector<int> Rcluster, 
			      double *LR_couplings,
			      double *phi_arr, Param p);

void clusterPossibleLR(int i, int *s, int cSpin,
		       bool *Pcluster, std::vector<int> &Rcluster, 
		       Param p);


int metropolisUpdateLR(double *phi_arr, int *s, Param &p,
		       double *LR_couplings_sq,
		       double & delta_mag_phi, int iter);


#endif
