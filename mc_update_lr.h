#ifndef MC_UPDATE_LR_H
#define MC_UPDATE_LR_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "mc_util.h"
#include "util.h"

//Long range cluster and Metropolis routines.
void wolffUpdateLR(double *phi, int *s, Param p,
		   double *LR_couplings, int iter, int n);

void wolffClusterAddLR(int i, int *s, int cSpin, 
		       double *LR_couplings, double *phi, Param p);

void swendsenWangUpdateLR(double *phi, int *s, Param p,
			  double *LR_couplings, int iter);

void swendsenWangClusterAddLR(int i, int *s, int cSpin, int clusterNum, 
			      int *clusterDef, std::vector<int> Rcluster, 
			      double *LR_couplings, double *phi,
			      Param p);

void clusterPossibleLR(int i, int *s, int cSpin,
		       bool *Pcluster, std::vector<int> &Rcluster,
		       Param p);

void init_connectivity(Param p);

#endif
