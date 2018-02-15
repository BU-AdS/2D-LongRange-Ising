#ifndef MONTE_CARLO_SQUARE_LOCAL_H
#define MONTE_CARLO_SQUARE_LOCAL_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "mc_util.h"
#include "util.h"


//Cluster and Metropolis routines.
void wolffUpdateSR(double *phi_arr, int *s, Param p,
		   double &delta_mag_phi, int iter);

void wolffClusterAddSR(int i, int *s, int clusterSpin,
		       bool *site_cluster, double *phi_arr, Param p);


void swendsenWangUpdateSR(double *phi_arr, int *s, Param p,
			  double &delta_mag_phi, int iter);

void swendsenWangClusterAddSR(int i, int *s, int cSpin, int clusterNum, 
			      int *clusterDef, double *phi_arr, Param p);

int metropolisUpdateSR(double *phi_arr, int *s, 
		       Param &p, double & delta_mag_phi, int iter);

//void runMonteCarloSqL(Param p);

#endif
