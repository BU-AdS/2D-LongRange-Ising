#ifndef MC_UPDATE_SR_H
#define MC_UPDATE_SR_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "mc_util.h"
#include "util.h"


//Cluster and Metropolis routines.
void wolffUpdateSR(double *phi, int *s, Param p, int iter);

void wolffClusterAddSR(int i, int *s, int clusterSpin,
		       bool *site_cluster, double *phi, Param p);

void swendsenWangUpdateSR(double *phi, int *s, Param p, int iter);

void swendsenWangClusterAddSR(int i, int *s, int cSpin, int clusterNum, 
			      int *clusterDef, double *phi, Param p);

int metropolisUpdateSR(double *phi, int *s, Param p, int iter);

#endif
