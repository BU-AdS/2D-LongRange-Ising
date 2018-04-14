#ifndef MCISING2D_H
#define MCISING2D_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "ising2D.h"
#include "util.h"

//Short range phi fourth Monte Carlo routines
void wolffUpdateISR(int *s, Param p, int iter);

void wolffClusterAddISR(int i, int *s, int clusterSpin,
		       bool *site_cluster, Param p);

void swendsenWangUpdateISR(int *s, Param p, int iter);

void swendsenWangClusterAddISR(int i, int *s, int cSpin, int clusterNum, 
			      int *clusterDef, Param p);

void metropolisUpdateISR(int *s, Param p, int iter);

double energyISR(int *s, Param p, double &KE);

//Long range phi fourth Monte Carlo routines
void wolffUpdateILR(int *s, Param p, double *LR_couplings, int iter);

void wolffClusterAddILR(int i, int *s, int cSpin, 
		       double *LR_couplings, Param p);

void swendsenWangUpdateILR(int *s, Param p, double *LR_couplings, int iter);

void swendsenWangClusterAddILR(int i, int *s, int cSpin, int clusterNum, 
			      int *clusterDef, std::vector<int> Rcluster, 
			      double *LR_couplings, Param p);

double energyILR(int *s, Param p, double *LR_couplings, double &KE);

void clusterPossibleILR(int i, int *s, int cSpin,
		       bool *Pcluster, std::vector<int> &Rcluster,
		       Param p);

void init_connectivity(Param p);


#endif
