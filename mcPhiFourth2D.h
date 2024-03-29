#ifndef MCPHIFOURTH2D_H
#define MCPHIFOURTH2D_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "phiFourth2D.h"
#include "util.h"


//Short range phi fourth Monte Carlo routines
void wolffUpdateSR(double *phi, int *s, Param p, int iter);

void wolffClusterAddSR(int i, int *s, int clusterSpin,
		       bool *site_cluster, double *phi, Param p);

void swendsenWangUpdateSR(double *phi, int *s, Param p, int iter);

void swendsenWangClusterAddSR(int i, int *s, int cSpin, int clusterNum, 
			      int *clusterDef, double *phi, Param p);

void metropolisUpdateSR(double *phi, int *s, Param p, int iter);

double actionSR(double *phi, int *s, Param p, double &KE, double &PE);

//Long range phi fourth Monte Carlo routines
void wolffUpdateLR(double *phi, int *s, Param p,
		   double *LR_couplings, int iter, int n);

void wolffClusterAddLR(int *spinStack, int *spinStackLR, int *s,
		       int &stackSize, int cSpin, 
		       double *LR_couplings, double *phi, Param p);

void swendsenWangUpdateLR(double *phi, int *s, Param p,
			  double *LR_couplings, int iter);

void swendsenWangClusterAddLR(int i, int *s, int cSpin, int clusterNum, 
			      int *clusterDef, std::vector<int> Rcluster, 
			      double *LR_couplings, double *phi,
			      Param p);

double actionLR(double *phi_arr, int *s, Param p,
		double *LR_couplings,
		double &KE, double &PE);

void clusterPossibleLR(int i, int *s, int cSpin,
		       bool *Pcluster, std::vector<int> &Rcluster,
		       Param p);

void init_connectivityPhi4(Param p);


#endif
