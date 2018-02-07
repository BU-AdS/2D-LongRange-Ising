#ifndef MONTE_CARLO_SQUARE_NONLOCAL_H
#define MONTE_CARLO_SQUARE_NONLOCAL_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "util.h"

//Square lattice, non-local r^{-(d+\sigma)} interactions
void clusterAddSqNL(int i, int *s, int cSpin, bool *Acluster,
		    std::vector<int> &Rcluster, double *LR_couplings,
		    double *phi, Param p);

void clusterPossibleSqNL(int i, int *s, int cSpin,
			 bool *Pcluster, std::vector<int> &Rcluster, Param p);

void wolffUpdateSqNL(double *phi_arr, int *s, Param p,
		     double *LR_couplings,		  
		     double &delta_mag_phi, int iter);

void swendsenWangUpdateSqNL(double *phi_arr, int *s, Param p,
			    double *LR_coupling,
			    double &delta_mag_phi, int iter);

void SWclusterAddSqNL(int i, int *s, int cSpin, int clusterNum, 
		      int *clusterDef, std::vector<int> Rcluster, 
		      double *LR_couplings,
		      double *phi_arr, Param p);

int metropolisUpdateSqNL(double *phi_arr, int *s, Param &p,
			 double *LR_couplings,
			 double & delta_mag_phi, int iter);

double actionSqNL(double *phi_arr, int *s, Param p,
		  double *LR_couplings,
		  double &KE, double &PE);

void thermaliseSqNL(double *phi, int *s, Param P,
		    double &delta_mag_phi, double *LR_couplings);

void runMonteCarloSqNL(std::vector<Vertex> &NodeList, Param p);

#endif
