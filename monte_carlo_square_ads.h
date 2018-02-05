#ifndef MONTE_CARLO_SQUARE_ADS_H
#define MONTE_CARLO_SQUARE_ADS_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "util.h"

//Square lattice, non-local AdS interactions
void clusterAddSqAdS(int i, int *s, int cSpin, bool *Acluster,
		   std::vector<int> &Rcluster, double *LR_couplings,
		   double *phi, Param p);

void clusterPossibleSqAdS(int i, int *s, int cSpin,
			bool *Pcluster, std::vector<int> &Rcluster, Param p);

void wolffUpdateSqAdS(double *phi_arr, int *s, Param p,
		    double *LR_couplings,		  
		    double &delta_mag_phi, int iter);

int metropolisUpdateSqAdS(double *phi_arr, int *s, Param &p,
			double *LR_couplings,
			double & delta_mag_phi, int iter);

double actionSqAdS(double *phi_arr, int *s, Param p,
		 double *LR_couplings,
		 double &KE, double &PE);

void thermaliseSqAdS(double *phi, int *s, Param P,
		   double &delta_mag_phi, double *LR_couplings);

void runMonteCarloSqAdS(double *LR_couplings, std::vector<Vertex> &NodeList,
			Param p);


#endif
