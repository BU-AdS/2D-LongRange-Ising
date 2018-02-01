#ifndef MONTE_CARLO_SQUARE_LOCAL_H
#define MONTE_CARLO_SQUARE_LOCAL_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "util.h"

//Square lattice
void clusterAddSqL(int i, int *s, int clusterSpin,
		   bool *site_cluster, double *phi_arr, Param p);

void wolffUpdateSqL(double *phi_arr, int *s, Param p,
		    double &delta_mag_phi, int iter);

int metropolisUpdateSqL(double *phi_arr, int *s, 
			Param &p, double & delta_mag_phi, int iter);

double actionSqL(double *phi_arr, int *s, 
		 Param p, double & KE, double & PE);

void thermaliseSqL(double *phi, int *s, Param P, double &delta_mag_phi);

void runMonteCarloSqL(std::vector<Vertex> &NodeList, Param p);

#endif
