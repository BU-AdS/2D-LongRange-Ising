#ifndef MONTE_CARLO_ADS_LOCAL_H
#define MONTE_CARLO_ADS_LOCAL_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "util.h"

//AdS lattice
void clusterAddAdSL(int i, int *s, int clusterSpin, bool *site_cluster, 
		   std::vector<Vertex> &NodeList, Param p);

void wolffUpdateAdSL(std::vector<Vertex> &NodeList, int *s, 
		    Param p, double &delta_mag_phi, int iter);

int metropolisUpdateAdSL(std::vector<Vertex> &NodeList, int *s,
			Param &p, double & delta_mag_phi, int iter);

double actionAdSL(std::vector<Vertex> &NodeList, int *s, 
		 Param p, double & KE, double & PE);

void thermaliseAdSL(std::vector<Vertex> &NodeList, int *s, Param p, 
		   double &delta_mag_phi);

void runMonteCarloAdSL(std::vector<Vertex> &NodeList, Param p);

#endif
