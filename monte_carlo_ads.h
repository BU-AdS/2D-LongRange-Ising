#ifndef MONTE_CARLO_ADS_H
#define MONTE_CARLO_ADS_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "util.h"

//AdS lattice
void clusterAddAdS(int i, int *s, int clusterSpin, bool *site_cluster, 
		   std::vector<Vertex> &NodeList, Param p);

void wolffUpdateAdS(std::vector<Vertex> &NodeList, int *s, 
		    Param p, double &delta_mag_phi, int iter);

int metropolisUpdateAdS(std::vector<Vertex> &NodeList, int *s,
			Param &p, double & delta_mag_phi, int iter);

double actionAdS(std::vector<Vertex> &NodeList, int *s, 
		 Param p, double & KE, double & PE);

void thermaliseAdS(std::vector<Vertex> &NodeList, int *s, Param p, 
		   double &delta_mag_phi);

void runMonteCarloAdS(std::vector<Vertex> &NodeList, Param p);

#endif
