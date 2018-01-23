#ifndef UPDATE_H
#define UPDATE_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "util.h"

//Square lattice
void clusterAddSqr(int i, int *s, int clusterSpin,
		   bool *site_cluster, double *phi_arr, Param p);

void wolffUpdateSqr(double *phi_arr, int *s, Param p,
		    double &delta_mag_phi, int iter);

int metropolisUpdateSqr(double *phi_arr, int *s, 
			Param &p, double & delta_mag_phi, int iter);

double actionSqr(double *phi_arr, int *s, 
		 Param p, double & KE, double & PE);

void thermaliseSqr(double *phi, int *s, Param P, double &delta_mag_phi);

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

void runMonteCarlo(std::vector<Vertex> &NodeList, Param p);



#endif
