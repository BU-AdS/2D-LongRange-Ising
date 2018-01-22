#ifndef UPDATE_H
#define UPDATE_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "util.h"

void clusterAddSqr(int i, std::vector<int> &s, int clusterSpin,
		   bool *site_cluster, std::vector<double> &phi_arr, Param p);

void clusterAddSqr(int i, std::vector<int> &s, int clusterSpin,
		   bool *cluster, std::vector<double> &phi_arr, Param p);

void wolffUpdatePhiSqr(std::vector<double> &phi_arr, std::vector<int> &s, Param p,
		       double &delta_mag_phi, int iter);

double actionPhiSqr(std::vector<double> &phi_arr, std::vector<int> &s, Param p,
		    double & KE, double & PE);

int metropolisUpdatePhiSqr(std::vector<double> &phi_arr, std::vector<int> &s, Param &p,
			   double & delta_mag_phi, int iter);


void clusterAddAdS(int i, std::vector<int> &s, int clusterSpin,
		   bool *cluster, std::vector<Vertex> &NodeList, Param p);

void clusterAddAdS(int i, std::vector<int> &s, int clusterSpin,
		   bool *cluster, std::vector<Vertex> &NodeList, Param p);

void wolffUpdatePhiAdS(std::vector<Vertex> &NodeList, std::vector<int> &s, Param p,
		       double &delta_mag_phi, int iter);

double actionPhiAdS(std::vector<Vertex> &NodeList, std::vector<int> &s, Param p,
		    double & KE, double & PE);

int metropolisUpdatePhiAdS(std::vector<Vertex> &NodeList, std::vector<int> &s, Param &p,
			   double & delta_mag_phi, int iter);

void runMonteCarloAdS(std::vector<Vertex> &NodeList, Param p);
void runMonteCarloSqr(std::vector<Vertex> &NodeList, Param p);

#endif
