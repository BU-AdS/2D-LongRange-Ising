#ifndef DATA_IO_H
#define DATA_IO_H

#include <complex>
#include <cstring>
#include <vector>

#include "util.h"
#include "hyp_util.h"

//Calcuate the Jack-Knife error and write the correlation functions
//to file
void writeObservables(double **ind_corr, double *run_corr, int *norm_corr,
		      double **ind_ft_corr, double *run_ft_corr,
		      int idx, observables obs, Param p);

void PrintNodeTables(const std::vector<Vertex> NodeList, Param P);
void PrintComplexPositions(const std::vector<Vertex> NodeList, Param P);

//One-Size-Fits-All data file for lattice/analytical propagator data.
void dataDump(std::vector<Vertex> NodeList, double *phi, Param p, int level,
	      int t_range, int shift);

void visualiserIsing(int *s, Param p);
void visualiserIsingCluster(int *s, bool *cluster, Param p);
void visualiserSqr(double *phi_cyl, double barr, Param p);
void visualiserAdS(std::vector<Vertex> NodeList, double barr, Param p);
void visualiserPhi2(double **phi_cyl, Param p, int iter);

#endif
