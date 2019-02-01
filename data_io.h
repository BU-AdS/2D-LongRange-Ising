#ifndef DATA_IO_H
#define DATA_IO_H

#include <complex>
#include <cstring>
#include <vector>

#include "util.h"

//Calcuate the Jack-Knife error and write the correlation functions
//to file
void writeObservables(double **ind_corr, double *run_corr, int *norm_corr,
		      double **ind_ft_corr, double *run_ft_corr,
		      int idx, observables obs, Param p);

void writePhi3(double **ind_corr_phi_phi3, double *run_corr_phi_phi3,
	       double **ind_corr_phi3_phi3, double *run_corr_phi3_phi3,
	       double **ind_ft_corr_phi_phi3, double *run_ft_corr_phi_phi3,
	       double **ind_ft_corr_phi3_phi3, double *run_ft_corr_phi3_phi3,
	       int *norm_corr, int idx, observables obs, Param p);

void visualiserIsing(int *s, Param p);
void visualiserIsingCluster(int *s, bool *cluster, Param p);
void visualiserSqr(double *phi_cyl, double barr, Param p);
void visualiserPhi2(double **phi_cyl, Param p, int iter);

#endif
