#ifndef CG_H 
#define CG_H 
#include <complex>
#include <cstring>
#include <vector>

#include "util.h"

//-------------------//
// Lin-Alg utilities //
//-------------------//

// v1 = 0
void zero_vector(double* v1, const int size);

// v1 = v2
void copy_vector(double* v1, double* v2, const int size);

// v2 += alpha v1
void caxpy(const double alpha, double* v1, double* v2, const int size);

// v2 = v1 + alpha * v2
void cxpay(double* v1, const double alpha, double* v2, const int size);

// v2 = alpha v1 + beta v2
void caxpby(const double alpha, double* v1, const double beta, double* v2, 
	    const int size);

// v1 dot v2
double dot(double* v1, double* v2, const int size);

// ||v1||^2
double norm2sq(double* v1, const int size);

// ||v1 - v2||^2
double diffnorm2sq(double* v1, double* v2, const int size);

//--------------------------------//
// Linear operators and inverters //
//--------------------------------//

void Mphi(double *phi, const double *phi0, 
	  const std::vector<Vertex> NodeList, Param P);

double Minv_phi(double *phi, double *b, 
		const std::vector<Vertex> NodeList, Param P);

//Wrapper for AdS code.
void Minv_phi_ms(double **phi, double *phi0, std::vector<Vertex> NodeList, 
		 Param p);

// Solves lhs = A^(-1) rhs using multishift CG as defined in
// http://arxiv.org/pdf/hep-lat/9612014.pdf
// Assumes there are n_shift values in "shifts".
// If they are sorted s.t. the smallest shift is the smallest
// (worst-conditioned solve), set worst_first = true. 
// resid_freq_check is how often to check the residual of other solutions.
// This lets us stop iterating on converged systems. 
void cg_multishift(double **phi, double *phi0, int n_shift, int size,
		   int resid_freq_check, int max_iter, double eps,
		   double* shifts, std::vector<Vertex> NodeList,
		   Param param);

//---------------------//
// AdS suite functions //
//---------------------//

//This function will calculate the M(x,x)^{-1} elements on the qth
//portion of the graph. From this information we can derive the
//one loop corrections, and the LR couplings.
void oneLoopCorrection(double *LR_couplings,
		       std::vector<Vertex> &NodeList,
		       Param &p);

//This function will calculate the scaling factors needed to 
//make the AdS lattice action match the analytcal action.
void latticeScaling(std::vector<Vertex> &NodeList, Param& p);

#endif
