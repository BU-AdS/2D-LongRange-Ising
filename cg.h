#ifndef CG_H 
#define CG_H 
#include <complex>
#include <cstring>
#include <vector>

#include "util.h"

void Mphi(double *phi, const double *phi0, 
	  const std::vector<Vertex> NodeList, Param P);

double Minv_phi(double *phi, double *b, 
		const std::vector<Vertex> NodeList, Param P);

//This function will calculate the M(x,x)^{-1} elements on the qth
//portion of the graph. 
void oneLoopCorrection(std::vector<Vertex> &NodeList, Param& p);

#endif
