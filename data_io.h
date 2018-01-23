#ifndef DATA_IO_H
#define DATA_IO_H

#include <complex>
#include <cstring>
#include <vector>

#include "util.h"
#include "hyp_util.h"

// Functions to print, dump, or visualise data.
//----------------------------------------------
void PrintNodeTables(const std::vector<Vertex> NodeList, Param P);
void PrintComplexPositions(const std::vector<Vertex> NodeList, Param P);

//One-Size-Fits-All data file for lattice/analytical propagator data.
void dataDump(std::vector<Vertex> NodeList, double *phi, Param p, int level,
	      int t_range, int shift);

void visualiserSqr(std::vector<double> phi_cyl, double barr, Param p);
void visualiserAdS(std::vector<Vertex> NodeList, double barr, Param p);
void visualiserPhi2(double **phi_cyl, Param p, int iter);

#endif
