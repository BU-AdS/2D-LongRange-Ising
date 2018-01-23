#ifndef DATA_PROC_H
#define DATA_PROC_H


// Functions to process data on the fly
//-------------------------------------

void correlators(double **corr, double **corr_ave, int corr_norm,
                 std::vector<Vertex> NodeList, double avePhi, Param p);

void correlators(double **corr, double **corr_ave, int corr_norm,
                 double *phi, double avePhi, Param p);

void autocorrelation(double *PhiAb_arr, double avePhiAbs, int meas);

void jackknife(double *array, int block, int length);

#endif
