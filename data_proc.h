#ifndef DATA_PROC_H
#define DATA_PROC_H


// Functions to process data on the fly
//-------------------------------------

void correlators(double **ind_corr, int meas, double *run_corr, bool dir,
		 std::vector<Vertex> NodeList, 
		 double avePhi, Param p);

void correlators(double **ind_corr, int meas, double *run_corr, bool dir,
		 double *phi, double avePhi, Param p);

void correlatorsImp(double **ind_corr, int meas, double *run_corr, 
		    bool dir, std::vector<Vertex> NodeList, 
		    double avePhi, int *s, Param p);

void correlatorsImp(double **ind_corr, int meas, double *run_corr, 
		    bool dir, double *phi, double avePhi, int *s, Param p);

void autocorrelation(double *PhiAb_arr, double avePhiAbs, 
		     int meas, double *auto_corr);

void jackknife(double **ind, double *run, double *jk_err, int block, 
	       int data_points, int arr_length, Param p);

void corr_wolffClusterAddLR(int i, int *s, int cSpin, double *LR_couplings,
			    double *phi, bool *cluster, Param p);

#endif
