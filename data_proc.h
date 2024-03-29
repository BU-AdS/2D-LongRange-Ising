#ifndef DATA_PROC_H
#define DATA_PROC_H


// Functions to process data on the fly
//-------------------------------------

void correlators(double **ind_corr, int meas, double *run_corr, bool dir,
		 double *phi, double avePhi, Param p);

void correlatorsImp(double **ind_corr, int meas, double *run_corr, 
		    bool dir, double *phi, double avePhi, int *s, Param p);

void FTcorrelation(double **ind_ft_corr, double *run_ft_corr,
		   double **ind_corr, int *norm_corr, 
		   int meas, Param p);

void autocorrelation(double *PhiAb_arr, double avePhiAbs, 
		     int meas, double *auto_corr);

void jackknife(double **ind, double *run, double *jk_err, int block, 
	       int data_points, int arr_length, int dth, Param p);

double jackknifeVar(double *x, double *xsq, double val, int block,
		    int data_points, double J);
double jackknifeBin(double *x, double *xsq, double val, int block,
		    int data_points);
void jackknifeFT(double **ind_ft_corr, double *run_ft_corr, double *jk_err,
		   int block, int data_points, int l, Param p);

void corr_wolffClusterAddLR(int i, int *s, int cSpin, double *LR_couplings,
			    double *phi, bool *cluster, Param p);

void corr_wolffClusterAddLRI(int i, int *s, int cSpin, double *LR_couplings,
			     bool *cluster, Param p);

void corr_wolffClusterAddSRI(int i, int *s, int cSpin, double prob,
                             bool *cpu_added, Param p);


#endif
