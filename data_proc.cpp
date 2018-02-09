#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>

#include "util.h"
#include "hyp_util.h"
#include "data_proc.h"

using namespace std;

//Overloaded version to handle AdS lattices
void correlators(double **ind_corr, int meas, double *run_corr, bool dir,
		 vector<Vertex> NodeList, 
		 double avePhi, Param p) {
  
  double *phi = (double*)malloc(p.S1*p.Lt*sizeof(double));
  int offset = endNode(p.Levels-1,p)+1;
  int disk   = p.AdSVol;  
  for(int i=0; i<p.S1*p.Lt; i++) {    
    phi[i] = NodeList[disk*(i/p.S1) + offset + i%p.S1].phi;
  }
  correlators(ind_corr, meas, run_corr, dir, phi, avePhi, p);
  free(phi);  
}


void correlators(double **ind_corr, int meas, double *run_corr, 
		 bool temporalDir, double *phi, double avePhi, Param p) {
  
  int s_size = p.S1;
  int t_size = p.Lt;
  int vol    = p.surfaceVol;
  int idx    = 0;

  double val = 0.0;

  //loop over all sources
  for(int i=0; i<vol; i++) {

    //loop over all sinks in the specified direction
    if(temporalDir) {
      for(int j=0; j<t_size; j++) {
	idx = abs(i/s_size-j);
	if(idx > t_size/2) idx = t_size - idx;
	val = ((phi[i] - avePhi) *
	       (phi[i%s_size + j*s_size] - avePhi));

	ind_corr[meas][idx] += val;
	run_corr[idx] += val;
      }
    } else {
      for(int j=0; j<s_size; j++) {
	
	idx = abs(i%s_size-j);
	if(idx > s_size/2) idx = s_size - idx;
	
	val = ((phi[i] - avePhi) *
	       (phi[s_size*(i/s_size) + j] - avePhi));
	
	ind_corr[meas][idx] += val;
	run_corr[idx] += val;
      }
    }
  }
}

//Calculate the autocorrelation of |phi|
void autocorrelation(double *PhiAb_arr, double avePhiAbs, 
		     int meas, double *auto_corr) {

  double auto_corr_t = 0;  
  //Calculate the variance function
  for(int i=0; i<meas; i++)
    auto_corr[0] += pow(PhiAb_arr[i] - avePhiAbs,2)/meas;
  
  //Calculate the ratios Ck/C0
  for(int k=1; k<meas; k++) {
    for(int i=0; i<k; i++) {
      auto_corr[k] += (PhiAb_arr[i]*PhiAb_arr[i+k] - 
		       avePhiAbs*avePhiAbs)/(meas-k);
    }
    cout<<"Autocorrelation "<<k<<" = "<<auto_corr[k]/auto_corr[0]<<endl;
    if(k<100) auto_corr_t -= k/log(auto_corr[k]);
    cout<<"Autocorrelation time at k = "<<1+2*auto_corr_t<<endl;
  }

}

//Data is entered in 'raw form.' normalisation of the data is done
//at the data write stage.
void jackknife(double **ind, double *run, double *jk_err, int block, 
	       int data_points, int arr_length, Param p) {
  
  int num_resamp = data_points/block;
  double coeff = (1.0*num_resamp - 1.0)/(1.0*num_resamp);
  double resamp_norm = 1.0/(data_points - block);
  
  double **resamp_ave = (double**)malloc(arr_length*sizeof(double*));
  for(int r=0; r<arr_length; r++) {
    jk_err[r] = 0.0;
    resamp_ave[r] = (double*)malloc(num_resamp*sizeof(double));
    for(int i=0; i<num_resamp; i++) resamp_ave[r][i] = 0.0;
  }
  
  //Get resampling averages
  for(int r=0; r<arr_length; r++) {
    for(int i=0; i<num_resamp; i++) {
      for(int j=0; j<data_points; j++) {   
	if(j < i*block || j >= (i+1)*block)  {
	  resamp_ave[r][i] += ind[j][r];
	}
      }
      resamp_ave[r][i] *= resamp_norm;
    }
  }
  
  //Get jk error
  //The run/data_ponts value is the arithmetic mean.
  for(int r=0; r<arr_length; r++) {
    for(int i=0; i<num_resamp; i++) {
      jk_err[r] += pow(run[r]/data_points - resamp_ave[r][i],2);
    }
    jk_err[r] = sqrt(coeff*jk_err[r]);    
  }  
}
