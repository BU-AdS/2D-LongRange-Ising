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
#include "mcPhiFourth2D.h"
#include "phiFourth2D.h"

using namespace std;

extern int seed;
extern mt19937 rng;
extern uniform_real_distribution<double> unif;

/*
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
*/

void PhiFourth2D::correlators(double **ind_corr, double *run_corr, bool temporalDir,
			      int meas, double avePhi, Param p) {
  
  int s_size = p.S1;
  int t_size = p.Lt;
  int idx    = 0;

  double val = 0.0;

  //loop over all sources
  for(int i=0; i<p.surfaceVol; i++) {
    
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


// Improved Correlation Functions.
//--------------------------------

/*
//Overloaded version to handle AdS lattices
void correlatorsImp(double **ind_corr, int meas, double *run_corr,
		    bool dir, vector<Vertex> NodeList, 
		    double avePhi, int *s, Param p) {
  
  double *phi = (double*)malloc(p.S1*p.Lt*sizeof(double));
  int offset = endNode(p.Levels-1,p)+1;
  int disk   = p.AdSVol;  
  for(int i=0; i<p.S1*p.Lt; i++) {    
    phi[i] = NodeList[disk*(i/p.S1) + offset + i%p.S1].phi;
  }
  correlatorsImp(ind_corr, meas, run_corr, dir, phi, avePhi, s, p);
  free(phi);  
}
*/

void PhiFourth2D::correlatorsImpSW(double **ind_corr, double *run_corr,
				   bool temporalDir, int meas, double avePhi,
				   Param p){
  
}

int corr_wc_size = 0;
int corr_wc_ave = 0;
int corr_wc_calls = 0;

void PhiFourth2D::correlatorsImpWolff(double **ind_corr_t, double *run_corr_t,
				      double **ind_corr_s, double *run_corr_s,
				      int meas, double avePhi, Param p){
  
  corr_wc_calls++;
  
  for(int a=0; a<p.surfaceVol; a++) {
    s_cpy[a] = s[a];
    phi_cpy[a] = phi[a];
    cpu_added[a] = false;
  }   
  
  //Choose a random spin.
  int i = int(unif(rng) * p.surfaceVol);
  int cSpin = s_cpy[i];
  
  // The site belongs to the cluster, so flip it.
  corr_wc_size = 1;
  s_cpy[i] *= -1;
  phi_cpy[i] *= -1;
  cpu_added[i] = true;
  
  //This function is recursive and will call itself
  //until all attempts to increase the cluster size
  //have failed.
  corr_wolffClusterAddLR(i, s_cpy, cSpin, LR_couplings, phi_cpy, cpu_added, p);
  
  corr_wc_ave += corr_wc_size;
  
  setprecision(4);
  cout<<"Average (CPU) corr. cluster size at iter "<<meas<<" = "<<corr_wc_ave<<"/"<<corr_wc_calls<<" = "<<1.0*corr_wc_ave/corr_wc_calls<<" = "<<100.0*corr_wc_ave/(corr_wc_calls*p.surfaceVol)<<"%"<<endl;

  int S1 = p.S1;
  int Lt = p.Lt;

  double clusterNorm = 1.0*p.surfaceVol/corr_wc_size;

  double val = 0.0;
  int idx = 0;

  //loop over all sources
  for(int i=0; i<p.surfaceVol; i++) {
    
    if(cpu_added[i]){ 
      
      //loop over all sinks in the specified direction
      for(int j=0; j<Lt; j++) {
	
	idx = abs(i/S1-j);
	if(idx > Lt/2) idx = Lt - idx;
	
	if(cpu_added[i%S1 + j*S1]) {
	  
	  val = ((phi_cpy[i] - avePhi) *
		 (phi_cpy[i%S1 + j*S1] - avePhi))*clusterNorm;
	  
	  ind_corr_t[meas][idx] += val;
	  run_corr_t[idx] += val;
	}
      }
      for(int j=0; j<S1; j++) {
	
	idx = abs(i%S1-j);
	if(idx > S1/2) idx = S1 - idx;
	
	if(cpu_added[S1*(i/S1) + j]) {
	  
	  val = ((phi_cpy[i] - avePhi) *
		 (phi_cpy[S1*(i/S1) + j] - avePhi))*clusterNorm;
	  
	  ind_corr_s[meas][idx] += val;
	  run_corr_s[idx] += val;
	}
      }
    }
  }
}

void Ising2D::correlatorsImpWolffI(double **ind_corr_t, double *run_corr_t,
				   double **ind_corr_s, double *run_corr_s,
				   int meas, double aveS, Param p){
  
  corr_wc_calls++;
  
  for(int a=0; a<p.surfaceVol; a++) {
    s_cpy[a] = s[a];
    cpu_added[a] = false;
  }   
  
  //Choose a random spin.
  int i = int(unif(rng) * p.surfaceVol);
  int cSpin = s_cpy[i];
  
  // The site belongs to the cluster, so flip it.
  corr_wc_size = 1;
  s_cpy[i] *= -1;
  cpu_added[i] = true;
  
  //This function is recursive and will call itself
  //until all attempts to increase the cluster size
  //have failed.
  corr_wolffClusterAddLRI(i, s_cpy, cSpin, LR_couplings, cpu_added, p);
  
  corr_wc_ave += corr_wc_size;
  
  setprecision(4);
  cout<<"Average (CPU) corr. cluster size at iter "<<meas<<" = "<<corr_wc_ave<<"/"<<corr_wc_calls<<" = "<<1.0*corr_wc_ave/corr_wc_calls<<" = "<<100.0*corr_wc_ave/(corr_wc_calls*p.surfaceVol)<<"%"<<endl;

  int S1 = p.S1;
  int Lt = p.Lt;

  double clusterNorm = 1.0*p.surfaceVol/corr_wc_size;

  double val = 0.0;
  int idx = 0;
  int t1,x1,t2,x2,dt,dx;
  
  //loop over all sources
  for(int i=0; i<p.surfaceVol; i++) {
    
    if(cpu_added[i]){ 

      t1 = i / S1;
      x1 = i % S1;
      
      //loop over all sinks 
      for(int j=0; j<p.surfaceVol; j++) {

	if(cpu_added[j]) {
	  
	  //Index divided by circumference, using the int floor feature/bug,
	  //gives the timeslice index.
	  t2 = j / S1;
	  dt = abs(t2-t1) > p.Lt/2 ? p.Lt - abs(t2-t1) : abs(t2-t1);
	  
	  //The index modulo the circumference gives the spatial index.
	  x2 = j % S1;            
	  dx = abs(x2-x1) > p.S1/2 ? p.S1 - abs(x2-x1) : abs(x2-x1);      
	  
	  val = (1 - aveS)*clusterNorm;
	  
	  ind_corr[meas][dx][dt] += val;
	  run_corr[dx][dt] += val;
	}
      }
    }
  }
}

void corr_wolffClusterAddLR(int i, int *s, int cSpin, double *LR_couplings,
			    double *phi, bool *cpu_added, Param p) {
  
  double phi_lc = phi[i];
  int S1 = p.S1;
  int x_len = S1/2 + 1;

  int t1,x1,t2,x2,dt,dx;
  t1 = i / S1;
  x1 = i % S1;
  
  double prob = 0.0;
  double rand = 0.0;
  //We now loop over the possible lattice sites, adding sites
  //(creating bonds) with the specified LR probablity.
  for(int j=0; j<p.surfaceVol; j++) {
    if(s[j] == cSpin && j != i) {
      
      //Index divided by circumference, using the int floor feature/bug,
      //gives the timeslice index.
      t2 = j / p.S1;
      dt = abs(t2-t1) > p.Lt/2 ? p.Lt - abs(t2-t1) : abs(t2-t1);
      
      //The index modulo the circumference gives the spatial index.
      x2 = j % p.S1;            
      dx = abs(x2-x1) > p.S1/2 ? p.S1 - abs(x2-x1) : abs(x2-x1);      
      
      prob = 1 - exp(2*phi_lc*phi[j]*LR_couplings[dx + dt*x_len]);
      rand = unif(rng);
      if(rand < prob) {
	corr_wc_size++;
	// The site belongs to the cluster, so flip it.
	s[j] *= -1;
	phi[j] *= -1;
	cpu_added[j] = true;
	corr_wolffClusterAddLR(j, s, cSpin, LR_couplings, phi, cpu_added, p);
      }
    }
  }
}

void corr_wolffClusterAddLRI(int i, int *s, int cSpin, double *LR_couplings,
			     bool *cpu_added, Param p) {
  
  int S1 = p.S1;
  int x_len = S1/2 + 1;
  double J = p.J;
  
  int t1,x1,t2,x2,dt,dx;
  t1 = i / S1;
  x1 = i % S1;
  
  double prob = 0.0;
  double rand = 0.0;
  //We now loop over the possible lattice sites, adding sites
  //(creating bonds) with the specified LR probablity.
  for(int j=0; j<p.surfaceVol; j++) {
    if(s[j] == cSpin && j != i) {
      
      //Index divided by circumference, using the int floor feature/bug,
      //gives the timeslice index.
      t2 = j / p.S1;
      dt = abs(t2-t1) > p.Lt/2 ? p.Lt - abs(t2-t1) : abs(t2-t1);
      
      //The index modulo the circumference gives the spatial index.
      x2 = j % p.S1;            
      dx = abs(x2-x1) > p.S1/2 ? p.S1 - abs(x2-x1) : abs(x2-x1);      
      
      prob = 1 - exp(-2*J*LR_couplings[dx + dt*x_len]);
      rand = unif(rng);
      if(rand < prob) {
	corr_wc_size++;
	// The site belongs to the cluster, so flip it.
	s[j] *= -1;
	cpu_added[j] = true;
	corr_wolffClusterAddLRI(j, s, cSpin, LR_couplings, cpu_added, p);
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
    if(k<10) {
      cout<<"Autocorrelation "<<k<<" = "<<auto_corr[k]/auto_corr[0]<<endl;
      auto_corr_t -= k/log(auto_corr[k]);
      cout<<"Autocorrelation time at k = "<<1+2*auto_corr_t<<endl;
    }
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
