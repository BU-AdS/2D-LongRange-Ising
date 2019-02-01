#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>

#include "util.h"
#include "data_proc.h"
#include "data_io.h"
#include "mcPhiFourth2D.h"
#include "phiFourth2D.h"

using namespace std;

extern int seed;
extern mt19937 rng;
extern uniform_real_distribution<double> unif;

//Basic utilites
// x = 0,1,..., p.S1-1  and
// t = 0,1,..., p.Lt-1
// i = x + p.S1*t so
// x = i % p.S1 and
// t = i/p.S1 
inline int xp(int i, Param p) {
  //if( (i+1)%p.S1 == 0 ) return (i/p.S1)*p.S1;
  //else return i+1;
  return (i + 1) % p.S1 + p.S1 * (i / p.S1);
}

inline int xm(int i, Param p) {
  //if( i%p.S1 == 0 ) return (i/p.S1)*p.S1 + (p.S1-1);
  //else return i-1;  
  return (i - 1 + p.S1) % p.S1 + p.S1 * (i / p.S1);
}

inline int tp(int i, Param p) {
  //if(i/p.S1 == p.Lt-1) return (i - (p.Lt-1)*p.S1);
  //else return i+p.S1;
  return i % p.S1 + p.S1 * ((i / p.S1 + 1) % p.Lt);
}

inline int ttm(int i, Param p){
  //if(i/p.S1 == 0) return (i + (p.Lt-1)*p.S1);
  //else return i-p.S1;
  return i % p.S1 + p.S1 * ((i / p.S1 - 1 + p.Lt) % p.Lt);
}


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

void PhiFourth2D::correlators(int meas, double avePhi, Param p) {
  
  int S1 = p.S1;
  int vol = p.surfaceVol;
  int x_len = S1/2 + 1;
  int idx = 0;
  double val = 0.0;
  double phi_lc = 0.0;

  int t1,x1,t2,x2,dt,dx;
  
  //loop over all sources
  for(int i=0; i<vol; i++) {

    t1 = i/S1;
    x1 = i%S1;
    phi_lc = phi[i];
    
    //loop over all sinks.
    for(int j=0; j<vol; j++) {
      
      //Index divided by circumference, using the int floor feature/bug,
      //gives the timeslice index.
      t2 = j / S1;
      dt = abs(t2-t1) > p.Lt/2 ? p.Lt - abs(t2-t1) : abs(t2-t1);
      
      //The index modulo the circumference gives the spatial index.
      x2 = j % S1;            
      dx = abs(x2-x1) > p.S1/2 ? p.S1 - abs(x2-x1) : abs(x2-x1);      

      idx = dx + x_len*dt;
      
      val = (phi_lc*phi[j] - avePhi*avePhi);

      ind_corr[meas][idx] += val;
      run_corr[idx] += val;
    }
  }
}

void PhiFourth2D::correlatorsPhi3(int meas, double avePhi, double phi3Ave, Param p) {
  
  int S1 = p.S1;
  int vol = p.surfaceVol;
  int x_len = S1/2 + 1;
  int idx = 0;
  double val1 = 0.0;
  double val2 = 0.0;
  double phi3_lc = 0.0;
  double phi_lc = 0.0;
  int t1,x1,t2,x2,dt,dx;
  
  //loop over all sources, cube the source and use symmetry.
  for(int i=0; i<vol; i++) {

    t1 = i / S1;
    x1 = i % S1;
    phi_lc  = phi[i]; 
    phi3_lc = pow(phi_lc,3);
    
    //loop over all sinks. 
    for(int j=0; j<vol; j++) {
      
      //Index divided by circumference, using the int floor feature/bug,
      //gives the timeslice index.
      t2 = j / S1;
      dt = abs(t2-t1) > p.Lt/2 ? p.Lt - abs(t2-t1) : abs(t2-t1);
      
      //The index modulo the circumference gives the spatial index.
      x2 = j % S1;            
      dx = abs(x2-x1) > p.S1/2 ? p.S1 - abs(x2-x1) : abs(x2-x1);      

      idx = dx + x_len*dt;
      
      val1 = (phi3_lc*phi[j] - phi3Ave*avePhi);
      val2 = (phi3_lc*phi[j]*phi[j]*phi[j] - phi3Ave);
      
      ind_corr_phi_phi3[meas][idx] += val1;
      run_corr_phi_phi3[idx] += val1;
      
      ind_corr_phi3_phi3[meas][idx] += val2;
      run_corr_phi3_phi3[idx] += val2;
      
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

void PhiFourth2D::correlatorsImpSW(int meas, double avePhi, Param p){
  cout<<"Error: PhiFourth2D::correlatorsImpSW not implemented."<<endl;
  exit(0);
}

int corr_wc_size = 0;
int corr_wc_ave = 0;
int corr_wc_calls = 0;

void PhiFourth2D::correlatorsImpWolff(int meas, double avePhi, Param p){
  
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
  cout<<"Average (CPU) Phi4 Wolff correlation cluster size at measurement "<<meas+1<<" = "<<corr_wc_ave<<"/"<<corr_wc_calls<<" = "<<(1.0*corr_wc_ave)/corr_wc_calls<<" = "<<(100.0*corr_wc_ave)/(corr_wc_calls*p.surfaceVol)<<"%"<<endl;

  int S1 = p.S1;
  int Lt = p.Lt;
  int vol = p.surfaceVol;
  int x_len = S1/2 + 1;
  int idx = 0;
  double val = 0.0;
  double phi_lc = 0.0;
  int t1,x1,t2,x2,dt,dx;
  double clusterNorm = 1.0*p.surfaceVol/corr_wc_size;

  //loop over all sources
  for(int i=0; i<vol; i++) {

    if(cpu_added[i]){ 
      t1 = i/S1;
      x1 = i%S1;
      phi_lc = phi[i];
      
      //loop over all sinks.
      for(int j=0; j<vol; j++) {

	if(cpu_added[j]) {
	
	  //Index divided by circumference, using the int floor feature/bug,
	  //gives the timeslice index.
	  t2 = j / S1;
	  dt = abs(t2-t1) > Lt/2 ? Lt - abs(t2-t1) : abs(t2-t1);
	  
	  //The index modulo the circumference gives the spatial index.
	  x2 = j % S1;            
	  dx = abs(x2-x1) > S1/2 ? S1 - abs(x2-x1) : abs(x2-x1);      
	  
	  idx = dx + x_len*dt;
	  
	  val = (phi_lc*phi[j] - avePhi*avePhi)*clusterNorm;

	  ind_corr[meas][idx] += val;
	  run_corr[idx] += val;
	}
      }
    }
  }
}

void Ising2D::correlatorsImpWolffI(int meas, double aveS, Param p){
  
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

  if(p.coupling_type == SR) {
    double prob = 1 - exp(-2*p.J);
    corr_wolffClusterAddSRI(i, s_cpy, cSpin, prob, cpu_added, p);
  }
  else corr_wolffClusterAddLRI(i, s_cpy, cSpin, isingProb, cpu_added, p);
  
  corr_wc_ave += corr_wc_size;
  
  setprecision(4);
  cout<<"Average (CPU) Ising Wolff correlation cluster size at measurement "<<meas+1<<" = "<<corr_wc_ave<<"/"<<corr_wc_calls<<" = "<<(1.0*corr_wc_ave)/(corr_wc_calls)<<" = "<<(100.0*corr_wc_ave)/(corr_wc_calls*p.surfaceVol)<<"%"<<endl;

  //visualiserIsingCluster(s, cpu_added, p);
  
  int S1 = p.S1;
  int Lt = p.Lt;
  int vol = p.surfaceVol;
  int x_len = S1/2 + 1;
  int idx = 0;
  int t1,x1,t2,x2,dt,dx;
  double clusterNorm = 1.0*p.surfaceVol/corr_wc_size;
  
  //loop over all sources
  for(int i=0; i<p.surfaceVol; i++) {
    
    if(cpu_added[i]){ 

      t1 = i / S1;
      x1 = i % S1;
      
      //loop over all sinks 
      for(int j=0; j<vol; j++) {

	if(cpu_added[j]) {
	  
	  //Index divided by circumference, using the int floor feature/bug,
	  //gives the timeslice index.
	  t2 = j / S1;
	  dt = abs(t2-t1) > Lt/2 ? Lt - abs(t2-t1) : abs(t2-t1);
	  
	  //The index modulo the circumference gives the spatial index.
	  x2 = j % S1;            
	  dx = abs(x2-x1) > S1/2 ? S1 - abs(x2-x1) : abs(x2-x1);      

	  idx = dx + x_len*dt;

	  ind_corr[meas][idx] += clusterNorm;
	  run_corr[idx] += clusterNorm;
	}
      }
    }
  }
}

void corr_wolffClusterAddLR(int i, int *s, int cSpin, double *LR_couplings,
			    double *phi, bool *cpu_added, Param p) {
  
  double phi_lc = phi[i];
  int S1 = p.S1;
  int Lt = p.Lt;
  int x_len = S1/2 + 1;

  int t1,x1,t2,x2,dt,dx;
  t1 = i / S1;
  x1 = i % S1;
  
  double prob = 0.0;
  double rand = 0.0;
  //We now loop over the possible lattice sites, adding sites
  //(creating bonds) with the specified LR probablity.
  for(int j=0; j<p.surfaceVol; j++) {
    if(s[j] == cSpin ) {
      
      //Index divided by circumference, using the int floor feature/bug,
      //gives the timeslice index.
      t2 = j / p.S1;
      dt = abs(t2-t1) > Lt/2 ? Lt - abs(t2-t1) : abs(t2-t1);
      
      //The index modulo the circumference gives the spatial index.
      x2 = j % p.S1;            
      dx = abs(x2-x1) > S1/2 ? S1 - abs(x2-x1) : abs(x2-x1);      
      
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

void corr_wolffClusterAddLRI(int i, int *s, int cSpin, double *isingProb,
			     bool *cpu_added, Param p) {
  
  int S1 = p.S1;
  int Lt = p.Lt;
  int x_len = S1/2 + 1;
  
  int t1,x1,t2,x2,dt,dx;
  t1 = i / S1;
  x1 = i % S1;
  
  double prob = 0.0;
  double rand = 0.0;
  //We now loop over the possible lattice sites, adding sites
  //(creating bonds) with the specified LR probablity.
  for(int j=0; j<p.surfaceVol; j++) {
    if(s[j] == cSpin ) {
      
      //Index divided by circumference, using the int floor feature/bug,
      //gives the timeslice index.
      t2 = j / p.S1;
      dt = abs(t2-t1) > Lt/2 ? Lt - abs(t2-t1) : abs(t2-t1);
      
      //The index modulo the circumference gives the spatial index.
      x2 = j % p.S1;            
      dx = abs(x2-x1) > S1/2 ? S1 - abs(x2-x1) : abs(x2-x1);      
      
      //prob = 1 - exp(-2*p.J*LR_couplings[dx + dt*x_len]);
      prob = isingProb[dx + dt*x_len];      
      rand = unif(rng);
      
      if(rand < prob) {
	corr_wc_size++;
	// The site belongs to the cluster, so flip it.
	s[j] *= -1;
	cpu_added[j] = true;
	corr_wolffClusterAddLRI(j, s, cSpin, isingProb, cpu_added, p);
	//corr_wolffClusterAddLRI(j, s, cSpin, LR_couplings, cpu_added, p);
      }
    }
  }
}

void corr_wolffClusterAddSRI(int i, int *s, int cSpin, double prob,
			     bool *cpu_added, Param p) {
  
  // The site belongs to the cluster, so flip it.
  cpu_added[i] = true;
  s[i] *= -1;

  // - If the (aligned) neighbour spin does not already belong to the
  // cluster, then try to add it to the cluster.
  // - If the site has already been added, then we may skip the test.

  //Forward in T
  if(!cpu_added[ tp(i,p) ] && s[tp(i,p)] == cSpin) {
    if(unif(rng) < prob) {
      corr_wc_size++;
      //cout<<"->tp";
      corr_wolffClusterAddSRI(tp(i,p), s, cSpin, prob, cpu_added, p);
    }
  }

  //Forward in X
  if(!cpu_added[ xp(i,p) ] && s[xp(i,p)] == cSpin) {
    if(unif(rng) < prob) {
      corr_wc_size++;
      //cout<<"->xp";
      corr_wolffClusterAddSRI(xp(i,p), s, cSpin, prob, cpu_added, p);
    }
  }

  
  //Backard in T 
  if(!cpu_added[ ttm(i,p) ] && s[ttm(i,p)] == cSpin) {  
    if(unif(rng) < prob) {
      corr_wc_size++;
      //cout<<"->tm";
      corr_wolffClusterAddSRI(ttm(i,p), s, cSpin, prob, cpu_added, p);
    }
  }

  //Backward in X
  if(!cpu_added[ xm(i,p) ] && s[xm(i,p)] == cSpin) {
    if (unif(rng) < prob) {
      corr_wc_size++;
      //cout<<"->xm";
      corr_wolffClusterAddSRI(xm(i,p), s, cSpin, prob, cpu_added, p);
    }
  }   
}

//Fourier transfrom the correlation functions
void FTcorrelation(double **ind_ft_corr, double *run_ft_corr,
		   double **ind_corr, int *norm_corr, 
		   int meas, Param p) {

  int x_len = p.S1/2 + 1;
  
  //loop over l projection
  for(int l=0; l<3; l++) {
    //loop over time slices
    for(int dt = 0; dt <p.Lt/2+1; dt++) {
      //loop over spatial index
      for(int dx = 0; dx<p.S1/2+1; dx++) {
	ind_ft_corr[3*dt+l][meas] += (cos(l*dx*2.0*M_PI/p.S1)*ind_corr[meas][dx + dt*x_len]/norm_corr[dx + dt*x_len]);
	//cout<<cos(l*dx*2.0*M_PI/p.S1)<<" "<<ind_corr[meas][dx + dt*(p.S1/2+1)]<<" ";
      }
      ind_ft_corr[3*dt+l][meas] /= p.S1;
      run_ft_corr[3*dt+l] += ind_ft_corr[3*dt+l][meas];
      //cout<<run_ft_corr[3*dt+l]<<" "<<ind_ft_corr[3*dt+l][meas]<<endl;
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
	       int data_points, int arr_length, int dth, Param p) {

  //Initialise
  int x_len = p.S1/2+1;
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
	  resamp_ave[r][i] += ind[j][dth + i*x_len];
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

double jackknifeVar(double *x, double *xsq, double val, int block, int data_points, double J) {
  
  //Initialise
  int num_resamp = data_points/block;
  double coeff = (1.0*num_resamp - 1.0)/(1.0*num_resamp);
  double resamp_norm = 1.0/(data_points - block);  
  double *resamp_x = (double*)malloc(num_resamp*sizeof(double));
  double *resamp_xsq = (double*)malloc(num_resamp*sizeof(double));
  double *resamp_val = (double*)malloc(num_resamp*sizeof(double));
  double jk_err = 0.0;
  for(int a=0; a<num_resamp; a++) {
    resamp_x[a]   = 0.0;
    resamp_xsq[a] = 0.0;
    resamp_val[a] = 0.0;
  }

  //Get resampled averages
  for(int i=0; i<num_resamp; i++) {
    for(int j=0; j<data_points; j++) {   
      if(j < i*block || j >= (i+1)*block)  {
	resamp_x[i] += x[j];
	resamp_xsq[i] += xsq[j];
      }
    }
    resamp_x[i] *= resamp_norm;
    resamp_xsq[i] *= resamp_norm;
  }

  //Use resampled averages to get resampled values of <E(x)^2> - <E(x)>^2.
  for(int i=0; i<num_resamp; i++) {   
    resamp_val[i] = (resamp_xsq[i] - resamp_x[i]*resamp_x[i]);
  }
  
  //Get jk_error
  for(int i=0; i<num_resamp; i++) {
    jk_err += pow(val - resamp_val[i],2);
  }

  return sqrt(coeff*jk_err);    
  
}

double jackknifeBin(double *xsq, double *x4, double val, int block, int data_points) {
  
  //Initialise
  int num_resamp = data_points/block;
  double coeff = (1.0*num_resamp - 1.0)/(1.0*num_resamp);
  double resamp_norm = 1.0/(data_points - block);  
  double *resamp_x4 = (double*)malloc(num_resamp*sizeof(double));
  double *resamp_xsq = (double*)malloc(num_resamp*sizeof(double));
  double *resamp_val = (double*)malloc(num_resamp*sizeof(double));
  double jk_err = 0.0;
  for(int a=0; a<num_resamp; a++) {
    resamp_x4[a] = 0.0;
    resamp_xsq[a] = 0.0;
    resamp_val[a]  = 0.0;
  }

  //Get resampled averages (do normalising in calculation of Binder C.)
  for(int i=0; i<num_resamp; i++) {
    for(int j=0; j<data_points; j++) {   
      if(j < i*block || j >= (i+1)*block)  {
	resamp_x4[i] += x4[j];
	resamp_xsq[i] += xsq[j];
      }
    }
  }
  
  //Use resampled averages to get resampled values of 1-<E(x)^4>/(3*<E(x)^2>^2)
  for(int i=0; i<num_resamp; i++) {   
    resamp_val[i] = (1.0 - resamp_x4[i]/(3*resamp_xsq[i]*resamp_xsq[i]*resamp_norm));
    //cout<<(1 - resamp_x4[i]/(3*resamp_xsq[i]*resamp_xsq[i]*resamp_norm))<<endl;
  }
  
  //Get jk_error
  for(int i=0; i<num_resamp; i++) {
    jk_err += pow(val - resamp_val[i],2);
  }
  
  return sqrt(coeff*jk_err);    
  
}


void jackknifeFT(double **ind_ft_corr, double *run_ft_corr, double *jk_err,
		 int block, int data_points, int l, Param p) {
  
  //Initialise
  int Lt = p.Lt;
  int num_resamp = data_points/block;
  double coeff = (1.0*num_resamp - 1.0)/(1.0*num_resamp);
  double resamp_norm = 1.0/(data_points - block);  
  double *resamp_ft_corr = (double*)malloc(num_resamp*sizeof(double));
  for(int a=0; a<num_resamp; a++) resamp_ft_corr[a] = 0.0;
  
  //Get resampled averages
  for(int dt=0; dt<Lt/2 + 1; dt++) {
    for(int i=0; i<num_resamp; i++) {
      for(int j=0; j<data_points; j++) {   
	if(j < i*block || j >= (i+1)*block)  {
	  resamp_ft_corr[i] += ind_ft_corr[3*dt + l][j];
	}
      }
      resamp_ft_corr[i] *= resamp_norm;
    }
    
    //Get jk error
    for(int i=0; i<num_resamp; i++) {
      jk_err[dt] += pow(run_ft_corr[3*dt+l]/data_points - resamp_ft_corr[i],2);
    }
    jk_err[dt] = sqrt(coeff*jk_err[dt]);      
  }
}
