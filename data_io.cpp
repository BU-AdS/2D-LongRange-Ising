#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>

#include "util.h"
#include "data_proc.h"

using namespace std;

// XXX Does idx represent the number of samples before the current one?
void writeObservables(double **ind_corr, double *run_corr, int *norm_corr,
		      double **ind_ft_corr, double *run_ft_corr,
		      int idx, observables obs, Param p){
 
  double norm = 1.0/idx;
  int S1 = p.S1;
  int Lt = p.Lt;
  int vol = S1*Lt;
  int x_len = S1/2 + 1;
  char fname[256];

  //The observables (spec heat, suscep, etc) are used only as a guide
  //to discern criticality, The real quantities of interest are
  //The critical exponents, especally \eta, and its dependence
  //on \sigma, hence MUST have proper error estimates of the
  //correlation function values.

  for(int dth = 0; dth <S1/2+1; dth++) {

    double *jk_err = (double*)malloc((Lt/2+1)*sizeof(double));
    jackknife(ind_corr, run_corr, jk_err, p.n_jkblock, idx, Lt/2+1, dth, p);
    
    sprintf(fname, "correlators_dth%d.dat", dth);    
    FILE *fp = fopen(fname, "w"); 
    for(int i=0; i<Lt/2+1; i++) {
      
      fprintf(fp, "%d %.15e %.15e\n", i, run_corr[dth + i*x_len]*norm/(norm_corr[dth + i*x_len]),
	      jk_err[i]*norm/(norm_corr[dth + i*x_len]));
      
    }
    fclose(fp);
    free(jk_err);
  }

  for(int dt = 0; dt <Lt/2+1; dt++) {

    double *jk_err = (double*)malloc((S1/2+1)*sizeof(double));
    jackknife(ind_corr, run_corr, jk_err, p.n_jkblock, idx, S1/2+1, dt, p);
    
    sprintf(fname, "correlators_dt%d.dat", dt);
    FILE *fp = fopen(fname, "w"); 
    ofstream file(fname);
    for(int i=0; i<S1/2+1; i++) {

      fprintf(fp, "%d %.15e %.15e\n", i, run_corr[i + dt*x_len]*norm/(norm_corr[i + dt*x_len]),
	      jk_err[i]*norm/(norm_corr[i + dt*x_len]));
      
    }
    fclose(fp);
    free(jk_err);
  }

  //Go to l=2
  for(int l=0; l<3; l++) {
    
    sprintf(fname, "correlators_FTl%d.dat", l);
    FILE *fp = fopen(fname, "w"); 
    double *jk_err = (double*)malloc((Lt/2+1)*sizeof(double));
    jackknifeFT(ind_ft_corr, run_ft_corr, jk_err, p.n_jkblock, idx, l, p);
    
    for(int dt=0; dt<Lt/2+1; dt++) {

      fprintf(fp, "%d %.15e %.15e\n", dt, run_ft_corr[3*dt + l]*norm, jk_err[dt]);
      
    }
    fclose(fp);
    free(jk_err);
  }
  
  
  double J = p.J;
  if(p.theory_type != ISING) J = 1.0;
  
  //Jacknife the second moments of means.
  double jkErrSuscep;
  jkErrSuscep = (vol*J)*jackknifeVar(obs.PhiAb_arr, obs.Phi2_arr,
				     obs.Suscep[idx-1]/(vol*J), p.n_jkblock, idx, J);
  
  double jkErrSpecHeat;
  jkErrSpecHeat = J*J*jackknifeVar(obs.E_arr, obs.E2_arr,
				   obs.SpecHeat[idx-1]/(J*J), p.n_jkblock, idx, J*J);
  
  //Jacknife the Binder cumumlant
  double jkErrBinder;
  jkErrBinder = jackknifeBin(obs.Phi2_arr, obs.Phi4_arr,
			     obs.Binder[idx-1], p.n_jkblock, idx);

  //Dump jackknifed observables
  sprintf(fname, "JKobservables.dat");
  FILE *fp = fopen(fname, "a"); 

  // Normalize Monte Carlo average by number of samples
  // XXX Check
  norm = 1.0/idx;

  fprintf(fp, "%d %.15e ", idx, p.sigma);
  if(p.theory_type == ISING) fprintf(fp, "%.15e ", p.J);
  if(p.theory_type == PHI4) fprintf(fp, "%.15e %.15e ", p.musqr, p.lambda);
  fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.15e\n",
	  obs.SpecHeat[idx-1], jkErrSpecHeat,
	  obs.Suscep[idx-1], jkErrSuscep, 
// Include estimate for the average magnetization XXX Check
      obs.avePhi*norm,
	  obs.Binder[idx-1], jkErrBinder);
  fclose(fp);

  /*
  //Dump raw observable data for any further offline analysis
  ofstream file_obs;
  file_obs.open("observables.dat", ios_base::app);
  int meas = 0;
  for(int i=0; i<10; i++) {

    std::cout.precision(15);
    std::cout<<std::scientific;
    
    meas = idx - 10 + i;
    file_obs<<meas<<" "<<p.sigma<<" "<<p.J<<" "<<obs.E_arr[meas]<<" ";
    file_obs<<obs.E2_arr[meas]<<" "<<obs.PhiAb_arr[meas]<<" "<<obs.Phi_arr[meas]<<" ";
    file_obs<<obs.Phi2_arr[meas]<<" "<<obs.Phi4_arr[meas]<<" "<<obs.Suscep[meas]<<" ";
    file_obs<<obs.SpecHeat[meas]<<" "<<obs.Binder[meas]<<endl;
  }
  file_obs.close();  
  */
}

void writePhi3(double **ind_corr_phi_phi3, double *run_corr_phi_phi3,
	       double **ind_corr_phi3_phi3, double *run_corr_phi3_phi3,
	       double **ind_ft_corr_phi_phi3, double *run_ft_corr_phi_phi3,
	       double **ind_ft_corr_phi3_phi3, double *run_ft_corr_phi3_phi3,
	       int *norm_corr, int idx, observables obs, Param p){
  
  double norm = 1.0/idx;
  int S1 = p.S1;
  int Lt = p.Lt;
  int x_len = S1/2 + 1;
  char fname[256];

  for(int dth = 0; dth <S1/2+1; dth++) {

    double *jk_err = (double*)malloc((Lt/2+1)*sizeof(double));
    
    //Phi Phi3
    jackknife(ind_corr_phi_phi3, run_corr_phi_phi3, jk_err, p.n_jkblock, idx, Lt/2+1, dth, p);
    sprintf(fname, "correlatorsPhiPhi3_dth%d.dat", dth);
    FILE *fp = fopen(fname, "a");
    for(int i=0; i<Lt/2+1; i++) {
      fprintf(fp, "%d %.15e %.15e\n", i, run_corr_phi_phi3[dth + i*x_len]*norm/(norm_corr[dth + i*x_len]),
	      jk_err[i]*norm/(norm_corr[dth + i*x_len]) );
    }
    fclose(fp);
    
    jackknife(ind_corr_phi3_phi3, run_corr_phi3_phi3, jk_err, p.n_jkblock, idx, Lt/2+1, dth, p);    
    sprintf(fname, "correlatorsPhi3Phi3_dth%d.dat", dth);
    fp = fopen(fname, "a");    
    for(int i=0; i<Lt/2+1; i++) {
      fprintf(fp, "%d %.15e %.15e\n", i, run_corr_phi3_phi3[dth + i*x_len]*norm/(norm_corr[dth + i*x_len]),
	      jk_err[i]*norm/(norm_corr[dth + i*x_len]) );      
    }
    fclose(fp);
    free(jk_err);
  }

  for(int dt = 0; dt <Lt/2+1; dt++) {

    double *jk_err = (double*)malloc((S1/2+1)*sizeof(double));

    
    jackknife(ind_corr_phi_phi3, run_corr_phi_phi3, jk_err, p.n_jkblock, idx, S1/2+1, dt, p);    
    sprintf(fname, "correlatorsPhiPhi3_dt%d.dat", dt);
    FILE *fp = fopen(fname, "a");
    for(int i=0; i<S1/2+1; i++) {
      fprintf(fp, "%d %.15e %.15e\n", i, run_corr_phi_phi3[i + dt*x_len]*norm/(norm_corr[i + dt*x_len]),
	      jk_err[i]*norm/(norm_corr[i + dt*x_len]) );      
    }
    fclose(fp);

    jackknife(ind_corr_phi3_phi3, run_corr_phi3_phi3, jk_err, p.n_jkblock, idx, S1/2+1, dt, p);
    sprintf(fname, "correlatorsPhi3Phi3_dt%d.dat", dt);
    fp = fopen(fname, "a");
    for(int i=0; i<S1/2+1; i++) {
      fprintf(fp, "%d %.15e %.15e\n", i, run_corr_phi3_phi3[i + dt*x_len]*norm/(norm_corr[i + dt*x_len]),
	      jk_err[i]*norm/(norm_corr[i + dt*x_len]) );      
    }
    fclose(fp);
    free(jk_err);    
  }

  //Go to l=2
  for(int l=0; l<3; l++) {

    double *jk_err = (double*)malloc((Lt/2+1)*sizeof(double));
    ofstream file;
    
    jackknifeFT(ind_ft_corr_phi_phi3, run_ft_corr_phi_phi3, jk_err, p.n_jkblock, idx, l, p);
    sprintf(fname, "correlatorsPhiPhi3_FTl%d.dat", l);
    FILE *fp = fopen(fname, "a");
    for(int dt=0; dt<Lt/2+1; dt++) {
      fprintf(fp, "%d %.15e %.15e\n", dt, run_ft_corr_phi_phi3[3*dt + l]*norm, jk_err[dt]);
    }
    fclose(fp);
    
    jackknifeFT(ind_ft_corr_phi3_phi3, run_ft_corr_phi3_phi3, jk_err, p.n_jkblock, idx, l, p);
    sprintf(fname, "correlatorsPhi3Phi3_FTl%d.dat", l);
    fp = fopen(fname, "a");
    for(int dt=0; dt<Lt/2+1; dt++) {
      fprintf(fp, "%d %.15e %.15e\n", dt, run_ft_corr_phi3_phi3[3*dt + l]*norm, jk_err[dt]);
    }
    fclose(fp);
    free(jk_err);
  }
}


void visualiserIsing(int *s, Param p) {  
  
  for(int i=0; i<p.S1; i++) {
    for(int j=0; j<p.Lt; j++) {
      if(s[i + p.S1*j] < 0) cout<<"\033[1;41m \033[0m";
      else cout<<"\033[1;44m \033[0m";
    }
    cout<<endl;
  }
}

void visualiserIsingCluster(int *s, bool *cluster, Param p) {  
  
  for(int i=0; i<p.S1; i++) {
    for(int j=0; j<p.Lt; j++) {
      if(cluster[i +p.S1*j]) cout<<"\033[1;42m \033[0m";
      else if(s[i + p.S1*j] < 0) cout<<"\033[1;41m \033[0m";
      else cout<<"\033[1;44m \033[0m";
      
    }
    cout<<endl;
  }
}

void visualiserSqr(double *phi_cyl, double barr, Param p) {  
  
  for(int i=0; i<p.S1; i++) {
    for(int j=0; j<p.Lt; j++) {
      if(phi_cyl[i + p.S1*j] < -1.5*barr) cout<<"\033[1;41m \033[0m";
      if(-1.5*barr < phi_cyl[i + p.S1*j] && phi_cyl[i + p.S1*j] < -0.5*barr) cout<<"\033[1;43m \033[0m";
      if(-0.5*barr < phi_cyl[i + p.S1*j] && phi_cyl[i + p.S1*j] <  0.5*barr) cout<<"\033[1;42m \033[0m";
      if( 0.5*barr < phi_cyl[i + p.S1*j] && phi_cyl[i + p.S1*j] <  1.5*barr) cout<<"\033[1;46m \033[0m";
      if( 1.5*barr < phi_cyl[i + p.S1*j]) cout<<"\033[1;44m \033[0m";
    }
    cout<<endl;
  }
  //usleep(250000);
}


void visualiserPhi2(double **phi_sq, Param p, int iter) {
  
  double barr_t[p.S1];
  double barr = 0.0;
  
  for(int i=0; i<p.S1; i++) {
    barr_t[i] = 0.0;
    for(int j=0; j<p.Lt; j++) {
      barr_t[i] += phi_sq[i][j];
    }
    barr_t[i] /= p.Lt;
    barr      += barr_t[i];
  }
  barr /= p.S1;
  
  for(int i=0; i<p.S1; i++) {
    for(int j=0; j<p.Lt; j++) {
      if( phi_sq[i][j] < barr ) cout<<"\033[1;41m \033[0m";
      if( barr < phi_sq[i][j] ) cout<<"\033[1;44m \033[0m";
    }
    cout<<" phi_sq_ave[s="<<i<<"] = "<<barr_t[i]/(iter+1)<<endl;
  }
}
