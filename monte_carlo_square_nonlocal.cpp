#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>
#include </usr/local/Cellar/gcc/4.9.2_1/lib/gcc/4.9/gcc/x86_64-apple-darwin14.3.0/4.9.2/include/omp.h>

#include "util.h"
#include "hyp_util.h"
#include "data_proc.h"
#include "data_io.h"
#include "monte_carlo_square_nonlocal.h"

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

//Declare variables to implement Wolff algorithm
int sqnl_wc_ave = 0;
int sqnl_wc_size = 0;
int sqnl_wc_calls = 0;
int sqnl_wc_poss = 0;
int sqnl_wc_t_size = 0;
int sqnl_wc_s_size = 0;

double actionSqNL(double *phi_arr, int *s, Param p,
		  double *LR_couplings, 
		  double & KE, double & PE) {  
  
  KE = 0.0;
  PE = 0.0;
  double phi_sq;
  double phi;
  double lambda_p = 0.25*p.lambda;
  double musqr_p  = 0.50*p.musqr;

  for (int i = 0; i < p.surfaceVol; i++)
    if (s[i] * phi_arr[i] < 0)
      printf("ERROR s and phi NOT aligned (actionPhi Square) ! \n");
  
  for (int i = 0; i < p.surfaceVol; i++) {    
    phi = phi_arr[i];
    phi_sq = phi*phi;

    //PE terms
    PE += lambda_p * phi_sq*phi_sq;
    PE += musqr_p  * phi_sq;

    //KE terms
#pragma omp parallel for
    for(int j=0; j<p.surfaceVol; j++) {
      int id = omp_get_thread_num();
      printf("Greetings from process %d!/n", id);
      
      //Here we are overcounting by a factor of two, hence the 0.25
      //coefficient.
      //FIXME
      KE += 0.25*(phi - phi_arr[j])*(phi - phi_arr[j])*LR_couplings[i+j*p.surfaceVol]*LR_couplings[i+j*p.surfaceVol];
    }
  }
  return PE + KE;
}

int sqnl_accept = 0;
int sqnl_tries  = 0;

int metropolisUpdateSqNL(double *phi_arr, int *s, Param &p,
			 double *LR_couplings,
			 double &delta_mag_phi, int iter) {

  delta_mag_phi = 0.0;
  
  int s_old     = 0;
  int delta_mag = 0;
  
  double phi_new = 0.0;
  double phi_new_sq = 0.0;
  double phi = 0.0;
  double phi_sq = 0.0;
  double lambda_p = 0.25*p.lambda;
  double musqr_p  = 0.50*p.musqr;
  
  double DeltaE = 0.0;
  
  for (int i = 0; i < p.surfaceVol; i++) {
    
    //Set some values we use a lot
    phi = phi_arr[i];
    DeltaE = 0.0;
    phi_new = phi + p.delta_phi * (2.0*unif(rng) - 1.0);
    phi_new_sq = phi_new*phi_new;
    phi_sq = phi*phi;
    
    if (s[i] * phi < 0) {
      printf("ERROR s and phi NOT aligned! (MUP)\n");
      exit(0);
    }
    
    //PE
    DeltaE += lambda_p*(phi_new_sq*phi_new_sq - phi_sq*phi_sq);
    DeltaE += musqr_p *(phi_new_sq            - phi_sq);
    
    //KE
    double pmpn = phi-phi_new;
    double pnsmps = 0.5*phi_new_sq-phi_sq;
#pragma omp parallel for
    for(int j=0; j<p.surfaceVol; j++) {
      DeltaE += (pmpn*phi_arr[j] + pnsmps)*LR_couplings[i+j*p.surfaceVol]*LR_couplings[i+j*p.surfaceVol];
    }
    
    sqnl_tries++;
    
    if(DeltaE < 0.0) {
      //  cout<< " Acepted  " << endl;
      s_old = s[i];
      delta_mag_phi += phi_new - phi_arr[i];
      phi_arr[i] = phi_new;
      sqnl_accept += 1;
      s[i] = (phi_new > 0) ? 1 : -1;
      delta_mag += s[i] - s_old;
    }
    else if ( unif(rng)  < exp(-DeltaE)) {
      //  cout<< " Acepted  " << endl;
      s_old = s[i];
      delta_mag_phi += phi_new - phi_arr[i];
      phi_arr[i] = phi_new;
      sqnl_accept += 1;
      s[i] = (phi_new > 0) ? 1 : -1;
      delta_mag += s[i] - s_old;
    }     
  }// end loop over lattice volume 
  
  // TUNING ACCEPTANCE 
  if (iter < p.n_therm/2 && (iter+1) % p.n_skip/10 == 0) {
    if ((double) sqnl_accept / (double) sqnl_tries < 0.5) {
      p.delta_phi -= 0.001;
    } else {
      p.delta_phi += 0.001;
    }
    if(p.n_wolff*1.0*sqnl_wc_ave/sqnl_wc_calls < p.surfaceVol && iter > p.n_skip) {
      p.n_wolff++;
    } else {
      p.n_wolff--;
      if(p.n_wolff < 3) p.n_wolff++;
    }
  }
  
  if( iter < p.n_therm ) {
    cout<<"At iter "<<iter<<" the Acceptance rate is "<<(double)sqnl_accept/(double)sqnl_tries<<endl;
    cout<<"and delta_phi is "<<p.delta_phi<<endl;
  }
  return delta_mag;
}


void wolffUpdateSqNL(double *phi_arr, int *s, Param p,
		     double *LR_couplings,
		     double &delta_mag_phi, int iter) {
  
  sqnl_wc_calls++;
  sqnl_wc_poss = 0;
  
  //Tracks which sites are acutally in the cluster.
  bool *Acluster = new bool[p.surfaceVol];
#pragma omp parallel for
  for (int i = 0; i < p.surfaceVol; i++)
    Acluster[i] = false;
  
  //Tracks which sites are potentially in the cluster.
  bool *Pcluster = new bool[p.surfaceVol];
#pragma omp parallel for
  for (int i = 0; i < p.surfaceVol; i++)
    Pcluster[i] = false;

  //Records which sites are potentially in the cluster.
  //This will have a modifiable size, hence the vector
  //style.
  vector<int> Rcluster;
  
  //Choose a random spin and identfy the possible cluster
  //with a flood fill.
  int i = int(unif(rng) * p.surfaceVol);
  int cSpin = s[i];
  clusterPossibleSqNL(i, s, cSpin, Pcluster, Rcluster, p);

  //cout<<"Maximum cluster size = "<<sqnl_wc_poss<<endl;
  
  //This function is recursive and will call itself
  //until all attempts to increase the cluster size
  //have failed.
  sqnl_wc_size = 1;
  clusterAddSqNL(i, s, cSpin, Acluster, Rcluster, LR_couplings, phi_arr, p);

  sqnl_wc_ave += sqnl_wc_size;

  if( iter%p.n_skip == 0 && iter < p.n_therm) {
    setprecision(4);
    cout<<"Using "<<p.n_wolff<<" Wolff hits."<<endl; 
    cout<<"Ave. cluster size at iter "<<iter<<" = "<<sqnl_wc_ave<<"/"<<sqnl_wc_calls<<" = "<<1.0*sqnl_wc_ave/sqnl_wc_calls<<endl;
    cout<<"S/T cluster growth ratio at iter "<<iter<<" = "<<1.0*sqnl_wc_s_size/sqnl_wc_ave<<":"<<1.0*sqnl_wc_t_size/sqnl_wc_ave<<endl;
  }
  
  delete[] Acluster;
  delete[] Pcluster;
}

void clusterPossibleSqNL(int i, int *s, int cSpin,
			 bool *Pcluster, vector<int> &Rcluster,
			 Param p) {
  //The site possibily belongs to the cluster...
  Pcluster[i] = true;
  //...so record it
  Rcluster.push_back(i);
  sqnl_wc_poss++;
  
  //If the neighbor's spin matches the cluster spin (cSpin)
  //then it is a possible cluster candidate. This is all we
  //need to identify at this point. The routine checks that
  //the candidate site is both not already identified, and
  //that it has the correct spin.

  //Forward in T
  if(!Pcluster[tp(i,p)] && s[tp(i,p)] == cSpin) {
    //cout<<"->tp";
    clusterPossibleSqNL(tp(i,p), s, cSpin, Pcluster, Rcluster, p);
  }
  
  //Forward in X
  if(!Pcluster[xp(i,p)] && s[xp(i,p)] == cSpin) {
    //cout<<"->xp";
    clusterPossibleSqNL(xp(i,p), s, cSpin, Pcluster, Rcluster, p);
  }
  
  //Backward in T 
  if(!Pcluster[ttm(i,p)] && s[ttm(i,p)] == cSpin) {  
    //cout<<"->tm";
    clusterPossibleSqNL(ttm(i,p), s, cSpin, Pcluster, Rcluster, p);
  }
  
  //Backward in X
  if(!Pcluster[xm(i,p)] && s[xm(i,p)] == cSpin) {
    //cout<<"->xm";
    clusterPossibleSqNL(xm(i,p), s, cSpin, Pcluster, Rcluster, p);
  }
}

void clusterAddSqNL(int i, int *s, int cSpin, bool *Acluster,
		    vector<int> &Rcluster, double *LR_couplings,
		    double *phi_arr, Param p) {

  // The site belongs to the cluster, so flip it.
  Acluster[i] = true;
  s[i] *= -1;
  phi_arr[i] *= -1;
  
  // ferromagnetic coupling
  const double J = 1.0;

  double prob = 0.0;
  double rand = 0.0;

  int idx = 0;
  
  //We now loop over the possible lattice sites, adding sites
  //(creating bonds) with the specified LR probablity.
  for(int j=0; j<sqnl_wc_poss; j++) {
    idx = Rcluster[j];
    if(LR_couplings[i + idx*p.surfaceVol] > 1e-7){
      //cout<<Rcluster[j]<<endl;
      prob = 1 - exp(2*J*phi_arr[i]*phi_arr[idx]*LR_couplings[i + idx*p.surfaceVol]);
      rand = unif(rng);
      //cout<<"Testing "<<i<<","<<Rcluster[j]<<" "<<rand<<" "<<prob<<" ";
      if(rand < prob) {
	//cout<<"hit"<<endl;
	sqnl_wc_size++;
	//next = Rcluster[j];
	//Rcluster.erase(Rcluster.begin() + j);
	clusterAddSqNL(idx, s, cSpin, Acluster, Rcluster,
		       LR_couplings, phi_arr, p);
      }    
      //cout<<"miss"<<endl;
    }
  }
  //cout<<endl<<endl<<"New Base"<<endl;
}

//------------- Monte Carlo Update  ----------------//
void runMonteCarloSqNL(vector<Vertex> &NodeList, Param p) {
  
  double KE = 0.0, PE = 0.0;
  double mag_phi = 0.0;
  double delta_mag_phi = 0.0;
  
  //Create the 1/r^{(2+\sigma)} LR couplings;
  double sigma = p.sigma;
  double *LR_couplings = new double[p.surfaceVol*p.surfaceVol];
  double r_x = 0.0;
  double r_y = 0.0;

  for(int i=0; i<p.surfaceVol; i++) {
#pragma omp parallel for
    for(int j=0; j<p.surfaceVol; j++) {
      
      r_x = abs(i%p.S1 - j%p.S1);
      if(r_x > p.S1/2) r_x = p.S1 - r_x;
      
      r_y = abs(i/p.Lt - j/p.Lt);
      if(r_y > p.Lt/2) r_y = p.Lt - r_y;
      
      if(i != j) LR_couplings[i + j*p.surfaceVol] = pow(sqrt(r_x*r_x + r_y*r_y), -(2+sigma));
      //Quick symmetry check...
      //if( (i == 10 && j == 20) || (i==20 && j == 10)) cout<<i<<" "<<j<<" "<<LR_couplings[i + j*p.surfaceVol]<<endl;
    }
  }
  
  int *s = (int*)malloc(p.surfaceVol*sizeof(int));
  double *phi = (double*)malloc(p.surfaceVol*sizeof(double));
  for(int i = 0;i < p.surfaceVol; i++) {
    phi[i] = 2.0*unif(rng) - 1.0;
    s[i] = (phi[i] > 0) ? 1:-1;
    mag_phi += phi[i];
  }  
  thermaliseSqNL(phi, s, p, delta_mag_phi, LR_couplings);
  
  //Arrays holding measurements for error analysis
  double E_arr[p.n_meas];
  double E2_arr[p.n_meas];
  double PhiAb_arr[p.n_meas];
  double Phi_arr[p.n_meas];
  double Phi2_arr[p.n_meas];
  double Phi4_arr[p.n_meas];
  for(int i=0; i<p.n_meas; i++) {
    E_arr[i]     = 0.0;
    E2_arr[i]    = 0.0;
    PhiAb_arr[i] = 0.0;
    Phi_arr[i]   = 0.0;
    Phi2_arr[i]  = 0.0;
    Phi4_arr[i]  = 0.0;
  }

  //Running averages
  double tmpE     = 0.0;  
  double aveE     = 0.0;
  double aveKE    = 0.0;
  double avePE    = 0.0;
  double aveE2    = 0.0;
  double avePhiAb = 0.0;
  double avePhi   = 0.0;
  double avePhi2  = 0.0;
  double avePhi4  = 0.0;
  double MagPhi   = 0.0;

  //Density of lattice sites
  double rhoVol  = 1.0/(double)p.surfaceVol;
  double rhoVol2 = rhoVol*rhoVol;

  //correlation function arrays. One keeps the running average,
  //one is temporary to pass around single measurent values. 
  double **corr_tmp = (double**)malloc((p.S1/2)*sizeof(double*));
  double **corr_ave = (double**)malloc((p.S1/2)*sizeof(double*));
  for(int i=0; i<p.S1/2; i++) {
    corr_tmp[i] = (double*)malloc((p.Lt/2)*sizeof(double));
    corr_ave[i] = (double*)malloc((p.Lt/2)*sizeof(double));
  }
  for(int i=0; i<p.S1/2; i++)
    for(int j=0; j<p.Lt/2; j++) 
      corr_ave[i][j] = 0.0;
  
  double **phi_sq_arr = (double**)malloc(p.S1*sizeof(double*));
  for(int i=0; i<p.S1; i++) {
    phi_sq_arr[i] = (double*)malloc(p.Lt*sizeof(double));
    for(int j=0; j<p.Lt; j++) phi_sq_arr[i][j] = 0.0;
  }

  int idx = 0;
  double norm;
  
  for(int iter = p.n_therm; iter < p.n_therm + p.n_skip*p.n_meas; iter++) {
    
    for(int i=0; i<p.n_wolff; i++) {
      wolffUpdateSqNL(phi, s, p, LR_couplings, delta_mag_phi, iter);      
    }
    metropolisUpdateSqNL(phi, s, p, LR_couplings, delta_mag_phi, iter);
    
    //Take measurements.
    if((iter+1) % p.n_skip == 0) {
      
      for(int i = 0;i < p.S1; i++) {
	for(int j=0; j<p.Lt; j++) 
	  phi_sq_arr[i][j] += pow(phi[i + p.S1*j],2);
      }
      
      tmpE   = actionSqNL(phi, s, p, LR_couplings, KE, PE);
      aveKE += rhoVol*KE;
      avePE += rhoVol*PE;
      aveE  += rhoVol*tmpE;
      aveE2 += rhoVol2*tmpE*tmpE;

      MagPhi = 0.0;
#pragma omp parallel for
      for(int i = 0;i < p.surfaceVol; i++) MagPhi += phi[i];
      MagPhi *= rhoVol;
      
      avePhiAb  += abs(MagPhi);
      avePhi    += MagPhi;
      avePhi2   += MagPhi*MagPhi;
      avePhi4   += MagPhi*MagPhi*MagPhi*MagPhi;
      
      E_arr[idx]     = rhoVol*tmpE;
      E2_arr[idx]    = rhoVol*tmpE*tmpE;
      PhiAb_arr[idx] = abs(MagPhi);
      Phi_arr[idx]   = MagPhi;
      Phi2_arr[idx]  = MagPhi*MagPhi;
      Phi4_arr[idx]  = MagPhi*MagPhi*MagPhi*MagPhi;
      
      idx++;
      
      cout<<setprecision(8);
      norm = 1.0/(idx);

      //Dump to stdout
      cout<<"Measurement "<<(iter+1 - p.n_therm)/p.n_skip<<" Sweep "<<iter+1<<endl;
      cout<<"Ave Energy= "<<aveE*norm<<endl;
      cout<<"Ave KE    = "<<aveKE*norm<<endl;
      cout<<"Ave PE    = "<<avePE*norm<<endl;
      cout<<"Ave |phi| = "<<avePhiAb*norm<<endl;
      cout<<"Ave phi   = "<<avePhi*norm<<endl;
      cout<<"Ave phi^2 = "<<avePhi2*norm<<endl;
      cout<<"Ave phi^4 = "<<avePhi4*norm<<endl;
      cout<<"Suscep    = "<<(avePhi2*norm-pow(avePhiAb*norm,2))/rhoVol<<endl;
      cout<<"Spec Heat = "<<(aveE2*norm-pow(aveE*norm,2))/rhoVol<<endl;
      cout<<"Binder    = "<<1.0-avePhi4/(3.0*avePhi2*avePhi2*norm)<<" ("<<norm*avePhi2*avePhi2/avePhi4<<")"<<endl;
      
      //Visualisation tools
      visualiserSqr(phi, avePhiAb*norm, p);
      //visualiserPhi2(phi_sq_arr, p, idx);      

      //Calculate correlaton functions and update the average.
      correlators(corr_tmp, corr_ave, idx, phi, avePhi*norm, p);
      
      ofstream filet("correlators_t.dat");
      for(int i=0; i<p.Lt/2; i++) {
	filet << i << " " << corr_tmp[0][i] << endl;
      }
      filet.close();
      
      ofstream files("correlators_s.dat");
      for(int i=0; i<p.S1/2; i++) {
	files << i << " " << corr_tmp[i][0] << endl;
      }
      files.close();
      //if(idx%50 == 0) corr_eigs(corr_run, p);
      
    }
  }

  //autocorrelation(PhiAb_arr, avePhiAb, p.n_meas);
  
  //correlators(corr_run, corr_ave, idx, phi_cyl, avePhi*norm, p);
  //corr_eigs(corr_run, p);
  
  //free(corr_run);
  //free(corr_ave);

  free(s);
  free(phi);
}

void thermaliseSqNL(double *phi, int *s, Param p,
		    double &delta_mag_phi, double *LR_couplings) {
  
  for(int iter = 0; iter < p.n_therm; iter++) {
    for(int i=0; i<p.n_wolff; i++) {
      wolffUpdateSqNL(phi, s, p, LR_couplings, delta_mag_phi, iter);      
    }
    metropolisUpdateSqNL(phi, s, p, LR_couplings, delta_mag_phi, iter);
    if((iter+1)%p.n_skip == 0) cout<<"Therm sweep "<<iter+1<<endl;
  }
}
