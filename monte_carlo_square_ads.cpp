#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>

#include "monte_carlo_square_ads.h"
#include "util.h"
#include "data_proc.h"
#include "data_io.h"

using namespace std;

extern int seed;
extern mt19937_64 rng;
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

// declare variables to implement Wolff algorithm
int sqads_wc_ave = 0;
int sqads_wc_size = 0;
int sqads_wc_poss = 0;
int sqads_wc_calls = 0;
int sqads_wc_t_size = 0;
int sqads_wc_s_size = 0;

double actionSqAdS(double *phi_arr, int *s, Param p,
		 double *LR_couplings, 
		 double &KE,  double &PE) {
  
  KE = 0.0;
  PE = 0.0;
  double phi_sq;
  double phi;
  double lambda_p = 0.25*p.lambda;
  double musqr_p  = 0.50*p.musqr;
  
  for (int i = 0; i < p.surfaceVol; i++)
    if (s[i] * phi_arr[i] < 0)
      printf("ERROR s and phi NOT aligned ! (actionPhi AdS)\n");
  

  for (int i = 0; i < p.surfaceVol; i++) {    
    phi = phi_arr[i];
    phi_sq = phi*phi;

    //PE terms
    PE += lambda_p * phi_sq*phi_sq;
    PE += musqr_p  * phi_sq;
    
    //KE terms
    //Here we use (a funtion of) the AdS Propagator to construct
    //the LR coupling. This requires us to know the explict value
    //of the AdS Propagator from i to j. We could exploit translational
    //symmetry in the temporal diretcion, and the symmetry in the qth
    //segment of the AdS disk, but this would require extra cycles to
    //make sense of the indices, hence the lookup table is large and
    //repetitive.
    
    for(int j=0; j<p.surfaceVol; j++) {
      //Here we are overcounting by a factor of two, hence the 0.25
      //coefficient.
      //FIXME
      KE += 0.25*(phi - phi_arr[j])*(phi - phi_arr[j])*LR_couplings[i+j*p.surfaceVol]*LR_couplings[i+j*p.surfaceVol];
    }
  }  
  return PE + KE;
}

int sqads_accept = 0;
int sqads_tries  = 0;

int metropolisUpdateSqAdS(double *phi_arr, int *s, Param &p,
			double *LR_couplings, 
			double &delta_mag_phi, int iter) {
  
  delta_mag_phi = 0.0;
  uniform_real_distribution<double> unif(0.0,1.0);

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
      cout<<"ERROR s and phi NOT aligned! (MetroUpdate AdS) ";
      cout<<"S["<<i<<"] = "<<s[i]<<" and phi["<<i<<"] = "<<phi<<endl;
    }
    
    //PE
    DeltaE += lambda_p*(phi_new_sq*phi_new_sq - phi_sq*phi_sq);
    DeltaE += musqr_p *(phi_new_sq            - phi_sq);
    
    //KE
    double pmpn = phi-phi_new;
    double pnsmps = 0.5*phi_new_sq-phi_sq;
    for(int j=0; j<p.surfaceVol; j++) {
      DeltaE += (pmpn*phi_arr[j] + pnsmps)*LR_couplings[i+j*p.surfaceVol]*LR_couplings[i+j*p.surfaceVol];
      //cout<<LR_couplings[i+j*p.surfaceVol]<<endl;
    }
    
    sqads_tries += 1;    

    if(DeltaE < 0.0) {
      //  cout<< " Acepted  " << endl;
      s_old = s[i];
      delta_mag_phi += phi_new - phi;
      phi_arr[i] = phi_new;
      sqads_accept += 1;
      s[i] = (phi_new > 0) ? 1 : -1;
      delta_mag += s[i] - s_old;
    }
    else if ( unif(rng)  < exp(-DeltaE)) {
      //  cout<< " Acepted  " << endl;
      s_old = s[i];
      delta_mag_phi += phi_new - phi;
      phi_arr[i] = phi_new;
      sqads_accept += 1;
      s[i] = (phi_new > 0) ? 1 : -1;
      delta_mag += s[i] - s_old;
    }     
  }// end loop over lattice volume 
  
  // TUNING ACCEPTANCE 
  if (iter < p.n_therm/2 && (iter+1) % p.n_skip/10 == 0) {
    if ((double)sqads_accept / (double) sqads_tries < 0.5) {
      p.delta_phi -= 0.001;
    } else {
      p.delta_phi += 0.001;
    }
    if(p.n_cluster*1.0*sqads_wc_ave/sqads_wc_calls < p.surfaceVol && iter > p.n_skip) {
      p.n_cluster++;
    } else {
      p.n_cluster--;
      if(p.n_cluster < 2) p.n_cluster++;
    }
  }
  
  if( iter%p.n_skip < 4) {
    cout<<"At iter "<<iter<<" the Acceptance rate is "<<(double)sqads_accept/(double)sqads_tries<<endl;
    cout<<"and delta_phi is "<<p.delta_phi<<endl;
  }
  
  return delta_mag;
}

void wolffUpdateSqAdS(double *phi_arr, int *s, Param p,
		    double *LR_couplings,
		    double &delta_mag_phi, int iter) {
  
  sqads_wc_calls++;
  sqads_wc_poss = 0;

  //Tracks which sites are acutally in the cluster.
  bool *Acluster = new bool[p.surfaceVol];
  for (int i = 0; i < p.surfaceVol; i++)
    Acluster[i] = false;
  
  //Tracks which sites are potentially in the cluster.
  bool *Pcluster = new bool[p.surfaceVol];
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
  clusterPossibleSqAdS(i, s, cSpin, Pcluster, Rcluster, p);
  
  //cout<<"Maximum cluster size = "<<sqnl_wc_poss<<endl;
  
  //This function is recursive and will call itself
  //until all q+2 attempts in the lattice directions
  //have failed to increase the cluster.
  sqads_wc_size = 1;
  clusterAddSqAdS(i, s, cSpin, Acluster, Rcluster, LR_couplings, phi_arr, p);

  sqads_wc_ave += sqads_wc_size;

  if( iter%p.n_skip == 0 && iter < p.n_therm) {
    setprecision(4);
    cout<<"Using "<<p.n_cluster<<" Wolff hits."<<endl; 
    cout<<"Ave. cluster size at iter "<<iter<<" = "<<sqads_wc_ave<<"/"<<sqads_wc_calls<<" = "<<1.0*sqads_wc_ave/sqads_wc_calls<<endl;
    cout<<"S/T cluster growth ratio at iter "<<iter<<" = "<<1.0*sqads_wc_s_size/sqads_wc_ave<<":"<<1.0*sqads_wc_t_size/sqads_wc_ave<<endl;
  }

  delete[] Acluster;
  delete[] Pcluster;
}

void clusterPossibleSqAdS(int i, int *s, int cSpin,
			bool *Pcluster, vector<int> &Rcluster,
			Param p) {
  //The site possibily belongs to the cluster...
  Pcluster[i] = true;
  //...so record it
  Rcluster.push_back(i);
  sqads_wc_poss++;
  
  //If the neighbor's spin matches the cluster spin (cSpin)
  //then it is a possible cluster candidate. This is all we
  //need to identify at this point. The routine checks that
  //the candidate site is both not already identified, and
  //that it has the correct spin.

  //Forward in T
  if(!Pcluster[tp(i,p)] && s[tp(i,p)] == cSpin) {
    //cout<<"->tp";
    clusterPossibleSqAdS(tp(i,p), s, cSpin, Pcluster, Rcluster, p);
  }
  
  //Forward in X
  if(!Pcluster[xp(i,p)] && s[xp(i,p)] == cSpin) {
    //cout<<"->xp";
    clusterPossibleSqAdS(xp(i,p), s, cSpin, Pcluster, Rcluster, p);
  }
  
  //Backward in T 
  if(!Pcluster[ttm(i,p)] && s[ttm(i,p)] == cSpin) {  
    //cout<<"->tm";
    clusterPossibleSqAdS(ttm(i,p), s, cSpin, Pcluster, Rcluster, p);
  }
  
  //Backward in X
  if(!Pcluster[xm(i,p)] && s[xm(i,p)] == cSpin) {
    //cout<<"->xm";
    clusterPossibleSqAdS(xm(i,p), s, cSpin, Pcluster, Rcluster, p);
  }
}

void clusterAddSqAdS(int i, int *s, int cSpin, bool *Acluster,
		   vector<int> &Rcluster, double *LR_couplings,
		   double *phi_arr, Param p) {
  
  // The site belongs to the cluster, so flip it.
  Acluster[i] = true;
  s[i] *= -1;
  phi_arr[i] *= -1;
  
  // ferromagnetic coupling
  const double J = +1.0;

  double prob = 0.0;
  double rand = 0.0;

  int idx = 0;
  
  //We now loop over the possible lattice sites, adding sites
  //(creating bonds) with the specified LR probablity.
  for(int j=0; j<sqads_wc_poss; j++) {
    idx = Rcluster[j];
    if(LR_couplings[i + idx*p.surfaceVol] > 1e-7 && !Acluster[idx]){
      //cout<<Rcluster[j]<<endl;
      prob = 1 - exp(2*J*phi_arr[i]*phi_arr[idx]*LR_couplings[i + idx*p.surfaceVol]);
      rand = unif(rng);
      //cout<<"Testing "<<i<<","<<Rcluster[j]<<" "<<rand<<" "<<prob<<" ";
      if(rand < prob) {
	//cout<<"hit"<<endl;
	sqads_wc_size++;
	//next = Rcluster[j];
	//Rcluster.erase(Rcluster.begin() + j);
	clusterAddSqAdS(idx, s, cSpin, Acluster, Rcluster,
		       LR_couplings, phi_arr, p);
      }    
      //cout<<"miss"<<endl;
    }
  }
  //cout<<endl<<endl<<"New Base"<<endl;
}


//------------- Monte Carlo Update  ----------------//
void runMonteCarloSqAdS(vector<Vertex> &NodeList, Param p) {
  
  double KE = 0.0, PE = 0.0;
  double mag_phi = 0.0;
  double delta_mag_phi = 0.0;

  double *LR_couplings = new double[p.surfaceVol*p.surfaceVol];  
  LRAdSCouplings(LR_couplings, NodeList, p);

  cout<<"1"<<endl;
  
  int *s = new int[p.surfaceVol];
  double *phi = new double[p.surfaceVol];
  for(int i=0; i<p.surfaceVol; i++) {
    phi[i] = 2.0*unif(rng) - 1.0;
    s[i] = (phi[i] > 0) ? 1:-1;
    mag_phi += phi[i];
  } 
  thermaliseSqAdS(phi, s, p, delta_mag_phi, LR_couplings);
  
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
    
    if((iter+1)%p.n_cluster == 0 && p.n_cluster != 0) {
      metropolisUpdateSqAdS(phi, s, p, LR_couplings, delta_mag_phi, iter);
    } else {
      wolffUpdateSqAdS(phi, s, p, LR_couplings, delta_mag_phi, iter);
    }
    
    //Take measurements.
    if((iter+1) % p.n_skip == 0) {
      
      for(int i = 0;i < p.S1; i++) {
	for(int j=0; j<p.Lt; j++) 
	  phi_sq_arr[i][j] += pow(phi[i + p.S1*j],2);
      }
      
      tmpE   = actionSqAdS(phi, s, p, LR_couplings, KE, PE);
      aveKE += rhoVol*KE;
      avePE += rhoVol*PE;
      aveE  += rhoVol*tmpE;
      aveE2 += rhoVol2*tmpE*tmpE;

      MagPhi = 0.0;
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
      cout<<"Binder    = "<<1.0-avePhi4/(3.0*avePhi2*avePhi2*norm)<<endl;
      
      //Visualisation tools
      visualiserSqr(phi, avePhiAb*norm, p);
      //visualiserPhi2(phi_sq_arr, p, idx);      

      //Calculate correlaton functions and update the average.
      //correlators(corr_tmp, corr_ave, idx, phi, avePhi*norm, p);
      
      ofstream myfilet("correlators_t.dat");
      for(int i=0; i<p.Lt/2; i++) {
	myfilet << i << " " << corr_tmp[0][i] << endl;
      }
      myfilet.close();
      
      ofstream myfiles("correlators_s.dat");
      for(int i=0; i<p.S1/2; i++) {
	myfiles << i << " " << corr_tmp[i][0] << endl;
      }
      myfiles.close();
      //if(idx%50 == 0) corr_eigs(corr_run, p);
      
    }
  }

  //autocorrelation(PhiAb_arr, avePhiAb, p.n_meas);
  
  //correlators(corr_run, corr_ave, idx, phi_cyl, avePhi*norm, p);
  //corr_eigs(corr_run, p);
  
  //free(corr_run);
  //free(corr_ave);
  
  delete[] s;
  delete[] phi;
  delete[] LR_couplings;;  
}

void thermaliseSqAdS(double *phi, int *s, Param p,
		     double &delta_mag_phi, double *LR_couplings) {
  
  for(int iter = 0; iter < p.n_therm; iter++) {
    for(int i=0; i<p.n_cluster; i++) {
      if(p.useWolff) wolffUpdateSqAdS(phi, s, p, LR_couplings, delta_mag_phi, iter);
      else if((iter+1)%p.n_skip == 0) cout<<"Therm sweep "<<iter+1<<endl;
    }
    metropolisUpdateSqAdS(phi, s, p, LR_couplings, delta_mag_phi, iter);    
    if((iter+1)%p.n_skip == 0) cout<<"Therm sweep "<<iter+1<<endl;
  }
}

void LRAdSCouplings(double *LR_couplings, vector<Vertex> NodeList, Param &p){

  int pos = endNode(p.Levels-1, p) + 1;  
  double rad = 0.0;
  for(int i=0; i<p.S1; i++)
    rad += abs(NodeList[pos + i%p.S1 + (i/p.S1)*p.AdSVol].z)/p.S1;

  rad = 0.99;
  
  cout<<rad<<endl;
  vector<complex<long double>> circum(p.S1);
  for(int i=0; i<p.S1; i++) {
    circum[i].real(rad*cos(2*i*(M_PI/p.S1)));
    circum[i].imag(rad*sin(2*i*(M_PI/p.S1)));
  }

  for(int i=0; i<p.S1; i++) {
    double delta_t = i > p.Lt/2 ? p.Lt - i : i;
    delta_t = delta_t*(p.t_weight_scale);
    //cout<<AdS2p1Prop(circum[0], circum[i], 0, p)/AdS2p1Prop(circum[0], circum[0], i*p.t_weight_scale, p)<<endl;
    cout<<i<<" "<<pow(AdS2p1Prop(circum[0], circum[i], 0, p)/AdS2p1Prop(circum[0], circum[0], delta_t, p),p.sigma)<<endl;
  }
  exit(0);
    
  //Here we scale the LR couplings so that nn interactions are
  //approximately unit valued. We do this by raking the ratio of
  //G(i,t;i+1,t) to G(i,t;i,t+1)
  double celeritas = (AdS2p1Prop(circum[0], circum[1], 0, p) /
		      AdS2p1Prop(circum[0], circum[0], 1, p));

  double c_tol = 1e-8;
  double kappa = p.t_weight_scale;

  while(abs(1-celeritas) > c_tol) {
    celeritas = (sigmaL(circum[0], circum[1], 0) /
		 sigmaL(circum[0], circum[0], kappa));
    //celeritas = (AdS2p1Prop(circum[0], circum[1], 0, p) /
    //AdS2p1Prop(circum[0], circum[0], kappa, p));
    if(celeritas > 1) kappa += 1e-10;
    if(celeritas < 1) kappa -= 1e-10;
    setprecision(10);
    cout<<kappa<<endl;
  }
  
  p.t_weight_scale = kappa;
  
  for(int i=0; i<p.surfaceVol; i++){
    for(int j=0; j<p.surfaceVol; j++){
      //index divided by disk size, using the int floor feature/bug,
      //gives the timeslice for each index.
      int t1 = i / p.S1;
      int t2 = j / p.S1;
      double delta_t = abs(t2-t1) > p.Lt/2 ? p.Lt - abs(t2-t1) : abs(t2-t1);
      delta_t = delta_t*(p.t_weight_scale);
      
      if(i!=j) {
	//LR_couplings[i+j*p.surfaceVol] = AdS2p1Prop(circum[i%p.S1],
	//circum[j%p.S1],
	//delta_t, p);
      }      
      if(i!=j) {
	LR_couplings[i+j*p.surfaceVol] = sigmaL(circum[i%p.S1],
						circum[j%p.S1],
						delta_t);
      }      
      else LR_couplings[i+j*p.surfaceVol] = 0.0;
    }
  }

  //double unit_scale = 1.0/(AdS2p1Prop(circum[0], circum[1], 0, p));
  double unit_scale = 1.0/(sigmaL(circum[0], circum[1], 0));
  cout<<"Unit Scale = "<<unit_scale<<endl;
  
  //exit(0);

  
  for(int i=0; i<p.surfaceVol; i++) {
    for(int j=0; j<p.surfaceVol; j++) {
      //Rescale temporally
      if(j == tp(i,p)) {
	//LR_couplings[i + j*p.surfaceVol] *= celeritas;
      }
      //Rescale G(i;nn) = 1;
      LR_couplings[i + j*p.surfaceVol] *= unit_scale;
      
      //Additional \sigma scaling (c.f. Slava)
      LR_couplings[i + j*p.surfaceVol] = pow(LR_couplings[i + j*p.surfaceVol], p.sigma);

      if(i==j) LR_couplings[i+j*p.surfaceVol] = 0.0;
    }
  }
  for(int i=0; i<p.surfaceVol; i++) {
    for(int j=0; j<p.surfaceVol; j++) {
      //Check      
      if(i==0 && j%32 == 0) {
	cout<<i<<" "<<j/32<<" x "<<LR_couplings[i + (j/32)*p.surfaceVol]<<endl;
	cout<<i<<" "<<j/32<<" t "<<LR_couplings[i + j*p.surfaceVol]<<endl;       
      }      
    }
  }
}
