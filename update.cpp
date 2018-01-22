#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>

#include "update.h"
#include "util.h"

using namespace std;

extern int seed;
extern mt19937_64 rng;
extern uniform_real_distribution<double> unif;

//Basic utilites
// x = 0,1,..., p.S1-1  and
// t = 0,1,..., p.Lt-1
// i = x + p.S1*r so
// x = i % p.S1  and
// t = i/p.S1 
inline int xp(int i, Param p) {
  return (i + 1) % p.S1 + p.S1 * (i / p.S1);
}

inline int xm(int i, Param p) {
  return (i - 1 + p.S1) % p.S1 + p.S1 * (i / p.S1);
}

/* jump to next t shell */
inline int tp(int i, Param p) {
  return i % p.S1 + p.S1 * ((i / p.S1 + 1) % p.Lt);
}

/* return to last t shell */
inline int ttm(int i, Param p){
  // Some compiler have tm as reserved!
  return i % p.S1 + p.S1 * ((i / p.S1 - 1 + p.Lt) % p.Lt);
}

double actionPhiSqr(vector<double> &phi_arr, vector<int> &s, Param p,
		    double & KE, double & PE) {  

  KE = 0.0;
  PE = 0.0;
  double phi_sq;
  double phi;
  
  for (int i = 0; i < p.surfaceVol; i++)
    if (s[i] * phi_arr[i] < 0)
      printf("ERROR s and phi NOT aligned (AP) ! \n");
  
  //PE terms
  for (int i = 0; i < p.surfaceVol; i++) {    
    phi = phi_arr[i];
    phi_sq = phi*phi;
    
    PE += 0.25 * p.lambda * phi_sq*phi_sq;
    PE += 0.50 * p.musqr  * phi_sq;
  }
  
  // KE terms 
  for (int i = 0; i <p.surfaceVol; i++) {
    KE += 0.5 * (phi - phi_arr[xp(i,p)]) * (phi - phi_arr[xp(i,p)]);
    KE += 0.5 * (phi - phi_arr[tp(i,p)]) * (phi - phi_arr[tp(i,p)]);
  }	  
  
  return PE + KE;
}

int accept = 0;
int tries  = 0;

int metropolisUpdatePhiSqr(vector<double> &phi_arr, vector<int> &s, Param &p,
			   double &delta_mag_phi, int iter) {

  delta_mag_phi = 0.0;
  
  int s_old     = 0;
  int delta_mag = 0;
  
  double phi_new = 0.0;
  double phi_new_sq = 0.0;
  double phi = 0.0;
  double phi_sq = 0.0;

  double DeltaE = 0.0;
  
  for (int i = 0; i < p.surfaceVol; i++) {
    
    //Set some values we use a lot
    phi = phi_arr[i];
    DeltaE = 0.0;
    phi_new = phi + p.delta_phi * (2.0*unif(rng) - 1.0);
    phi_new_sq = phi_new*phi_new;
    phi_sq = phi*phi;
    
    if (s[i] * phi_arr[i] < 0) {
      printf("ERROR s and phi NOT aligned! (MUP)\n");
      exit(0);
    }

    //PE
    DeltaE += 0.25*p.lambda*(phi_new_sq*phi_new_sq - phi_sq*phi_sq);
    DeltaE += 0.50*p.musqr *(phi_new_sq            - phi_sq);
    
    //KE

    DeltaE += 2.0 * (phi_new_sq-phi_sq);
    DeltaE += 1.0 * (phi-phi_new)*(phi_arr[xp(i, p)] + phi_arr[xm(i, p)] +
				   phi_arr[tp(i, p)] + phi_arr[ttm(i, p)]);
    
    tries++;
    
    if(DeltaE < 0.0) {
      //  cout<< " Acepted  " << endl;
      s_old = s[i];
      delta_mag_phi += phi_new - phi_arr[i];
      phi_arr[i] = phi_new;
      accept += 1;
      s[i] = (phi_new > 0) ? 1 : -1;
      delta_mag += s[i] - s_old;
    }
    else if ( unif(rng)  < exp(-DeltaE)) {
      //  cout<< " Acepted  " << endl;
      s_old = s[i];
      delta_mag_phi += phi_new - phi_arr[i];
      phi_arr[i] = phi_new;
      accept += 1;
      s[i] = (phi_new > 0) ? 1 : -1;
      delta_mag += s[i] - s_old;
    }     
  }// end loop over lattice volume 
  
  // TUNING ACCEPTANCE 
  if (iter < p.n_therm/2 && iter % 500 == 0) {    
    if ((double) accept / (double) tries < 0.5) {
      p.delta_phi -= 0.001;
    } else {
      p.delta_phi += 0.001;
    }
  }

  if( iter < p.n_therm ) {
    //cout<<"At iter "<<iter<<" the Acceptance rate is "<<(double)accept/(double)tries<<endl;
    //cout<<"and delta_phi is "<<p.delta_phi<<endl;
  }
  return delta_mag;
}


//------------------------// 
// Wolff Cluster Routines //
//------------------------//
// declare functions and variables to implement Wolff algorithm
int wc_ave = 0;
int wc_size = 0;
int wc_calls = 0;

void wolffUpdatePhiSqr(vector<double> &phi, vector<int> &s, Param p,
		       double &delta_mag_phi, int iter) {

  wc_calls++;

  bool *site_cluster = new bool[p.surfaceVol];
  for (int i = 0; i < p.surfaceVol; i++)
    site_cluster[i] = false;

  // choose a random spin and grow a cluster
  int i = int(unif(rng) * p.surfaceVol);
  int clusterSpin = s[i];
  
  //This function is recursive and will call itself
  //until all four attempts in the lattice directions
  // (+x, -x, +t, -t) have failed to ncrease the cluster.
  wc_size = 1;
  clusterAddSqr(i, s, clusterSpin, site_cluster, phi, p);
  
  wc_ave += wc_size;

  if( (iter+1)%(p.n_skip) == 0 && iter < p.n_therm) {
    setprecision(4);
    cout<<"Ave. cluster size at iter "<<iter+1<<" = "<<wc_ave<<"/"<<wc_calls<<" = "<<1.0*wc_ave/wc_calls<<endl;
    //cout<<"S/T cluster growth ratio at iter "<<iter+1<<" = "<<1.0*wc_s_size/wc_ave<<":"<<1.0*wc_t_size/wc_ave<<endl;
  }
  
  delete[] site_cluster;    
}

void clusterAddSqr(int i, vector<int> &s, int clusterSpin,
		   bool *site_cluster, vector<double> &phi, Param p) {
  
  // The site belongs to the cluster, so flip it.
  site_cluster[i] = true;
  s[i] *= -1;
  phi[i] *= -1;

  // ferromagnetic coupling
  const double J = 1.0;
  
  // - If the (aligned) neighbor spin does not already belong to the
  // cluster, then try to add it to the cluster.
  // - If the site has already been added, then we may skip the test.

  //Backward in X
  if(!site_cluster[ xm(i,p) ] && s[xm(i,p)] == clusterSpin) {
    if (unif(rng) < 1 - exp(2*J*phi[i]*phi[xm(i,p)])) {
      wc_size++;
      //cout<<"->xm";
      clusterAddSqr(xm(i,p), s, clusterSpin, site_cluster, phi, p);
    }
  }

  //Forward in X
  if(!site_cluster[ xp(i,p) ] && s[xp(i,p)] == clusterSpin) {
    if(unif(rng) < 1 - exp(2*J*phi[i]*phi[xp(i,p)])) {
      wc_size++;
      //cout<<"->xp";
      clusterAddSqr(xp(i,p), s, clusterSpin, site_cluster, phi, p);
    }
  }

  //Forward in T
  if(!site_cluster[ tp(i,p) ] && s[tp(i,p)] == clusterSpin) {
    if(unif(rng) < 1 - exp(2*J*phi[i]*phi[tp(i,p)])) {
      wc_size++;
      //cout<<"->tp";
      clusterAddSqr(tp(i,p), s, clusterSpin, site_cluster, phi, p);
    }
  }

  //Backard in T 
  if(!site_cluster[ ttm(i,p) ] && s[ttm(i,p)] == clusterSpin) {  
    if(unif(rng) < 1 - exp(2*J*phi[i]*phi[ttm(i,p)])) {
      wc_size++;
      //cout<<"->tm";
      clusterAddSqr(ttm(i,p), s, clusterSpin, site_cluster, phi, p);
    }
  }  
}








int wc_t_size = 0;
int wc_s_size = 0;

double actionPhiAdS(vector<Vertex> &NodeList, vector<int> &s, Param p,
		    double &KE,  double &PE) {
  
  KE = 0.0;
  PE = 0.0;

  int interior = endNode(p.Levels-1,p);
  int disk     = p.AdSVol;
  double phi_sq;
  double phi;

  for (int i = 0; i < p.latVol; i++)
    if (s[i] * NodeList[i].phi < 0)
      printf("ERROR s and phi NOT aligned ! (AP AdS)\n");
  
  for (int i = 0; i < p.latVol; i++) {

    phi = NodeList[i].phi;
    phi_sq = phi*phi;

    //PE terms 
    if( i%disk > interior) {
      PE += 0.25 * p.lambda * phi_sq*phi_sq;
      PE += 0.5  * p.musqr  * phi_sq;
    }
    PE += 0.5*p.msqr*phi_sq;
    
    //KE terms
    //N.B. We must assess a link once, and once only.
    //To do this, we take the link to the 0th
    //neighbour, and all the forward links.
    //On the boundary, the fwdLinks value is set to 0, so only
    //the links connecting boundary points to other boundary
    //points are counted. 
    
    //spatial
    for(int q=0; q<NodeList[i].fwdLinks+1; q++) {      
      KE += 0.5*((phi - NodeList[NodeList[i].nn[q]].phi)*
		 (phi - NodeList[NodeList[i].nn[q]].phi)  );
    }    
    //temporal (the qth neighbour is always the forward time link.
    KE += 0.5*((phi - NodeList[NodeList[i].nn[p.q]].phi)*
	       (phi - NodeList[NodeList[i].nn[p.q]].phi)*
	       (NodeList[i].temporal_weight*NodeList[i].temporal_weight));
  }
  
  return PE + KE;
}

int metropolisUpdatePhiAdS(vector<Vertex> &NodeList, vector<int> &s,
			   Param &p, double &delta_mag_phi, int iter) {
  
  delta_mag_phi = 0.0;
  uniform_real_distribution<double> unif(0.0,1.0);

  int interior  = endNode(p.Levels-1,p);
  int disk      = p.AdSVol;  
  int s_old     = 0;
  int delta_mag = 0;
  int accept    = 0;
  int tries     = 0;

  double phi_new = 0.0;
  double phi_new_sq = 0.0;
  double phi = 0.0;
  double phi_sq = 0.0;

  double t_weight_sq = 0.0;
  double DeltaE = 0.0;
  
  for (int i = 0; i < p.latVol; i++) {

    //Set some values we use a lot
    t_weight_sq = NodeList[i].temporal_weight*NodeList[i].temporal_weight;
    phi = NodeList[i].phi;
    DeltaE = 0.0;
    phi_new = phi + p.delta_phi * (2.0*unif(rng) - 1.0);
    phi_new_sq = phi_new*phi_new;
    phi_sq = phi*phi;
    

    if (s[i] * phi < 0) {
      printf("ERROR s and phi NOT aligned! (MUP AdS)\n");
    }
    
    //PE
    if( i%disk > interior) {
      DeltaE += 0.25*p.lambda*(phi_new_sq*phi_new_sq - phi_sq*phi_sq);
      DeltaE += 0.50*p.musqr *(phi_new_sq            - phi_sq);
    }
    DeltaE += 0.50*p.msqr*(phi_new_sq - phi_sq);
    
    
    //KE
    for(int q=0; q<p.q; q++) {
      if(NodeList[i].nn[q] != -1) {
	DeltaE += 0.5 * (phi_new_sq - phi_sq + 2*NodeList[NodeList[i].nn[q]].phi*(phi - phi_new));
      }
    }
    for(int q=p.q; q<p.q+2; q++) {
      DeltaE += 0.5 * t_weight_sq * (phi_new_sq - phi_sq + 2*NodeList[NodeList[i].nn[q]].phi*(phi - phi_new));
    }
    
    tries += 1;    
    if(DeltaE < 0.0) {
      //  cout<< " Acepted  " << endl;
      s_old = s[i];
      delta_mag_phi += phi_new - NodeList[i].phi;
      NodeList[i].phi = phi_new;
      accept += 1;
      s[i] = (phi_new > 0) ? 1 : -1;
      delta_mag += s[i] - s_old;
    }
    else if ( unif(rng)  < exp(-DeltaE)) {
      //  cout<< " Acepted  " << endl;
      s_old = s[i];
      delta_mag_phi += phi_new - NodeList[i].phi;
      NodeList[i].phi = phi_new;
      accept += 1;
      s[i] = (phi_new > 0) ? 1 : -1;
      delta_mag += s[i] - s_old;
    }     
  }// end loop over lattice volume 
  
  // TUNING ACCEPTANCE 
  if (iter < p.n_therm/2) {
    if ((double) accept / (double) tries < 0.5) {
      p.delta_phi -= 0.01;
    } else {
      p.delta_phi += 0.01;
    }
  }

  if( iter < p.n_therm && (iter+1)%(100*p.n_wolff) == 0 ) {
    //cout<<"At iter "<<iter<<" the Acceptance rate is "<<(double)accept/(double)tries<<endl;
    //cout<<"and delta_phi is "<<p.delta_phi<<endl;
  }
  
  return delta_mag;
}

void wolffUpdatePhiAdS(vector<Vertex> &NodeList, vector<int> &s, Param p,
		       double &delta_mag_phi, int iter) {
  
  wc_calls++;
  
  bool *cluster = new bool[p.latVol];
  for (int i = 0; i < p.latVol; i++)
    cluster[i] = false;
  
  // choose a random spin and grow a cluster
  int i = int(unif(rng) * p.latVol);
  int cSpin = s[i];
  
  //This function is recursive and will call itself
  //until all q+2 attempts in the lattice directions
  //have failed to increase the cluster.
  wc_size = 1;
  clusterAddAdS(i, s, cSpin, cluster, NodeList, p);

  wc_ave += wc_size;

  if( (iter+1)%(p.n_skip) == 0 && iter < p.n_therm) {
    setprecision(4);
    cout<<"Ave. cluster size at iter "<<iter+1<<" = "<<wc_ave<<"/"<<wc_calls<<" = "<<1.0*wc_ave/wc_calls<<endl;
    //cout<<"S/T cluster growth ratio at iter "<<iter+1<<" = "<<1.0*wc_s_size/wc_ave<<":"<<1.0*wc_t_size/wc_ave<<endl;
  }
  
  delete[] cluster;    
}

void clusterAddAdS(int i, vector<int> &s, int cSpin,
		   bool *cluster, vector<Vertex> &NodeList, Param p) {
  
  // Mark the spin as belonging to the cluster and flip it
  cluster[i] = true;
  s[i] *= -1;
  NodeList[i].phi *= -1;
  
  // ferromagnetic coupling
  const double J = +1.0;
  double t_weight = 1.0;
  
  // if the neighbor spin does not already belong to the
  // cluster, then try to add it to the cluster
  for(int q=0; q<p.q+2; q++) {
    
    //qth nearest neighbour. If boundary node, nn_q = -1
    int nn_q = NodeList[i].nn[q];    
    if(q >= p.q) t_weight = NodeList[i].temporal_weight;
    
    //Check if the candidate node is both not in the cluster, and not
    //beyond the boundary.    
    if(nn_q != -1) {
      if(!cluster[NodeList[nn_q].pos]) {
	if(s[NodeList[nn_q].pos] == cSpin &&
	   unif(rng) < 1 - exp(2*J*t_weight*NodeList[i].phi*NodeList[nn_q].phi)) {
	  wc_size++;
	  //if(q<p.q)  wc_s_size++;
	  //if(q>=p.q) wc_t_size++;	  
	  clusterAddAdS(NodeList[nn_q].pos, s, cSpin, cluster, NodeList, p);
	}
      }
    }
  }
}

//-------------AdS Monte Carlo Update  ----------------//
void runMonteCarloAdS(vector<Vertex> &NodeList, Param p) {
  
  p.AdSVol = endNode(p.Levels,p)+1;
  p.Lt = p.t;
  p.latVol = p.AdSVol * p.Lt;
  p.surfaceVol = endNode(p.Levels,p) - endNode(p.Levels-1,p);
  
  double KE = 0.0, PE = 0.0;
  double Etot = 0.0;
  double mag_phi = 0.0;
  double delta_mag_phi = 0.0;

  vector<int> s(p.latVol,0.0);
  
  for(int i = 0;i < p.latVol; i++) {
    NodeList[i].phi = 2.0*unif(rng) - 1.0;
    s[i] = (NodeList[i].phi > 0) ? 1:-1;
    mag_phi += NodeList[i].phi;
  }
  
  cout<<"Levels = "<<p.Levels<<" Lattice Size "<<p.AdSVol<<" x "<<p.Lt<<" = ";
  cout<<p.latVol<<" Initial Mag = "<<mag_phi<<endl;

  Etot = actionPhiAdS(NodeList, s, p, KE, PE);

  cout<<"K = "<<KE<<" U = "<<PE<<" K + U = "<<KE+PE<<endl;
  cout<<"Etot = "<<Etot<<endl;
  cout<<"Etot - (K+U) = "<<Etot-(KE+PE)<<endl;

  //Thermalisation
  
  for(int iter = 0;iter < p.n_therm; iter++) {
    //if p.n_wolff is set to zero, no Wolff steps occur.
    if((iter+1)%p.n_wolff == 0 && p.n_wolff != 0) {
      metropolisUpdatePhiAdS(NodeList, s, p, delta_mag_phi, iter);
    }
    wolffUpdatePhiAdS(NodeList, s, p, delta_mag_phi, iter);

    if((iter+1)%p.n_skip == 0) cout<<"Therm sweep "<<iter+1<<endl;
  }

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
  double aveE2    = 0.0;
  double avePhiAb = 0.0;
  double avePhi   = 0.0;
  double avePhi2  = 0.0;
  double avePhi4  = 0.0;
  double MagPhi   = 0.0;

  //Density of lattice sites
  double rhoVol  = 1.0/(double)p.latVol;
  double rhoVol2 = rhoVol*rhoVol;

  //correlation function arrays. One keeps the running average,
  //one is temporary to pass around single measurent values. 
  double **corr_tmp = (double**)malloc((p.surfaceVol/2)*sizeof(double*));
  double **corr_ave = (double**)malloc((p.surfaceVol/2)*sizeof(double*));
  for(int i=0; i<p.surfaceVol/2; i++) {
    corr_tmp[i] = (double*)malloc((p.Lt/2)*sizeof(double));
    corr_ave[i] = (double*)malloc((p.Lt/2)*sizeof(double));
  }
  for(int i=0; i<p.surfaceVol/2; i++)
    for(int j=0; j<p.Lt/2; j++) 
      corr_ave[i][j] = 0.0;
  
  double **phi_sq_arr = (double**)malloc(p.surfaceVol*sizeof(double*));
  for(int i=0; i<p.surfaceVol; i++) {
    phi_sq_arr[i] = (double*)malloc(p.Lt*sizeof(double));
    for(int j=0; j<p.Lt; j++) phi_sq_arr[i][j] = 0.0;
  }

  int idx = 0;
  double norm;
  
  for(int iter = p.n_therm; iter < p.n_therm + p.n_skip*p.n_meas; iter++) {
    
    if((iter+1)%p.n_wolff == 0 && p.n_wolff != 0)
      metropolisUpdatePhiAdS(NodeList, s, p, delta_mag_phi, iter);
    else wolffUpdatePhiAdS(NodeList, s, p, delta_mag_phi, iter);
    
    //Take measurements.
    if((iter+1) % p.n_skip == 0) {
      
      int offset = endNode(p.Levels-1,p) + 1;
      for(int i = 0;i < p.surfaceVol; i++) {
	for(int j=0; j<p.Lt; j++) 
	  phi_sq_arr[i][j] += pow(NodeList[i + offset + p.AdSVol*j].phi,2);
      }
      
      tmpE     = actionPhiAdS(NodeList, s, p, KE, PE);
      aveE    += rhoVol*tmpE;
      aveE2   += rhoVol2*tmpE*tmpE;

      MagPhi = 0.0;
      for(int i = 0;i < p.latVol; i++) MagPhi += NodeList[i].phi;
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
      cout<<"Measurement "<<(iter+1)/p.n_skip<<" Sweep "<<iter+1<<endl;
      cout<<"Ave Energy= "<<aveE*norm<<endl;
      cout<<"Ave |phi| = "<<avePhiAb*norm<<endl;
      cout<<"Ave phi   = "<<avePhi*norm<<endl;
      cout<<"Ave phi^2 = "<<avePhi2*norm<<endl;
      cout<<"Ave phi^4 = "<<avePhi4*norm<<endl;
      cout<<"Suscep    = "<<(avePhi2*norm-pow(avePhiAb*norm,2))/rhoVol<<endl;
      cout<<"Spec Heat = "<<(aveE2*norm-pow(aveE*norm,2))/rhoVol<<endl;
      cout<<"Binder    = "<<1.0-avePhi4/(3.0*avePhi2*avePhi2*norm)<<endl;

      //Visualisation tools
      //visualiserAdS(NodeList, avePhiAb*norm, p);
      //visualiserPhi2(phi_sq_arr, p, idx);      

      //Calculate correlaton functions and update the average.
      correlatorsAdS(corr_tmp, corr_ave, idx, NodeList, avePhi*norm, p);
      ofstream myfile("correlators.dat");
      for(int i=0; i<p.Lt/2; i++) {
	myfile << i << " " << corr_tmp[0][i] << endl;
      }
      myfile.close();

      
      
      //if(idx%50 == 0) corr_eigs(corr_run, p);
      
    }
  }

  //autocorrelation(PhiAb_arr, avePhiAb, p.n_meas);
  
  //correlators(corr_run, corr_ave, idx, phi_cyl, avePhi*norm, p);
  //corr_eigs(corr_run, p);
  
  //free(corr_run);
  //free(corr_ave);
}

//-------------Sqr Monte Carlo Update  ----------------//
void runMonteCarloSqr(vector<Vertex> &NodeList, Param p) {
  
  p.Lt = p.t;
  p.S1 = endNode(p.Levels,p) - endNode(p.Levels-1,p);
  p.surfaceVol = p.S1 * p.Lt;

  double KE = 0.0, PE = 0.0;
  double Etot = 0.0;
  double mag_phi = 0.0;
  double delta_mag_phi = 0.0;

  vector<double> phi(p.surfaceVol,0.0);
  vector<int> s(p.surfaceVol,0.0);
  
  for(int i = 0;i < p.surfaceVol; i++) {
    phi[i] = 2.0*unif(rng) - 1.0;
    s[i] = (phi[i] > 0) ? 1:-1;
    mag_phi += phi[i];
  }
  
  cout<<"Levels = "<<p.Levels<<" Lattice Size "<<p.S1<<" x "<<p.Lt<<" = ";
  cout<<p.surfaceVol<<" Initial Mag = "<<mag_phi<<endl;

  Etot = actionPhiSqr(phi, s, p, KE, PE);
  
  cout<<"K = "<<KE<<" U = "<<PE<<" K + U = "<<KE+PE<<endl;
  cout<<"Etot = "<<Etot<<endl;
  cout<<"Etot - (K+U) = "<<Etot-(KE+PE)<<endl;
  
  //Thermalisation
  for(int iter = 0;iter < p.n_therm; iter++) {
    //if p.n_wolff is set to zero, no Wolff steps occur.
    if((iter+1)%p.n_wolff == 0 && p.n_wolff != 0) {
      metropolisUpdatePhiSqr(phi, s, p, delta_mag_phi, iter);
    }
    wolffUpdatePhiSqr(phi, s, p, delta_mag_phi, iter);

    if((iter+1)%p.n_skip == 0) cout<<"Therm sweep "<<iter+1<<endl;
  }

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
    
    if((iter+1)%p.n_wolff == 0 && p.n_wolff != 0) {
      metropolisUpdatePhiSqr(phi, s, p, delta_mag_phi, iter);
    }
    wolffUpdatePhiSqr(phi, s, p, delta_mag_phi, iter);
    
    
    //Take measurements.
    if((iter+1) % p.n_skip == 0) {

      for(int i = 0;i < p.S1; i++) {
	for(int j=0; j<p.Lt; j++) 
	  phi_sq_arr[i][j] += pow(phi[i + p.S1*j],2);
      }
      
      tmpE     = actionPhiSqr(phi, s, p, KE, PE);
      aveE    += rhoVol*tmpE;
      aveE2   += rhoVol2*tmpE*tmpE;

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
      cout<<"Measurement "<<(iter+1)/p.n_skip<<" Sweep "<<iter+1<<endl;
      cout<<"Ave Energy= "<<aveE*norm<<endl;
      cout<<"Ave |phi| = "<<avePhiAb*norm<<endl;
      cout<<"Ave phi   = "<<avePhi*norm<<endl;
      cout<<"Ave phi^2 = "<<avePhi2*norm<<endl;
      cout<<"Ave phi^4 = "<<avePhi4*norm<<endl;
      cout<<"Suscep    = "<<(avePhi2*norm-pow(avePhiAb*norm,2))/rhoVol<<endl;
      cout<<"Spec Heat = "<<(aveE2*norm-pow(aveE*norm,2))/rhoVol<<endl;
      cout<<"Binder    = "<<1.0-avePhi4/(3.0*avePhi2*avePhi2*norm)<<endl;

      //Visualisation tools
      //visualiserSqr(phi, avePhiAb*norm, p);
      //visualiserPhi2(phi_sq_arr, p, idx);      

      //Calculate correlaton functions and update the average.
      correlatorsSqr(corr_tmp, corr_ave, idx, phi, avePhi*norm, p);
      ofstream myfile("correlators.dat");
      for(int i=0; i<p.Lt/2; i++) {
	myfile << i << " " << corr_tmp[0][i] << endl;
      }
      myfile.close();
      //if(idx%50 == 0) corr_eigs(corr_run, p);
      
    }
  }

  //autocorrelation(PhiAb_arr, avePhiAb, p.n_meas);
  
  //correlators(corr_run, corr_ave, idx, phi_cyl, avePhi*norm, p);
  //corr_eigs(corr_run, p);
  
  //free(corr_run);
  //free(corr_ave);
}

