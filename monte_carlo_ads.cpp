#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>

#include "monte_carlo_ads.h"
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
int ads_wcave = 0;
int ads_wcsize = 0;
int ads_wccalls = 0;
int ads_wct_size = 0;
int ads_wcs_size = 0;

double actionAdS(vector<Vertex> &NodeList, int *s, Param p,
		 double &KE,  double &PE) {
  
  KE = 0.0;
  PE = 0.0;
  double phi_sq;
  double phi;
  int interior = endNode(p.Levels-1,p);
  int disk     = p.AdSVol;

  
  for (int i = 0; i < p.latVol; i++)
    if (s[i] * NodeList[i].phi < 0)
      printf("ERROR s and phi NOT aligned ! (actionPhi AdS)\n");
  
  for (int i = 0; i < p.latVol; i++) {
    phi = NodeList[i].phi;
    phi_sq = phi*phi;

    //2D Ising terms
    if( i%disk > interior) {
      PE += 0.25 * p.lambda * phi_sq*phi_sq;
      PE += 0.5  * p.musqr  * phi_sq;
      
      //Spatial: q=0 and q=fwdLinks+1 are on the same
      //level. These are 'Ising' terms and should not be
      //rescaled?
      for(int q=0; q<NodeList[i].fwdLinks+1; q++) {
	if(NodeList[NodeList[i].nn[q]].pos != -1) {
	  KE += 0.5*((phi - NodeList[NodeList[i].nn[q]].phi)*
		     (phi - NodeList[NodeList[i].nn[q]].phi));
	}
      }
      //temporal (the qth neighbour is always the forward time link)
      KE += 0.5*((phi - NodeList[NodeList[i].nn[p.q]].phi)*
		 (phi - NodeList[NodeList[i].nn[p.q]].phi))*NodeList[i].temporal_weight;
      
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
		 (phi - NodeList[NodeList[i].nn[q]].phi));
    }    
    //temporal (the qth neighbour is always the forward time link)
    KE += 0.5*((phi - NodeList[NodeList[i].nn[p.q]].phi)*
	       (phi - NodeList[NodeList[i].nn[p.q]].phi)*
	       (NodeList[i].temporal_weight*NodeList[i].temporal_weight));
  }
  
  return PE + KE;
}

int ads_accept = 0;
int ads_tries  = 0;

int metropolisUpdateAdS(vector<Vertex> &NodeList, int *s,
			Param &p, double &delta_mag_phi, int iter) {
  
  delta_mag_phi = 0.0;
  uniform_real_distribution<double> unif(0.0,1.0);

  int interior  = endNode(p.Levels-1,p);
  int disk      = p.AdSVol;  
  int s_old     = 0;
  int delta_mag = 0;

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
      //cout<<"outer:"<<DeltaE<<endl;      
    }
    DeltaE += 0.50*p.msqr*(phi_new_sq - phi_sq);
    //cout<<"inner="<<DeltaE<<endl;
    
    
    //KE
    for(int q=0; q<p.q; q++) {
      if(NodeList[i].nn[q] != -1) {
	DeltaE += 0.5 * (phi_new_sq - phi_sq + 2*NodeList[NodeList[i].nn[q]].phi*(phi - phi_new));
      }
    }
    for(int q=p.q; q<p.q+2; q++) {
      DeltaE += 0.5 * t_weight_sq * (phi_new_sq - phi_sq + 2*NodeList[NodeList[i].nn[q]].phi*(phi - phi_new));
    }
    
    ads_tries += 1;    
    if(DeltaE < 0.0) {
      //  cout<< " Acepted  " << endl;
      s_old = s[i];
      delta_mag_phi += phi_new - NodeList[i].phi;
      NodeList[i].phi = phi_new;
      ads_accept += 1;
      s[i] = (phi_new > 0) ? 1 : -1;
      delta_mag += s[i] - s_old;
    }
    else if ( unif(rng)  < exp(-DeltaE)) {
      //  cout<< " Acepted  " << endl;
      s_old = s[i];
      delta_mag_phi += phi_new - NodeList[i].phi;
      NodeList[i].phi = phi_new;
      ads_accept += 1;
      s[i] = (phi_new > 0) ? 1 : -1;
      delta_mag += s[i] - s_old;
    }     
  }// end loop over lattice volume 
  
  // TUNING ACCEPTANCE 
  if (iter < p.n_therm/2 && (iter+1) % p.n_skip/10 == 0) {
    if ((double) ads_accept / (double) ads_tries < 0.5) {
      p.delta_phi -= 0.001;
    } else {
      p.delta_phi += 0.001;
    }
    if(p.n_wolff*1.0*ads_wcave/ads_wccalls < p.latVol && iter > p.n_skip) {
      p.n_wolff++;
    } else {
      p.n_wolff--;
      if(p.n_wolff < 3) p.n_wolff++;
    }
  }

  if( iter < p.n_therm && (iter+1)%(100*p.n_wolff) == 0 ) {
    //cout<<"At iter "<<iter<<" the Acceptance rate is "<<(double)ads_accept/(double)ads_tries<<endl;
    //cout<<"and delta_phi is "<<p.delta_phi<<endl;
  }
  
  return delta_mag;
}

void wolffUpdateAdS(vector<Vertex> &NodeList, int *s, Param p,
		    double &delta_mag_phi, int iter) {
  
  ads_wccalls++;
  
  bool *cluster = new bool[p.latVol];
  for (int i = 0; i < p.latVol; i++)
    cluster[i] = false;
  
  // choose a random spin and grow a cluster
  int i = int(unif(rng) * p.latVol);
  int cSpin = s[i];
  
  //This function is recursive and will call itself
  //until all q+2 attempts in the lattice directions
  //have failed to increase the cluster.
  ads_wcsize = 1;
  clusterAddAdS(i, s, cSpin, cluster, NodeList, p);
  //cout<<"ads_wcsize="<<ads_wcsize<<endl;
  ads_wcave += ads_wcsize;

  if( iter%p.n_skip == 0 && iter < p.n_therm) {
    setprecision(4);
    cout<<"Using "<<p.n_wolff<<" Wolff hits."<<endl; 
    cout<<"Ave. cluster size at iter "<<iter<<" = "<<ads_wcave<<"/"<<ads_wccalls<<" = "<<1.0*ads_wcave/ads_wccalls<<endl;
    cout<<"S/T cluster growth ratio at iter "<<iter<<" = "<<1.0*ads_wcs_size/ads_wcave<<":"<<1.0*ads_wct_size/ads_wcave<<endl;
  }
  
  delete[] cluster;    
}

void clusterAddAdS(int i, int *s, int cSpin,
		   bool *cluster, vector<Vertex> &NodeList, Param p) {
  
  // Mark the spin as belonging to the cluster and flip it
  cluster[i] = true;
  s[i] *= -1;
  NodeList[i].phi *= -1;
  
  // ferromagnetic coupling
  const double J = +1.0;
  double t_weight = 1.0;

  double prob = 0.0;
  double rand = 0.0;
  
  // if the neighbor spin does not already belong to the
  // cluster, then try to add it to the cluster
  //cout<<"Base="<<i<<" ";
  for(int q=p.q+1; q >= 0; q--) {  
    //for(int q=0; q<p.q+2; q++) {
    
    //qth nearest neighbour. If boundary node, nn_q = -1
    int nn_q = NodeList[i].nn[q];
    //cout<<nn_q<<" "<<NodeList[nn_q].pos<<" ";

    if(q >= p.q) t_weight = NodeList[i].temporal_weight;
    
    //Check if the candidate node is both not in the cluster, and not
    //beyond the boundary.    
    if(nn_q != -1) {
      if(!cluster[NodeList[nn_q].pos]) {
	rand = unif(rng);
	prob = 1 - exp(2*J*t_weight*NodeList[i].phi*NodeList[nn_q].phi);
	if(s[NodeList[nn_q].pos] == cSpin && prob < rand) {
	  //cout<<prob<<" "<<rand<<" ";
	  ads_wcsize++;
	  if(q<p.q)  ads_wcs_size++;
	  if(q>=p.q) ads_wct_size++;
	  //cout<<"->"<<q<<" ";
	  clusterAddAdS(NodeList[nn_q].pos, s, cSpin, cluster, NodeList, p);
	}
      }
    }
  }
  //cout<<"End"<<endl;
}

//------------- Monte Carlo Update  ----------------//
void runMonteCarloAdS(vector<Vertex> &NodeList, Param p) {
  
  double KE = 0.0, PE = 0.0;
  double mag_phi = 0.0;
  double delta_mag_phi = 0.0;

  int *s = (int*)malloc(p.latVol*sizeof(int));
  for(int i = 0;i < p.latVol; i++) {
    NodeList[i].phi = 2.0*unif(rng) - 1.0;
    s[i] = (NodeList[i].phi > 0) ? 1:-1;
    mag_phi += NodeList[i].phi;
  }
  thermaliseAdS(NodeList, s, p, delta_mag_phi);
  
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
  double rhoVol  = 1.0/(double)p.latVol;
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
      metropolisUpdateAdS(NodeList, s, p, delta_mag_phi, iter);
    } else {
      wolffUpdateAdS(NodeList, s, p, delta_mag_phi, iter);
    }
    
    //Take measurements.
    if((iter+1) % p.n_skip == 0) {
      
      int offset = endNode(p.Levels-1,p) + 1;
      for(int i = 0;i < p.S1; i++) {
	for(int j=0; j<p.Lt; j++) 
	  phi_sq_arr[i][j] += pow(NodeList[i + offset + p.AdSVol*j].phi,2);
      }
      
      tmpE   = actionAdS(NodeList, s, p, KE, PE);
      aveKE += rhoVol*KE;
      avePE += rhoVol*PE;
      aveE  += rhoVol*tmpE;
      aveE2 += rhoVol2*tmpE*tmpE;

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
      for(int i=endNode(p.Levels-1,p)+1; i<endNode(p.Levels,p)+1; i++) cout<<NodeList[i].phi<<" ";
      cout<<endl;
      
      //Visualisation tools
      //visualiserAdS(NodeList, avePhiAb*norm, p);
      //visualiserPhi2(phi_sq_arr, p, idx);      

      //Calculate correlaton functions and update the average.
      correlators(corr_tmp, corr_ave, idx, NodeList, avePhi*norm, p);
      
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

  free(s);
}

void thermaliseAdS(vector<Vertex> &NodeList, int *s, 
		   Param p, double &delta_mag_phi) {
  
  for(int iter = 0; iter < p.n_therm; iter++) {
    for(int i=0; i<p.n_wolff; i++) {
      wolffUpdateAdS(NodeList, s, p, delta_mag_phi, iter);
      if((iter+1)%p.n_skip == 0) cout<<"Therm sweep "<<iter+1<<endl;
      iter++;
    }
    metropolisUpdateAdS(NodeList, s, p, delta_mag_phi, iter);    
    if((iter+1)%p.n_skip == 0) cout<<"Therm sweep "<<iter+1<<endl;
  }
}
