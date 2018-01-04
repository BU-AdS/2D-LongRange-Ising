#ifndef UPDATE_2DADS_H
#define UPDATE_2DADS_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "graph.h"

double action_phi_AdS(vector<Vertex> &NodeList, vector<int> &s, Param p, double &KE,  double &PE) {
  
  int i;
  KE = 0.0;
  PE = 0.0;
  
  p.latVol = p.SurfaceVol;
  
  for (i = 0; i < p.latVol; i++)
    if (s[i] * NodeList[i].phi < 0)
      printf("ERROR s and phi NOT aligned ! (AP AdS)\n");
  
  // PE  terms 
  for (i = 0; i < p.latVol; i++) {
    PE += 0.25 * p.lambda * NodeList[i].phi*NodeList[i].phi*NodeList[i].phi*NodeList[i].phi;
    PE -= 0.5  * p.musqr  * NodeList[i].phi*NodeList[i].phi;
  }
  
  // KE Terms on the AdS and along the cylinder
  for (i = 0; i <p.latVol; i++) {
    for(int q=0; q<NodeList[i].fwdLinks+1; q++) {
      if(NodeList[NodeList[i].nn[q]].pos != -1) {
	KE += 0.5*((NodeList[i].phi - NodeList[NodeList[i].nn[q]].phi)*
		   (NodeList[i].phi - NodeList[NodeList[i].nn[q]].phi)  );
      }
    }
    KE += 0.5*((NodeList[i].phi - NodeList[NodeList[i].nn[p.q]].phi)*
	       (NodeList[i].phi - NodeList[NodeList[i].nn[p.q]].phi)*
	       (NodeList[i].temporal_weight*NodeList[i].temporal_weight));
  }
  
  return PE + KE;
}


// Prob = EXP[ - E]  AcceptProb = Min[1, exp[ - DeltaE]]

int metropolis_update_phi_AdS(vector<Vertex> &NodeList, vector<int> &s, Param p,
			      double & delta_mag_phi, int iter) {
  int  s_old;
  uniform_real_distribution<double> unif(0.0,1.0);
  int delta_mag = 0;
  int accept = 0;
  int tries = 0;
  // static double delta_max = 10.; 
  // static double delta_min = 0.0; 
  
  static double delta_phi = 1.5;
  
  delta_mag_phi = 0.;
  
  for (int i = 0; i < p.latVol; i++) {
    if (s[i] * NodeList[i].phi < 0) {
      printf("ERROR s and phi NOT aligned! (MUP AdS)\n");
      exit(0);
    }
    
    double DeltaE = 0.0;
    double phi_new = NodeList[i].phi + delta_phi * (2.0*unif(rng) - 1.0);
    
    if(p.verbosity) {
      // Put into test.h later ? << FIXME
      printf("\n Entering metropolis for i = %d \n", i);
      double KE = 0.0;
      double PE = 0.0;
      double DeltaKE = 0.0;
      double DeltaPE = 0.0;
      double Energy = action_phi_AdS(NodeList, s, p, KE,  PE);
      double DeltaEnergy = 0.0;
      DeltaKE = KE;
      DeltaPE = PE; 
      printf(" Energy = %f \n",Energy);
      
      double phi_tmp = NodeList[i].phi;
      NodeList[i].phi = phi_new;
      DeltaEnergy = action_phi_AdS(NodeList, s,  p, KE,  PE) - Energy;
      NodeList[i].phi = phi_tmp;
      
      DeltaKE = KE - DeltaKE;
      DeltaPE = PE - DeltaPE;
      printf(" DeltaE = %f DeltaKE = %f DeltaPE = %f DeltaKE+DeltaPE = %f \n",
	     DeltaEnergy, DeltaKE, DeltaPE, DeltaKE + DeltaPE);
    }
    
    //PE    
    DeltaE += 0.25*p.lambda*(phi_new*phi_new*phi_new*phi_new - NodeList[i].phi*NodeList[i].phi*NodeList[i].phi*NodeList[i].phi);
    DeltaE += -0.5*p.musqr *(phi_new*phi_new                 - NodeList[i].phi*NodeList[i].phi);
    
    //KE
    for(int q=0; q<p.q; q++) {
      DeltaE += 0.5*((phi_new - NodeList[NodeList[i].nn[q]].phi)*
		     (phi_new - NodeList[NodeList[i].nn[q]].phi) -
		     (NodeList[i].phi - NodeList[NodeList[i].nn[q]].phi)*
		     (NodeList[i].phi - NodeList[NodeList[i].nn[q]].phi));
    }
    for(int q=p.q; q<p.q+2; q++) {
      DeltaE += (0.5*((phi_new - NodeList[NodeList[i].nn[q]].phi)*
		      (phi_new - NodeList[NodeList[i].nn[q]].phi) -
		      (NodeList[i].phi - NodeList[NodeList[i].nn[q]].phi)*
		      (NodeList[i].phi - NodeList[NodeList[i].nn[q]].phi))*
		 NodeList[i].temporal_weight*NodeList[i].temporal_weight);
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
  if (iter < p.n_therm) {
    if ((double) accept / (double) tries < .5) {
      delta_phi -= 0.01;
    } else {
      delta_phi += 0.01;
    }
  }
  
  return delta_mag;
}

//------------------------// 
// Wolff Cluster Routines //
//------------------------//

// declare functions to implement Wolff algorithm
void cluster_add_AdS(int i, vector<int> &s, int clusterSpin,
		     bool *cluster, vector<Vertex> &NodeList, Param p);

void wolff_update_phi_AdS(vector<Vertex> &NodeList, vector<int> &s, Param p,
			  double &delta_mag_phi, int iter) {
  
  bool *cluster = new bool[p.SurfaceVol];
  for (int i = 0; i < p.SurfaceVol; i++)
    cluster[i] = false;
  
  // choose a random spin and grow a cluster
  int i = int(unif(rng) * p.latVol);
  int clusterSpin = s[i];

  //This function is recursive and will call itself
  //until all four attempts in the lattice direactions
  // (+x, -x, +t, -t) have failed to ncrease the cluster.  
  cluster_add_AdS(i, s, clusterSpin, cluster, NodeList, p);
  
  int size = 0;
  for(int i=0; i<p.SurfaceVol; i++) if(cluster[i] == true) size++;
  if(p.verbosity) {
    cout<<"cluster size at iter "<<iter+1<<" = "<<size;
    cout<<" = "<<100.00*(double)size/p.SurfaceVol<<"%"<<endl;
  }  
  delete[] cluster;  
}

void cluster_add_AdS(int i, vector<int> &s, int clusterSpin,
		     bool *cluster, vector<Vertex> &NodeList, Param p) {
  
  // Mark the spin as belonging to the cluster and flip it
  cluster[i] = true;
  s[i] *= -1;
  NodeList[i].phi *= -1;
  
  // ferromagnetic coupling
  double J = +1;
  
  // if the neighbor spin does not already belong to the
  // cluster, then try to add it to the cluster
  for(int q=0; q<p.q+2; q++)
    if(!cluster[NodeList[NodeList[i].nn[q]].pos]) {
      if(s[NodeList[NodeList[i].nn[q]].pos] == clusterSpin && unif(rng) < 1 - exp(2*J*NodeList[i].phi*NodeList[NodeList[i].nn[q]].phi)) {
	cluster_add_AdS(NodeList[NodeList[i].nn[q]].pos, s, clusterSpin, cluster, NodeList, p);
      }
    }
}

#endif
