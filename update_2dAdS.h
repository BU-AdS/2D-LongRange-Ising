#ifndef UPDATE_2DADS_H
#define UPDATE_2DADS_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "graph.h"

double action_phi_AdS(vector<Vertex> &NodeList, vector<int> &s, Param p,
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

    // PE  terms 
    if( i%disk > interior) {
      PE += 0.25 * p.lambda * phi_sq*phi_sq;
      PE += 0.5  * p.musqr  * phi_sq;
    }
    PE += 0.5*p.msqr*phi_sq;
    
    
    // KE Terms
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


int metropolis_update_phi_AdS(vector<Vertex> &NodeList, vector<int> &s,
			      Param p, double & delta_mag_phi, int iter) {
  
  delta_mag_phi = 0.0;
  uniform_real_distribution<double> unif(0.0,1.0);

  int interior  = endNode(p.Levels-1,p);
  int disk      = p.AdSVol;  
  int s_old     = 0;
  int delta_mag = 0;
  int accept    = 0;
  int tries     = 0;
  
  double delta_phi = 1.5;
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
    phi_new = phi + delta_phi * (2.0*unif(rng) - 1.0);
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
    DeltaE += 0.5*p.msqr*(phi_new_sq - phi_sq);
    
    
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
      delta_phi -= 0.01;
    } else {
      delta_phi += 0.01;
    }
  }

  if( (iter+1)%p.n_therm == 0 ) {
    cout<<"At Thermalisation, the Acceptance rate is "<<(double)accept/(double)tries<<endl;
    cout<<"and delta_phi is "<<delta_phi<<endl;
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

void cluster_add_AdS(int i, vector<int> &s, int clusterSpin,
		     bool *cluster, vector<Vertex> &NodeList, Param p);

void wolff_update_phi_AdS(vector<Vertex> &NodeList, vector<int> &s, Param p,
			  double &delta_mag_phi, int iter) {
  
  wc_calls++;
  
  bool *cluster = new bool[p.latVol];
  for (int i = 0; i < p.latVol; i++)
    cluster[i] = false;
  
  // choose a random spin and grow a cluster
  int i = int(unif(rng) * p.latVol);
  int cSpin = s[i];
  //Here we flip the spin of this node because the
  //Cluster add routines flip the spin automatically.
  //Hence the spin starts of as +/-, we flip the spin
  //to -/+, and the true spin +/- is restotred in the
  //course of the algorithm.
  s[i] *= -1;
  NodeList[i].phi *= -1;

  //This function is recursive and will call itself
  //until all q+2 attempts in the lattice directions
  //have failed to increase the cluster.
  wc_size = 0;
  cluster_add_AdS(i, s, cSpin, cluster, NodeList, p);

  wc_ave += wc_size;

  //cout<<"iter="<<iter+1<<endl;

  if( (iter+1) % p.n_skip == 0) {
    setprecision(4);
    cout<<"Ave. cluster size at iter "<<iter+1<<" = "<<((100.00*wc_ave)/wc_calls)/p.latVol<<" %"<<endl;
  }  
  delete[] cluster;    
}

void cluster_add_AdS(int i, vector<int> &s, int cSpin,
		     bool *cluster, vector<Vertex> &NodeList, Param p) {
  
  // Mark the spin as belonging to the cluster and flip it
  cluster[i] = true;
  s[i] *= -1;
  NodeList[i].phi *= -1;
  
  // ferromagnetic coupling
  double J = +1.0;
  
  // if the neighbor spin does not already belong to the
  // cluster, then try to add it to the cluster
  for(int q=0; q<p.q+2; q++) {
    
    //qth nearest neighbour. If boundary node, nn_q = -1
    int nn_q = NodeList[i].nn[q];
    
    //Check if the candidate node is both not in the cluster, and not
    //beyond the boundary.    
    if(nn_q != -1) {
      if(!cluster[NodeList[nn_q].pos]) {
	if(s[NodeList[nn_q].pos] == cSpin &&
	   unif(rng) < 1 - exp(-2*J*NodeList[i].phi*NodeList[nn_q].phi)) {
	  wc_size++;
	  cluster_add_AdS(NodeList[nn_q].pos, s, cSpin, cluster, NodeList, p);
	}
      }
    }
  }
}

/*
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
*/

#endif
