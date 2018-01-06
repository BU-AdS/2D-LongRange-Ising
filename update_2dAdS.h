#ifndef UPDATE_2DADS_H
#define UPDATE_2DADS_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "graph.h"

double action_phi_AdS(vector<Vertex> &NodeList, vector<int> &s, Param p, double &KE,  double &PE) {
  
  KE = 0.0;
  PE = 0.0;

  int interior = endNode(p.Levels-1,p);
  int disk     = endNode(p.Levels,p)+1;
  double phi_sq;
  
  for (int i = 0; i < p.latVol; i++)
    if (s[i] * NodeList[i].phi < 0)
      printf("ERROR s and phi NOT aligned ! (AP AdS)\n");
  

  for (int i = 0; i < p.latVol; i++) {
    // PE  terms 
    phi_sq = NodeList[i].phi*NodeList[i].phi;
    if( i%disk > interior) {
      PE += 0.25 * p.lambda * phi_sq*phi_sq;
      PE += 0.5  * p.musqr  * phi_sq;
    }
    PE += 0.5*p.msqr*phi_sq;
    
    
    // KE Terms
    //N.B. We must assess a link once, and once only.
    //To do this, we take the link the 0th
    //neighbour, and all the forward links.
    //On the interior, a loop from .nn[0] to .nn[fwdLinks]
    //will give is all links connected to internal nodes
    //On the boundary, the fwdLinks value is set to 0, so only
    //the links connecting the boundary point to another boundary
    //point are counted. 

    //spatial
    for(int q=0; q<NodeList[i].fwdLinks+1; q++) {      
      KE += 0.5*((NodeList[i].phi - NodeList[NodeList[i].nn[q]].phi)*
		 (NodeList[i].phi - NodeList[NodeList[i].nn[q]].phi)  );
    }

    
    //temporal (the qth neighbour is always the forward time link.
    KE += 0.5*((NodeList[i].phi - NodeList[NodeList[i].nn[p.q]].phi)*
	       (NodeList[i].phi - NodeList[NodeList[i].nn[p.q]].phi)*
	       (NodeList[i].temporal_weight*NodeList[i].temporal_weight));
  }
  
  return PE + KE;
}


// Prob = EXP[ - E]  AcceptProb = Min[1, exp[ - DeltaE]]

int metropolis_update_phi_AdS(vector<Vertex> &NodeList, vector<int> &s, Param p,
			      double & delta_mag_phi, int iter) {

  uniform_real_distribution<double> unif(0.0,1.0);

  int interior = endNode(p.Levels-1,p)+1;
  int disk     = endNode(p.Levels,p)+1;  
  int s_old;
  int delta_mag = 0;
  int accept = 0;
  int tries = 0;  
  double delta_phi = 1.5;
  double phi_new = 0.0;
  double DeltaE = 0.0;
  double phi_sq;

  delta_mag_phi = 0.0;
  
  for (int i = 0; i < p.latVol; i++) {
    if (s[i] * NodeList[i].phi < 0) {
      printf("ERROR s and phi NOT aligned! (MUP AdS)\n");
    }
    
    DeltaE = 0.0;
    phi_new = NodeList[i].phi + delta_phi * (2.0*unif(rng) - 1.0);
    phi_sq = NodeList[i].phi*NodeList[i].phi;
    
    //PE
    if( i%disk > interior) {
      DeltaE += 0.25*p.lambda*(phi_new*phi_new*phi_new*phi_new - phi_sq*phi_sq);
      DeltaE += 0.5*p.musqr *(phi_new*phi_new - phi_sq);
    }
    DeltaE += 0.5*p.musqr *(phi_new*phi_new  - phi_sq);
      
    
    //KE
    for(int q=0; q<p.q; q++) {
      if(NodeList[i].nn[q] != -1) {
	DeltaE += 0.5*((phi_new - NodeList[NodeList[i].nn[q]].phi)*
		       (phi_new - NodeList[NodeList[i].nn[q]].phi) -
		       (NodeList[i].phi - NodeList[NodeList[i].nn[q]].phi)*
		       (NodeList[i].phi - NodeList[NodeList[i].nn[q]].phi));
      }
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

  if( (iter+1)%100 == 0 ) {
    //printf(" At iter = %d the Acceptance = %d with %d tries and rate %g with delta phi %g\n",
    //iter,accept, tries, (double)accept/(double)tries, delta_phi);
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
  double J = +1.0;
  
  // if the neighbor spin does not already belong to the
  // cluster, then try to add it to the cluster
  for(int q=0; q<p.q+2; q++)
    if(!cluster[NodeList[NodeList[i].nn[q]].pos] && NodeList[NodeList[i].nn[q]].pos != -1) {
      if(s[NodeList[NodeList[i].nn[q]].pos] == clusterSpin && unif(rng) < 1 - exp(2*J*NodeList[i].phi*NodeList[NodeList[i].nn[q]].phi)) {
	cluster_add_AdS(NodeList[NodeList[i].nn[q]].pos, s, clusterSpin, cluster, NodeList, p);
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
