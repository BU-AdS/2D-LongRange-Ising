#ifndef UPDATE_H
#define UPDATE_H

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "graph.h"

//  x = 0,1,.., p.S1-1  and  t = 0, 1, ..., p.Lt
// i =  x + p.S1* r so  x = i % p.S1  and t = i/p.S1 
// 

inline int
xp(int i, Param p)
{
  return (i + 1) % p.S1 + p.S1 * (i / p.S1);
}

inline int
xm(int i,  Param  p)
{
  return (i - 1 + p.S1) % p.S1 + p.S1 * (i / p.S1);
}

/* jump to next t shell */
inline int
tp(int i,  Param  p)
{
  return i % p.S1 + p.S1 * ((i / p.S1 + 1) % p.Lt);
}

/* return to last t shell */
inline int
ttm(int i,  Param  p)   // Some compiler have tm as resurved!
{
  return i % p.S1 + p.S1 * ((i / p.S1 - 1 + p.Lt) % p.Lt);
}

double
action_phi(vector<double> &phi, vector<int> &s, Param p, double & KE,  double & PE)
{
  int i;
  KE = 0.0;
  PE = 0.0;
  //  int p.latVol = p.SurfaceVol;
 
  for (i = 0; i < p.latVol; i++)
    {  
      if (s[i] * phi[i] < 0) printf("ERROR s and phi NOT aligned ! \n");
    };

  //Shaich parameter  p.lambaa = 4 lambda_Schaich
  // p.musqr =  2*musqr_Schaich  + 4
  
  /* PE  terms */
  for (i = 0; i < p.latVol; i++)
    {
      PE += 0.25 * p.lambda* phi[i]*phi[i]*phi[i]*phi[i] - 0.5* p.msqr * phi[i] * phi[i];
    }
  
  /* KE Terms on the sphere AdS and along the cylinder */
  for (i = 0; i <p.latVol; i++)
    {
      KE += 0.5 * (phi[i] - phi[xp(i,p)]) * (phi[i] - phi[xp(i,p)]);
      KE += 0.5 * (phi[i] - phi[tp(i,p)]) * (phi[i] - phi[tp(i,p)]);
    }	  
  
  return PE + KE;
}


/******* Prob = EXP[ - E]  AcceptProb = Min[1, exp[ - DeltaE]]    ********/

int
metropolis_update_phi(vector<double> & phi, vector<int> & s, Param p,
		      double & delta_mag_phi)
{
  int  s_old;
  /* int x; */
   uniform_real_distribution<double> unif;
  int delta_mag = 0;
  int accept = 0;
  int tries = 0;
  static int iter = 0;
  /* static double delta_max = 10.; */
  /* static double delta_min = 0.0; */
  
  static double delta_phi = 1.5;
  
  delta_mag_phi = 0.;
  
  for (int i = 0; i < p.latVol; i++) {
    
    if (s[i] * phi[i] < 0)
      {
	printf("ERROR s and phi NOT aligned ! \n");
      }
    
    double DeltaE = 0.0;
    double phi_new = phi[i] + delta_phi * (2.0*unif(rng) - 1.0);
    
#if 0
    // Put into test.h later
    printf("\n Entering metropolisr for i = %d \n", i);
    double KE = 0.0;
    double PE = 0.0;
    double DeltaKE = 0.0;
    double DeltaPE = 0.0;
    double Energy = action_phi(phi, s,  p, KE,  PE);
    double DeltaEnergy = 0.0;
    DeltaKE = KE;
    DeltaPE = PE; 
    printf(" Energy = %f \n",Energy);
    
    double phi_tmp = phi[i];
    phi[i] = phi_new;
    DeltaEnergy = action_phi(phi, s,  p, KE,  PE) - Energy;
    phi[i] = phi_tmp;
    
    DeltaKE = KE - DeltaKE;
    DeltaPE = PE - DeltaPE;
    printf(" DeltaEnergy = %f DeltaKE = %f DeltaPE = %f  DeltaKE +DeltaPE = %f \n",
	   DeltaEnergy, DeltaKE, DeltaPE, DeltaKE + DeltaPE);
    
#endif
    
    /**********   Potential Term    **************/
    
    DeltaE += 0.25*  p.lambda* ( phi_new*phi_new*phi_new*phi_new - phi[i]*phi[i]*phi[i]*phi[i] );
    DeltaE += - 0.5 * p.msqr* ( phi_new * phi_new  - phi[i]*phi[i] );
    
    /**********   Kinetic  Term  **************/
    
    DeltaE += 0.5*( (phi_new - phi[xp(i, p)]) * (phi_new - phi[xp(i, p)])
		    -  (phi[i] - phi[xp(i, p)]) * (phi[i] - phi[xp(i, p)]) );
    DeltaE += 0.5* ( (phi_new - phi[xm(i, p)]) * (phi_new - phi[xm(i, p)])
		     -  (phi[i] - phi[xm(i, p)]) * (phi[i] - phi[xm(i, p)]) );
    
    DeltaE += 0.5* ( (phi_new - phi[tp(i, p)]) * (phi_new - phi[tp(i, p)])
		     -   (phi[i] - phi[tp(i, p)]) * (phi[i] - phi[tp(i, p)]) );
    DeltaE += 0.5* ( (phi_new - phi[ttm(i, p)]) * (phi_new - phi[ttm(i, p)]) 
		     -   (phi[i] - phi[ttm(i, p)]) * (phi[i] - phi[ttm(i, p)]) );
    
    //    cout<<"i= "<< i <<"  "<< xp(i, p) <<"  " << xm(i, p) <<"  " << tp(i, p) <<"  " << tm(i, p) << endl;
    tries += 1;    
    if(DeltaE < 0.0)
      {
	//  cout<< " Acepted  " << endl;
	s_old = s[i];
	delta_mag_phi += phi_new - phi[i];
	phi[i] = phi_new;
	accept += 1;
	s[i] = (phi_new > 0) ? 1 : -1;
	delta_mag += s[i] - s_old;
      }
    else if ( unif(rng)  < exp(-DeltaE))
      {
	//  cout<< " Acepted  " << endl;1
	s_old = s[i];
	delta_mag_phi += phi_new - phi[i];
	phi[i] = phi_new;
	accept += 1;
	s[i] = (phi_new > 0) ? 1 : -1;
	delta_mag += s[i] - s_old;
      }     
  }                                /* end  for (i = 0; i < p.latVol; i++)  */
  
  /**** TUNING ACCEPTANCE *****/
#if 0
  if (iter < 1000) {
    if ((double) accept / (double) tries < .50) {
      delta_phi -= 0.05;
      /*  delta_max = delta_phi;
	  delta_phi = (delta_max + delta_min) / 2.; */
    } else {
      delta_phi += 0.05;
      /*      delta_min = delta_phi;
	      delta_phi = (delta_max + delta_min) / 2.; */
    }
  }
  
  if(iter < 100 ){
    printf(" At iter = %d the Acceptance  =  %d  with  %d tries and rate  %g \n", iter,accept, tries,
  	   (double)accept/(double)tries);
    //  printf(" delta_phi = %.8f  (delta_min , delta_max) = ( %.8f, %.8f ) \n \n", delta_phi, delta_min, delta_max);
  };
  
  if (iter > 990 && iter < 1000)
    {
      fprintf(stdout, "#RTPHI: delta_phi %.10g accept %d tries %d rate %.3f\n",
	      delta_phi, accept, tries, (double) accept / (double) tries);
    }
#endif
  
  iter++;
  return delta_mag;
}

#if 0
int
swendsen_wang_update(douvle phi, int *s, Param p)
{
  int num_sw_clusters = 0;

  make_sw_clusters_phi(phi,s, p);

  num_sw_clusters = flatten_sw_clusters(p);

  flip_sw_clusters(phi,s, p);

  return num_sw_clusters;
}

int
flatten_sw_clusters(param_t p)
{
  int i, lead, number_cluster = 0, root, trail;

  for (i = 0; i < p->latVol; i++) {
    if (p->cluster[i] == i) {
      number_cluster++;
    } else {
      /* p->cluster[i] = find(i, p) ; */
      /*  inline the collapsing find */
      for (root = i; p->cluster[root] != root; root = p->cluster[root]);

      for (trail = i; trail != root; trail = lead) {
        lead = p->cluster[trail];
        p->cluster[trail] = root;
      }
    }
  }

  return number_cluster;
}

void
flip_sw_clusters_phi(double *phi, int *s, Param p)
{
  int i;
  int *flip = NULL;

  flip = (int *) malloc(p->latVol * sizeof(int));
  if (flip == NULL) {
    fprintf(stderr, "ERROR: flip_sw_clusters: flip allocation failed.\n");
    perror("sw flip");
    exit(-1);
  }
  
  for (i = 0; i < p->latVol; i++) {
    /*
     * find the roots of all cluster trees.
     *  (p->cluster[root] == root) is always true, even if the
     * sw_cluster has not been flattened.
     */
    
    if (p->cluster[i] == i) {
      /* set flip for each cluster */
      flip[i] = (ISING_RANDOM() < 0.5) ? 1 : 0;
#ifdef _DEBUG_
    } else {
      flip[i] = 10;
#endif
    }
#ifdef _DEBUG_
    if (flip[find(i, p)] == 10) {
      fprintf(stderr, "\n ERROR SPIN %d  NOT IN CLUSTER! \n", i);
    }
#endif
    /*
     * For efficiency, assume that the table has been fully flattened,
     * otherwise use:
     *
     *   if (flip[find(i, p)] == 0) {
     */
    if (flip[p.cluster[i]] == 0) {
      s[i] = -s[i];
      phi[i] = - phi[i];
      }
  }      /* end for (i =n 0; i < p.latVol; i++) */

  free(flip);
}


void
make_sw_clusters_phi(double *phi, int *s, param_t p)
{
  int i, j, r, x;
  int i_root, j_root;

  for (i = 0; i < p.latVol; i++) {
    p.cluster[i] = i;                /* everthing is in its own cluster */
  }

  /*
   * THIS FIRST LOOP AT FIXED RADIUS CODE WILL BE DRIVEN BY A NEIGHBORHOOD
   * TABLE FOR EACH SLIDE IN HIGH DIMENISION
   */
  /* printf("NOW WORKING ON INDEVIDUAL r SLICES \n"); */

  /* BEGIN  */
  p.NumClusters = p.latVol;        /* could count mass of cluster as well. */

  for (r = 0; r < p.Lt; r++) {
    int rLx = r * p.g->nv;
    for (x = 0; x < p.g->nv; x++) {
      int b;
      i = x + rLx;
      for (b = p.g->first[x]; b < p.g->first[x + 1]; b++) {
        j = p.g->nn[b] + rLx;
        int nu = b - p.g->first[x] + 1;        /* RCB: change this */
        /* Opposite of George  in force  (j > i) */
        if (i >= j) {
          continue;
        }
        /*  DEBUG */
        /*
         * if (r == 0) {
         *   printf("DEBUG: s[%d] = %d , s[%d] = %d\n", i, s[i], j, s[j]);
         * }
         */
        /* BEGIN: add_bond(s, i, j, p) ; */
        /* bonds++ ; */
        if (s[i] == s[j] && ISING_RANDOM() < 1. -
            exp(-2. * p.beta * p.g->Ksphere[nu][x] * fabs(phi[i] * phi[j]))) { // RCB Fix error in code
          /* p.perq_bond[bonds-1] = 1 ; */
          i_root = find(i, p);
          j_root = find(j, p);
          if (j_root < i_root) {
            p.cluster[i_root] = j_root;
            p.NumClusters--;
          } else if (j_root > i_root) {
            p.cluster[j_root] = i_root;
            p.NumClusters--;
          }
        }
        /* END: add_bond(s, i, j, p) ; */
      }
    }
  }
  /* END  */

  /*
   * SECOND PASS CONNECT THE FIXED R SLICES
   */

  /* couple next r slice to last r slice. */

  for (j = p.g->nv; j < p.latVol; j++) {
    i = j - p.g->nv;
    x = i % p.g->nv;
    /* BEGIN: add_bond(s, i, j, p) ; */
    /* bonds++ ; */
    if (s[i] == s[j] && ISING_RANDOM() <1. -
          exp(-2. * p.beta * p.g->Karea[x] * fabs(phi[i] * phi[j]))) {
      /* p.perq_bond[bonds-1] = 1 ; */
      i_root = find(i, p);
      j_root = find(j, p);
      if (j_root < i_root) {
        p.cluster[i_root] = j_root;
        p.NumClusters--;
      } else if (j_root > i_root) {
        p.cluster[j_root] = i_root;
        p.NumClusters--;
      }
    }
    /* END: add_bond(s, i, j, p) ; */
  }

  // PERIODIC:
    for (i = 0; i < p.g->nv; i++) {
      j = i + p.latVol - p.g->nv;
       x = i % p.g->nv;
      /* BEGIN: add_bond(s, i, j, p) ; */
     /* bonds++ ; */
      if (s[i] == s[j] && ISING_RANDOM() < 1. -
          exp(-2. * p.beta * p.g->Karea[x] * fabs(phi[i] * phi[j]))) {
        /* p.perq_bond[bonds-1] = 1 ; */
        i_root = find(i, p);
        j_root = find(j, p);
        if (j_root < i_root) {
          p.cluster[i_root] = j_root;
          p.NumClusters--;
        } else if (j_root > i_root) {
          p.cluster[j_root] = i_root;
          p.NumClusters--;
        }
      }
      /* END: add_bond(s, i, j, p) ; */
    }
 
  return;
}

#endif

#endif   //End of UPDATE_H

