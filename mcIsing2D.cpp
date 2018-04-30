#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>
#include <omp.h>
#include <chrono>
#include <algorithm>

#include "util.h"
#include "mcIsing2D.h"

using namespace std;

extern int seed;
extern mt19937 rng;
extern uniform_real_distribution<double> unif;
extern int CACHE_LINE_SIZE;

int **addedI;
vector<int> toCheckI;
int sitesToCheckI;
int clusterIterI;

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

//----------------------//
// Short Range Routines //
//----------------------//

// declare variables to implement Cluster algorithms
//Declare variables to implement Cluster algorithms
int ising_wc_ave = 0;
int ising_wc_size = 0;
int ising_wc_calls = 0;
int ising_wc_poss = 0;
int ising_wc_t_size = 0;
int ising_wc_s_size = 0;

int ising_sw_ave = 0;
int ising_sw_size = 0;
int ising_sw_calls = 0;
int ising_sw_t_size = 0;
int ising_sw_s_size = 0;

int ising_accept = 0;
int ising_tries  = 0;

void metropolisUpdateISR(int *s, Param p, int iter) {

  double DeltaE = 0.0;
  double extern_h = p.h;
  double J = p.J;
  
  for (int i = 0; i < p.surfaceVol; i++) {

    //External field
    DeltaE += -extern_h*s[i];
    
    //Interaction
    DeltaE += -J*(s[xp(i, p)] + s[xm(i, p)] +
		  s[tp(i, p)] + s[ttm(i, p)]);
    
    ising_tries++;
    
    if(DeltaE < 0.0) {
      //  cout<< " Acepted  " << endl;
      ising_accept += 1;
      s[i] *= -1;
    }
    else if ( unif(rng)  < exp(-DeltaE)) {
      //  cout<< " Acepted  " << endl;
      ising_accept += 1;
      s[i] *= -1;
    }     
  }// end loop over lattice volume

  if( iter < p.n_therm && iter%p.n_skip == 0 ) {
    cout<<"At iter "<<iter<<" the Metro acceptance rate is "<<(double)ising_accept/(double)ising_tries<<endl;
  }
  
}

void swendsenWangUpdateISR(int *s, Param p, int iter) {
  

}

void swendsenWangClusterAddISR(int i, int *s, int cSpin, int clusterNum, 
			       int *clusterDef, Param p) {
  
}

void wolffUpdateISR(int *s, Param p, int iter) {

  ising_wc_calls++;
  
  bool *cluster = new bool[p.surfaceVol];
  for (int i = 0; i < p.surfaceVol; i++)
    cluster[i] = false;

  // choose a random spin and grow a cluster
  int i = int(unif(rng) * p.surfaceVol);
  int cSpin = s[i];
  
  //This function is recursive and will call itself
  //until all four attempts in the lattice directions
  // (+x, -x, +t, -t) have failed to ncrease the cluster.
  ising_wc_size = 1;
  double prob = 1 - exp(-2*p.J);
  wolffClusterAddISR(i, s, cSpin, prob, cluster, p);

  ising_wc_ave += ising_wc_size;

  if(iter%p.n_skip == 0) {
    setprecision(4);
    cout<<"Average (CPU) cluster size at iter "<<iter<<" = "<<ising_wc_ave<<"/"<<ising_wc_calls<<" = "<<(1.0*ising_wc_ave)/ising_wc_calls<<" = "<<(100.0*ising_wc_ave)/(ising_wc_calls*p.surfaceVol)<<"%"<<endl;
  }
  
  delete[] cluster;
  
}

void wolffClusterAddISR(int i, int *s, int cSpin, double prob,
			bool *cluster, Param p) {

  // The site belongs to the cluster, so flip it.
  cluster[i] = true;
  s[i] *= -1;

  // - If the (aligned) neighbour spin does not already belong to the
  // cluster, then try to add it to the cluster.
  // - If the site has already been added, then we may skip the test.

  //Forward in T
  if(!cluster[ tp(i,p) ] && s[tp(i,p)] == cSpin) {
    if(unif(rng) < prob) {
      ising_wc_size++;
      ising_wc_t_size++;
      //cout<<"->tp";
      wolffClusterAddISR(tp(i,p), s, cSpin, prob, cluster, p);
    }
  }

  //Forward in X
  if(!cluster[ xp(i,p) ] && s[xp(i,p)] == cSpin) {
    if(unif(rng) < prob) {
      ising_wc_size++;
      ising_wc_s_size++;
      //cout<<"->xp";
      wolffClusterAddISR(xp(i,p), s, cSpin, prob, cluster, p);
    }
  }

  
  //Backard in T 
  if(!cluster[ ttm(i,p) ] && s[ttm(i,p)] == cSpin) {  
    if(unif(rng) < prob) {
      ising_wc_size++;
      ising_wc_t_size++;
      //cout<<"->tm";
      wolffClusterAddISR(ttm(i,p), s, cSpin, prob, cluster, p);
    }
  }

  //Backward in X
  if(!cluster[ xm(i,p) ] && s[xm(i,p)] == cSpin) {
    if (unif(rng) < prob) {
      ising_wc_size++;
      ising_wc_s_size++;
      //cout<<"->xm";
      wolffClusterAddISR(xm(i,p), s, cSpin, prob, cluster, p);
    }
  }   
}


//---------------------//
// Long Range Routines //
//---------------------//

void init_connectivityI(Param p) {
  
  addedI = (int**)malloc(p.surfaceVol*sizeof(int*));  
  for(int i=0; i<p.surfaceVol; i++) {
    addedI[i] = (int*)malloc((p.surfaceVol+1)*sizeof(int));
#ifdef USE_OMP
#pragma omp parallel for
#endif
    for(int j=0; j<p.surfaceVol+1; j++) {
      addedI[i][j] = -1;
    }
  }
}



void Ising2D::metropolisUpdateILR(Param p, int iter) {
  
  
  
}

void Ising2D::swendsenWangUpdateILR(Param p, int iter) {
  

}

void swendsenWangClusterAddILR(int i, int *s, int cSpin, int clusterNum, 
			       int *clusterDef, vector<int> Rcluster, 
			       double *LR_couplings, Param p) {
  
}

void clusterPossibleILR(int i, int *s, int cSpin,
		       bool *Pcluster, vector<int> &Rcluster,
		       Param p) {


}

void wolffUpdateILR(int *s, Param p, double *isingProb, int iter) {

  
  // #ifdef USE_OMP
  //   init_connectivityI(p);
  // #endif

  ising_wc_calls++;

  if(iter == p.n_therm) {
    cout<<"Resetting Cluster Stats."<<endl;    
    ising_wc_calls = 1;
    ising_wc_ave = 0;
  }
  
  //Choose a random spin.
  int i = int(unif(rng) * p.surfaceVol);
  int cSpin = s[i];
  
  // The site belongs to the cluster, so flip it.
  ising_wc_size = 1;
  s[i] *= -1;
  
  //This function is recursive and will call itself
  //until all attempts to increase the cluster size
  //have failed.
  wolffClusterAddILR(i, s, cSpin, isingProb, p);
  
  ising_wc_ave += ising_wc_size;

  if(iter%p.n_skip == 0) {
    setprecision(4);
    cout<<"Average (CPU) cluster size at iter "<<iter<<" = "<<ising_wc_ave<<"/"<<ising_wc_calls<<" = "<<(1.0*ising_wc_ave)/ising_wc_calls<<" = "<<(100.0*ising_wc_ave)/(ising_wc_calls*p.surfaceVol)<<"%"<<endl;
  }

// #ifdef USE_OMP
//   for(int a=0; a<p.surfaceVol; a++) free(addedI[a]);
//   free(addedI);
//   toCheckI.clear();
// #endif
  
}

void wolffClusterAddILR(int i, int *s, int cSpin, double *isingProb, Param p) {
  
  int S1 = p.S1;
  int Lt = p.Lt;
  int x_len = S1/2 + 1;
  int vol = S1*Lt;

  int t1,x1,t2,x2,dt,dx;
  t1 = i / p.S1;
  x1 = i % p.S1;
  
  double prob = 0.0;
  double rand = 0.0;
  //We now loop over the possible lattice sites, adding sites
  //(creating bonds) with the specified LR probablity.
  for(int j=0; j<vol; j++) {
    if(s[j] == cSpin) {
      
      //Index divided by circumference, using the int floor feature/bug,
      //gives the timeslice index.
      t2 = j / p.S1;
      dt = abs(t2-t1) > Lt/2 ? Lt - abs(t2-t1) : abs(t2-t1);
      
      //The index modulo the circumference gives the spatial index.
      x2 = j % p.S1;            
      dx = abs(x2-x1) > S1/2 ? S1 - abs(x2-x1) : abs(x2-x1);      
      
      //prob = 1 - exp(-2*J*LR_couplings[dx + dt*x_len]);
      prob = isingProb[dx + dt*x_len];
      rand = unif(rng);
      if(rand < prob) {
	ising_wc_size++;
	// The site belongs to the cluster, so flip it.
	s[j] *= -1;
	wolffClusterAddILR(j, s, cSpin, isingProb, p);
      }
    }
  }
}

//This routine is a prototype for the GPU parallel routine. Is has not yet been
//optimised.

void wolffUpdateILRProto(int *s, Param p, double *isingProb, int iter) {
  
  ising_wc_calls++;

  if(iter == p.n_therm) {
    cout<<"Resetting Cluster Stats."<<endl;    
    ising_wc_calls = 1;
    ising_wc_ave = 0;
  }

  int S1 = p.S1;
  int Lt = p.Lt;
  int vol = S1*Lt;
  
  int *spinStack = (int*)malloc(vol*sizeof(int));
  int *spinStackLC = (int*)malloc(vol*sizeof(int));
  int stackSize = 0;
  
  //Choose a random spin.
  int i = int(unif(rng) * p.surfaceVol);
  int cSpin = s[i];
  
  // The site belongs to the cluster. Flip it, add to the stack.
  ising_wc_size = 1;
  s[i] *= -1;
  spinStack[0] = i;
  stackSize++;
  while(stackSize > 0) {
    //This function is NOT recursive. It performs a single sweep of the lattice
    //to ascertain which sites (if any) will be added.
    wolffClusterAddILRProto(spinStack, spinStackLC, s, stackSize, cSpin, isingProb, p);
  }
  
  ising_wc_ave += ising_wc_size;
  
  if(iter%p.n_skip == 0) {
    setprecision(4);
    cout<<"Average (CPU) cluster size at iter "<<iter<<" = "<<ising_wc_ave<<"/"<<ising_wc_calls<<" = "<<(1.0*ising_wc_ave)/ising_wc_calls<<" = "<<(100.0*ising_wc_ave)/(ising_wc_calls*p.surfaceVol)<<"%"<<endl;
  }  
  free(spinStack);
  free(spinStackLC);
}

int wolffClusterAddILRProto(int *spinStack, int *spinStackLC, int *s, int &stackSize,
			    int cSpin, double *isingProb, Param p){
  
  //Here we construct the cumulative probability function.
  //We have n sites in the spin stack. The probability of a
  //candidate site k being added from site i=1...n is p(add i). The
  //probability of k being added is:
  //
  //P(add k) = 1 - p(!add 1)p(!add 2)...p(!add n).
  //
  //We know that p(!add i) = 1 - ( 1 - exp(-2*J_{ik})) 
  //                       = exp(-2*J_{ik})
  //hence
  //
  //P(add k) = 1 - exp(-2*\sum_i J_{ik})
  //
  //We add site k if a random number is less than
  //P(add k).
  
  int t1,x1,t2,x2,dt,dx;
  int S1 = p.S1;
  int Lt = p.Lt;
  int x_len = S1/2 + 1;
  int vol = S1*Lt;
  
  double prob = 0.0;
  double rand = 0.0;

  int stackSizeLC = stackSize;
  stackSize = 0;
  
  //We now loop over the possible lattice sites, adding sites
  //(creating bonds) with the specified LR probablity.
  
  for(int j=0; j<vol; j++) {
    if(s[j] == cSpin) {
      
      //Index divided by circumference, using the int floor feature/bug,
      //gives the timeslice index.
      t2 = j / S1;
      
      //The index modulo the circumference gives the spatial index.
      x2 = j % S1;            
      
      //Now we may loop over the spinStack to calculate P(add k)        
      prob = 1.0;
      for(int a=0; a<stackSizeLC; a++) {

	t1 = spinStack[a]/S1;
	x1 = spinStack[a]%S1;	
	dt = abs(t2-t1) > Lt/2 ? Lt - abs(t2-t1) : abs(t2-t1);
	dx = abs(x2-x1) > S1/2 ? S1 - abs(x2-x1) : abs(x2-x1);

	prob *= (1-isingProb[dx + x_len*dt]);
      }
      
      prob = 1 - prob;
      rand = unif(rng);
      
      if(rand < prob) {
	// The site belongs to the cluster, so flip it, and add to the stack.	
	ising_wc_size++;
	s[j] *= -1;
	spinStackLC[stackSize] = j;
	stackSize++;
      }
    }
  }
  for(int a=0; a<stackSize; a++) spinStack[a] = spinStackLC[a];
  return stackSize;
}

double energyISR(int *s, Param p, double &KE) {  

  KE = 0.0;
  
  int S1 = p.S1;
  int Lt = p.Lt;
  int vol = S1*Lt;  
  int s_lc = 1;
  
  for (int i = 0; i < vol; i++) {    
    s_lc = s[i];    
    KE += s_lc*(s[xp(i,p)] + s[xm(i,p)] + s[tp(i,p)] + s[ttm(i,p)]);  
  }
  return -0.5*p.J*KE;
}

double energyILR(int *s, Param p, double *LR_couplings, double &KE) {  

  KE = 0.0;
  int S1 = p.S1;
  int Lt = p.Lt;
  int vol = S1*Lt;  
  int x_len = S1/2 + 1;
  int s_lc = 1;
  
  for (int i = 0; i < vol; i++) {    
    s_lc = s[i];

    int t1,x1;
    t1 = i / S1;
    x1 = i % S1;

    //KE terms
#ifdef USE_OMP
    
    int chunk = CACHE_LINE_SIZE/sizeof(double);
    
#pragma omp parallel 
    {
      double local = 0.0;
      double val = 0.0;
      int t2,dt,x2,dx;
#pragma omp for nowait 
      for(int j=0; j<vol; j+=chunk) {
	for(int k=0; k<chunk; k++) { 

	  if( (j+k) != i ) {
	    //Index divided by circumference, using the int floor feature/bug,
	    //gives the timeslice index.
	    t2 = (j+k) / S1;
	    dt = abs(t2-t1) > Lt/2 ? Lt - abs(t2-t1) : abs(t2-t1);
	    
	    //The index modulo the circumference gives the spatial index.
	    x2 = (j+k) % S1;            
	    dx = abs(x2-x1) > S1/2 ? S1-abs(x2-x1) : abs(x2-x1);      
	    
	    val = s_lc*s[j]*LR_couplings[dx+dt*x_len];
	    local += val;
	  }
	}
      }
#pragma omp atomic 
      KE += local;
    }    
#else

    int t2,dt,x2,dx;
    for(int j=0; j<p.surfaceVol; j++) {
      
      if( j != i ) {
	//Index divided by circumference, using the int floor feature/bug,
	//gives the timeslice index.
	t2 = j / p.S1;
	dt = abs(t2-t1) > p.Lt/2 ? p.Lt - abs(t2-t1) : abs(t2-t1);
	
	//The index modulo the circumference gives the spatial index.
	x2 = j % p.S1;            
	dx = abs(x2-x1) > p.S1/2 ? p.S1-abs(x2-x1) : abs(x2-x1);      

	KE += s_lc*s[j]*LR_couplings[dx+dt*x_len];
	
      }
    }
#endif
  }
  
  return -0.5*p.J*KE;
}

