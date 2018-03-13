#ifdef USE_GPU
#include <cuda.h>
#include <stdio.h>

voud GPU_wolffUpdateLR(double *phi, int *s, Param p,
		       double *LR_couplings,
		       double &delta_mag_phi, int iter) {
  
  //The GPU Wolff needs a few extra arrays
  
  double phi_lc = phi[i];
  int newSites = 0;
  int t1,x1;
  t1 = i / p.S1;
  x1 = i % p.S1;
  
  //We now loop over the possible lattice sites, adding sites
  //(creating bonds) with the specified LR probablity.
  
  double prob = 0.0;
  double rand = 0.0;
  double added_local[lc_sze];
  for(int k=0; k<lc_sze; k++) {
    added_local[k] = -1;
  }
  int t2,dt,x2,dx;
  
  for(int j=0; j<p.surfaceVol; j+=lc_sze) {
    for(int k=0; k<lc_sze; k++) {
      int idx = j+k;
      if(s[idx] == cSpin) {
	
	//Index divided by circumference, using the int floor feature/bug,
	//gives the timeslice index.
	t2 = (idx) / p.S1;
	dt = abs(t2-t1) > p.Lt/2 ? p.Lt - abs(t2-t1) : abs(t2-t1);
	
	//The index modulo the circumference gives the spatial index.
	x2 = (idx) % p.S1;            
	dx = abs(x2-x1) > p.S1/2 ? p.S1 - abs(x2-x1) : abs(x2-x1);      
	
	prob = 1 - exp(2*phi_lc*phi[idx]*LR_couplings[dt+dx*LC_arr_len]);
	rand = unif(rng);
	if(rand < prob) {
	  sqnl_wc_size++;
	  added_local[k] = 1;
	  ++newSites_local;
	}
      }
    }
  }
  
  if(newSites > 0) {
    //'added' now contains all of the spins on the wave front. Each
    //of these new sites must be explored, but no new exploration 
    //need test these sites. First, sort the added array so that all
    //the hits are at the beginning, order is unimportant.
    int pos = 0;
    int l = 0;
    for(l=0; l<p.surfaceVol; l++) {
      if(added[i][l] > 0) {
	added[i][pos] = l;
	pos++;
      }
    }

    // These sites belongs to the cluster, so flip them.
    for(int k=0; k<newSites; k++) {
      s[added[i][k]] *= -1;
      phi[added[i][k]] *= -1;
    }    
    for(int k=0; k<newSites; k++)
      wolffClusterAddLR(added[i][k], s, cSpin, LR_couplings, phi, p);
  }  
}














#endif
