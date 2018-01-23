#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>

#include "util.h"
#include "hyp_util.h"
#include "data_proc.h"

using namespace std;

//Overloaded version to handle AdS lattices
void correlators(double **corr, double **corr_ave, int corr_norm,
		 vector<Vertex> NodeList, double avePhi, Param p) {

  double *phi = (double*)malloc(p.S1*p.Lt*sizeof(double));
  int offset = endNode(p.Levels-1,p)+1;
  int disk   = p.AdSVol;  
  for(int i=0; i<p.S1*p.Lt; i++) {    
    phi[i] = NodeList[disk*(i/p.S1) + offset + i%p.S1].phi;
  }

  correlators(corr, corr_ave, corr_norm, phi, avePhi, p);
  
}


void correlators(double **corr, double **corr_ave, int corr_norm,
		 double *phi, double avePhi, Param p) {
  
  int s_idx = 0;
  int t_idx = 0;

  int s_size = p.S1;
  int t_size = p.Lt;

  //Each value of \delta t and \delta \theta
  //must be normalised seperately.
  int norm[s_size/2][t_size/2];
  for(int i=0; i<s_size/2; i++)
    for(int j=0; j<t_size/2; j++) {
      norm[i][j] = 0;
      corr[i][j] = 0.0;
    }
  
  //loop over sink/source *theta*
  for(int is=0; is<s_size; is++)
    for(int js=0; js<s_size; js++) {
      
      s_idx = abs(is-js);
      if(s_idx != s_size/2) {
	if(s_idx >= s_size/2) s_idx = s_size - s_idx - 1;
      
	//loop over sink/source *temporal*
	for(int il=0; il<t_size; il++) 
	  for(int jl=0; jl<t_size; jl++) {
	    
	    t_idx = abs(il-jl);
	    if(t_idx != t_size/2) {
	      if(t_idx >= t_size/2) t_idx = t_size - t_idx - 1;
	      
	      corr[s_idx][t_idx] += ((phi[is + p.S1*il] - avePhi) *
				     (phi[js + p.S1*jl] - avePhi));
	      norm[s_idx][t_idx]++;
	    }
	  }
      }
    }
  
  //Normalise, add to running average
  for(int i=0; i<s_size/2; i++)
    for(int j=0; j<t_size/2; j++) {
      corr[i][j] /= norm[i][j];
      corr_ave[i][j] += corr[i][j];
    }

  //Corr dump to stdout and collect running average.
  if(p.verbosity) cout<<setprecision(4);      
  if(p.verbosity) cout<<"Corr Sample: "<<endl;
  for(int i=0; i<s_size/2; i++) {
    if(p.verbosity) cout<<"theta "<<2*M_PI*i/(double)(s_size)<<":";
    for(int j=0; j<t_size/2; j++) {
      if(p.verbosity) cout<<" "<<corr_ave[i][j]/corr_norm;
      corr[i][j] = corr_ave[i][j]/corr_norm;
    }
    if(p.verbosity) cout<<endl;
  }

  if(!p.verbosity) {
    cout<<setprecision(4);      
    cout<<"Corr Sample:";
    for(int j=0; j<t_size/2; j++) {
      cout<<" "<<corr[0][j];
    }
    cout<<endl;
  }
  
  //corr now contains the current, normalised, running average correlation
  //matrix.  
}

void autocorrelation(double *PhiAb_arr, double avePhiAbs, int meas) {
  
}

void jackknife(double *array, int block, int length) {

}
