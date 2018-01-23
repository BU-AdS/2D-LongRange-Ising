#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>

#include "cg.h" 
#include "util.h"

using namespace std;

void Mphi(double *phi, const double *phi0, 
	  const vector<Vertex> NodeList, Param P) { 

  int Levels = P.Levels; 
  int q = P.q; 
  int T = P.t; 
  double msqr = P.msqr; 
  double C_msqr = P.C_msqr; 
  
  bool bc = P.bc; 
  int TotNumber = endNode(Levels,P)+1; 
  int T_offset = 0; 
  T == 1 ? T_offset = 0 : T_offset = 2; 
  
  //loop over nodes on all disks     
  for(int i = 0; i < T*TotNumber; i++) { 
    if(NodeList[i].pos != -1) { 
      
      //mass term 
      phi[i] = C_msqr*msqr*phi0[i];     

      //links 
      for(int mu = 0; mu < q; mu++) { 

	//interior nodes 
	if(NodeList[i].nn[mu] != -1) phi[i] += (phi0[i] - phi0[NodeList[i].nn[mu]]); 

	//boundary nodes 
	else { 
	  if(bc == true) { 
	    //Apply Dirchlet BC. 
	    phi[i] += phi0[i]; 
	  } else { 
	    //Apply Neumann 
	    phi[i] += 0.0; 
	  }  
	} 
      } 
      
      //Temporal links 
      for(int mu = q; mu < q+T_offset; mu++) { 
	phi[i] += (phi0[i] - phi0[NodeList[i].nn[mu]])*(NodeList[i].temporal_weight); 
      }     
    } 
  } 
} 

double Minv_phi(double *phi, double *b, 
	       const vector<Vertex> NodeList, Param P) { 
  // CG solutions to Mphi = b  
  //  see http://en.wikipedia.org/wiki/Conjugate_gradient_method 
  int Levels = P.Levels; 
  int diskN = endNode(Levels,P) + 1; 
  int TotNumber = P.t*diskN; 
  
  double *res, *resNew, *pvec, *Mpvec, *pvec_tmp;   
  res      = new double[TotNumber]; 
  resNew   = new double[TotNumber]; 
  pvec     = new double[TotNumber]; 
  Mpvec    = new double[TotNumber]; 
  pvec_tmp = new double[TotNumber]; 

  for(int i=0; i<TotNumber; i++) { 
    res[i]      = 0.0; 
    resNew[i]   = 0.0; 
    pvec[i]     = 0.0;   
    Mpvec[i]    = 0.0; 
    pvec_tmp[i] = 0.0; 
  } 
  
  double alpha = 0.0, beta = 0.0, denom = 0.0; 
  double rsq = 0.0, rsqNew = 0.0, bsqrt = 0.0, truersq = 0.0; 
  int  i; 
  
  for(i = 0; i<TotNumber; i++){ 
    res[i] = b[i]; 
    pvec[i] = res[i]; 
    bsqrt += b[i]*b[i]; 
  } 
  bsqrt = sqrt(bsqrt); 
  
  int maxIter = P.MaxIter; 
  double resEpsilon = P.tol; 
  // iterate until convergence 
  rsqNew = 100.0; 
  int k = 0; 
  
  while( (k<maxIter)&&(sqrt(rsqNew) > resEpsilon*bsqrt) ){ 
    
    k++; 
    rsq = 0; 
    for (int i = 0; i < TotNumber; i++) { 
      if(NodeList[i].pos != -1) { 
	rsq += res[i]*res[i]; 
      } 
    } 
    //if(k%10 == 0) cout<< "In CG at iteration = "<< k << " Montonicity "<<(beta  < 1.0 ? "True" : "False")<<" Residue Squared  = " << rsq <<endl;

    //Mat-Vec operation
    Mphi(Mpvec, pvec, NodeList, P);  
    
    denom = 0;
    for(i=0; i<TotNumber; i++) {
      if(NodeList[i].pos != -1) {
	denom += pvec[i]*Mpvec[i];
      }
    }
    
    //cout<<"Denom "<<k<<" = "<<denom<<endl;
    //exit(0);
    
    alpha = rsq/denom;
    
    for(i=0; i < TotNumber; i++) {
      if(NodeList[i].pos != -1) {
	phi[i] +=  alpha * pvec[i];
      }
    }
    for(i=0; i < TotNumber; i++) {
      if(NodeList[i].pos != -1) {
	resNew[i] = res[i] - alpha*Mpvec[i];
      }
    }
    
    // Exit if new residual is small enough
    rsqNew = 0;
    for (i = 0; i < TotNumber; i++) {
      if(NodeList[i].pos != -1) {
	rsqNew += resNew[i]*resNew[i];
      }
    }
    // Update vec using new residual
    beta = rsqNew / rsq;
    for (i = 0; i < TotNumber; i++) {
      if(NodeList[i].pos != -1) {      
	pvec[i] = resNew[i] + beta * pvec[i];
	res[i] = resNew[i];
      }
    }
  }
  
  if(k == maxIter) {
    printf("CG: Failed to converge iter = %d, rsq = %e\n", k, (double)rsq); 
    //  Failed convergence 
  }
  
  Mphi(Mpvec, phi, NodeList, P);  
  for(int i=0; i < TotNumber ; i++)
    if(NodeList[i].pos != -1) {
      truersq += (Mpvec[i] - b[i])*(Mpvec[i] - b[i]);
    }
  
  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,(double)rsq,(double)truersq);
  
  return truersq; // Convergence 
}

//This function will calculate the M(x,x)^{-1} elements on the qth
//portion of the graph. 
void oneLoopCorrection(std::vector<Vertex> &NodeList, Param& p){
  
  int n_shift = p.n_shift;
  int latVol  = p.latVol;
  int sources = (endNode(p.Levels,p) - endNode(p.Levels-1,p)) / p.q ;  

  cout<<"Solving for "<<sources<<" sources."<<endl;
  
  double** phi = new double*[n_shift];
  for(int i=0; i<n_shift; i++) {
    phi[i] = new double[latVol];
    for(int j=0; j<latVol; j++) phi[i][j] = 0.0;
  }
  double* b = new double[latVol];  
  for(int s=0; s<sources; s++) { 
    for(int i=0; i<latVol; i++) {
      b[i] = 0.0;
      phi[0][i] = 0.0;
    }
    b[endNode(p.Levels-1,p) + 1 + s] = 1.0;    
    Minv_phi(phi[0], b, NodeList, p);
    cout<<phi[0][endNode(p.Levels-1,p) + 1 + s]<<endl;
    NodeList[endNode(p.Levels-1,p) + 1 + s].oneLoopCorr = phi[0][endNode(p.Levels-1,p) + 1 + s];
  }

  //Propagate the one loop correction through the lattice.
  int pos = endNode(p.Levels-1,p) + 1;
  for(int i=0; i<p.t; i++) 
    for(int q=0; q<p.q; q++)     
      for(int s=0; s<sources; s++) 
	NodeList[pos + s + q*sources + i*p.AdSVol].oneLoopCorr = NodeList[pos + s].oneLoopCorr;
  
  
  //Minv_phi_ms(phi, b, NodeList, p);
  //DataDump(NodeList, phi[i], p, p.Levels, p.t/2, i);

  //Mphi_ev(NodeList, p);
  
  for(int i=0; i<n_shift; i++) delete[] phi[i];
  delete[] phi;
  delete[] b;
}
