#ifndef CG_H
#define CG_H
#include <complex>
#include <cstring>
//#include <util.h>

using namespace std;

int Mphi(vector<double> &phi, const vector<double> phi0,
	 vector<Vertex> NodeList, Param P) {

  int Levels = P.Levels;
  int q = P.q;
  int T = P.t;
  double msqr = P.msqr;
  bool bc = P.bc;
  int InternalNodes = endNode(Levels-1,P)+1;
  int TotNumber = endNode(Levels,P)+1;
  int offset = TotNumber;
  int T_offset = 0;
  T == 1 ? T_offset = 0 : T_offset = 2;
  
  for(int t = 0; t<T; t++) {

    //loop over interior nodes on all disks    
    for(int i = t*offset; i < t*offset + InternalNodes; i++) {
      
      //mass term
      phi[i] = msqr * phi0[i];    
      
      for(int mu = 0; mu < q+T_offset; mu++) {
	phi[i] += (phi0[i] - phi0[NodeList[i].nn[mu]]);
      }
    }
  
    // Dirichlet or Neuman at Exterior Nodes.  
    for(int i = t*offset + InternalNodes; i < t*offset + TotNumber; i++){
      
      //cout<<"Exterior i="<<i<<" t="<<t<<endl;
      //mass term
      phi[i] = msqr * phi0[i];
      //if(bc == true) phi[i] += q * phi0[i];
      
      //the zeroth link is always on the same level.
      phi[i] += phi0[i] - phi0[NodeList[i].nn[0]];
      
      //The q-1 (and q-2) link(s) go back one level,
      //the q-3 (or q-2) link is on the same level. These
      //links are always computed in the same way.
      for(int mu = q-1; mu > (NodeList[i].fwdLinks); mu--) {
	phi[i] += phi0[i] - phi0[NodeList[i].nn[mu]];
      }      

      //We use the member data fwdLinks to apply the boundary
      //condition. For Dirichlet, the field value is 0. For
      //Neumann, the derivative is zero.
      for(int mu = NodeList[i].fwdLinks; mu > 0; mu--) {
	if(bc == true) {
	  //Apply Dirchlet BC.
	  phi[i] += phi0[i];
	} else {
	  //Apply Neumann
	  phi[i] += 0.0;
	}
      }
      
      //Temporal links at exterior
      if(T>1) {
	for(int mu = q; mu < q+2; mu++) {
	  phi[i] += (phi0[i] - phi0[NodeList[i].nn[mu]]);
	}
      }    
    }
  }
  return 0;  
}


double Minv_phi(vector<double> &phi, const vector<double> phi0,
		const vector<double> b, const vector<Vertex> NodeList,
		Param P)
{
  // CG solutions to Mphi = b 
  //  see http://en.wikipedia.org/wiki/Conjugate_gradient_method
  int Levels = P.Levels;
  int diskN = endNode(Levels,P) + 1;
  int N = P.t*diskN;
  
  vector<double> res(N,0.0), resNew(N,0.0),  pvec(N,0.0), Mpvec(N,0.0), pvec_tmp(N,0.0);
  double alpha, beta, denom;
  double rsq = 0, rsqNew = 0, bsqrt = 0, truersq = 0.0;
  int  i;
  
  for(i = 0; i<N; i++){
    res[i] = b[i];
    pvec[i] = res[i];
    bsqrt += b[i]*b[i];
  }
  bsqrt = sqrt(bsqrt);
  
  int maxIter = P.MaxIter;
  double resEpsilon = P.tol;
  // iterate till convergence
  rsqNew = 100.0;
  int k = 0;
  
  while( (k<maxIter)&&(sqrt(rsqNew) > resEpsilon*bsqrt) ){
    
    k++;
    rsq = 0;
    for (int i = 0; i < N; i++) rsq += res[i]*res[i];
    
    cout << endl << "In CG at iteration = "<<k <<" Residue Squared  = " << rsq << endl;

    //Mat-Vec operation
    Mphi(Mpvec, pvec, NodeList, P);  

    denom = 0;
    for(i=0; i< N; i++) denom += pvec[i]*Mpvec[i];
    alpha = rsq/denom;
    
    for(i=0; i < N; i++) phi[i] +=  alpha * pvec[i];
    for(i=0; i < N; i++) resNew[i] = res[i]- alpha*Mpvec[i];
    
    // Exit if new residual is small enough
    rsqNew = 0;
    for (i = 0; i < N; i++) rsqNew += resNew[i]*resNew[i];
    
    // Update vec using new residual
    beta = rsqNew / rsq;
    for (i = 0; i < N; i++) {
      pvec[i] = resNew[i] + beta * pvec[i];
      res[i] = resNew[i];
    }
  }
  
  if(k == maxIter) {
    printf("CG: Failed to converge iter = %d, rsq = %e\n", k, rsq); 
    //  Failed convergence 
  }
  
  Mphi(Mpvec, phi, NodeList, P);  
  for(int i=0; i < N ; i++) truersq += (Mpvec[i] - b[i])*(Mpvec[i] - b[i]);
  
  printf("# CG: Converged iter = %d, rsq = %e, truesq = %e\n",k,rsq,truersq);

  return truersq; // Convergence 
}

#endif
