#ifndef CG_H
#define CG_H
#include <complex>
#include <cstring>
//#include <util.h>

using namespace std;

int Mphi(Float *phi, const Float *phi0,
	 const vector<Vertex> NodeList, Param P) {

  int Levels = P.Levels;
  int q = P.q;
  int T = P.t;
  Float msqr = P.msqr;
  Float C_msqr = P.C_msqr;
  bool bc = P.bc;
  int TotNumber = endNode(Levels,P)+1;
  int T_offset = 0;
  T == 1 ? T_offset = 0 : T_offset = 2;
  
  //loop over nodes on all disks    
  for(int i = 0; i < T*TotNumber; i++) {
    
    if(NodeList[i].pos != -1) {
      //cout<<NodeList[i].pos<<endl;
      //mass term
      phi[i] = C_msqr*msqr * phi0[i];    
      
      //Spatial links
      for(int mu = 0; mu < q; mu++) {
	if(NodeList[i].nn[mu] != -1) phi[i] += (phi0[i] - phi0[NodeList[i].nn[mu]]);
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
	phi[i] += (phi0[i] - phi0[NodeList[i].nn[mu]])/(NodeList[i].temporal_weight);
      }
    }
  }
  
  /*
    //Dirichlet or Neuman at Exterior Nodes.  
    for(int i = t*offset + InternalNodes; i < t*offset + TotNumber; i++){
      
      if(NodeList[i].pos != -1) {
      
	//cout<<"Exterior i="<<i<<" t="<<t<<endl;
	//mass term
	phi[i] = C_msqr*msqr * phi0[i];
	
	//Spatial links
	//the zeroth link is always on the same level.
	phi[i] += phi0[i] - phi0[NodeList[i].nn[0]];
	
	//The q-1 (and q-2) link(s) go back one level,
	//the q-3 (or q-2) link is on the same level. These
	//links are always computed in the same way.
	for(int mu = q-1; mu > (NodeList[i].fwdLinks); mu--) {
	  //cout<<"i="<<i<<" mu="<<mu<<" phi0="<<phi0[i]<<" phi0["<<NodeList[i].nn[mu]<<"]="<<phi0[NodeList[i].nn[mu]]<<endl;
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
	  for(int mu = q; mu < q+T_offset; mu++) {
	    //cout<<"i="<<i<<" mu="<<mu<<" phi0="<<phi0[i]<<" phi0["<<NodeList[i].nn[mu]<<"]="<<phi0[NodeList[i].nn[mu]]<<endl;
	    phi[i] += (phi0[i] - phi0[NodeList[i].nn[mu]])/(NodeList[i].temporal_weight);
	  }
	}    
      }
    }
  */
  return 0;  
}

Float Minv_phi(Float *phi, Float *b,
	       const vector<Vertex> NodeList, Param P) {
  // CG solutions to Mphi = b 
  //  see http://en.wikipedia.org/wiki/Conjugate_gradient_method
  int Levels = P.Levels;
  int diskN = endNode(Levels,P) + 1;
  int TotNumber = P.t*diskN;
  
  Float *res, *resNew, *pvec, *Mpvec, *pvec_tmp;  
  res      = new Float[TotNumber];
  resNew   = new Float[TotNumber];
  pvec     = new Float[TotNumber];
  Mpvec    = new Float[TotNumber];
  pvec_tmp = new Float[TotNumber];

  for(int i=0; i<TotNumber; i++) {
    res[i]      = 0.0;
    resNew[i]   = 0.0;
    pvec[i]     = 0.0;  
    Mpvec[i]    = 0.0;
    pvec_tmp[i] = 0.0;
    //cout<<"phi "<<phi[i]<<" b "<<b[i]<<endl;
  }
  
  Float alpha = 0.0, beta = 0.0, denom = 0.0;
  Float rsq = 0.0, rsqNew = 0.0, bsqrt = 0.0, truersq = 0.0;
  int  i;
  
  for(i = 0; i<TotNumber; i++){
    res[i] = b[i];
    pvec[i] = res[i];
    bsqrt += b[i]*b[i];
  }
  bsqrt = sqrt(bsqrt);
  
  int maxIter = P.MaxIter;
  Float resEpsilon = P.tol;
  // iterate till convergence
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
    cout << endl << "In CG at iteration = "<< k << " Montonicity "<<(beta < 1.0 ? "True" : "False")<<" Residue Squared  = " << rsq <<endl;

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

#endif
