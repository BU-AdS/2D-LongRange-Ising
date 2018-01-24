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

// Solves lhs = A^(-1) rhs using multishift CG as defined in
// http://arxiv.org/pdf/hep-lat/9612014.pdf
// Assumes there are n_shift values in "shifts".
// If they are sorted s.t. the smallest shift is the smallest
// (worst-conditioned solve), set worst_first = true. 
// resid_freq_check is how often to check the residual of other solutions.
// This lets us stop iterating on converged systems. 
void cg_multishift(double **phi, double *phi0, int n_shift, int size,
		   int resid_freq_check, int max_iter, double eps,
		   double* shifts, vector<Vertex> NodeList,
		   Param param) {
  
  // Initialize vectors.
  double *r, *p, *Ap;
  double **p_s;
  double alpha, beta, beta_prev, rsq, rsqNew, bsqrt, tmp; 
  double *alpha_s, *beta_s, *zeta_s, *zeta_s_prev;
  int k,i,n;
  int n_shift_rem = n_shift; // number of systems to still iterate on. 

  // Allocate memory.
  alpha_s = new double[n_shift];
  beta_s = new double[n_shift];
  zeta_s = new double[n_shift];
  zeta_s_prev = new double[n_shift];
  
  p_s = new double*[n_shift];
  for (n = 0; n < n_shift; n++)
  {
    p_s[n] = new double[size];
  }

  r = new double[size];
  p = new double[size];
  Ap = new double[size];

  // Initialize values.
  rsq = 0.0; rsqNew = 0.0; bsqrt = 0.0; k=0;
  for (n = 0; n < n_shift; n++)
  {
    // beta_0, zeta_0, zeta_-1
    beta_s[n] = zeta_s[n] = zeta_s_prev[n] = 1.0;
    // alpha_0. 
    alpha_s[n] = 0.0;
  }
  beta = 1.0; alpha = 0.0;

  // Zero vectors;
  zero_vector(r, size); 
  zero_vector(p, size); zero_vector(Ap, size);
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq(phi0, size));
  
  // There can't be an initial guess... though it is sort of possible, in reference to:
  // http://arxiv.org/pdf/0810.1081v1.pdf
  
  // 1. x_sigma = 0, r = p_sigma = b.
  for (n = 0; n < n_shift; n++)
  {
    copy_vector(p_s[n], phi0, size);
    zero_vector(phi[n], size);
  }
  copy_vector(p, phi0, size);
  copy_vector(r, phi0, size);
  
  // Compute Ap.
  zero_vector(Ap, size);
  Mphi(Ap, p, NodeList, param);
  
  // Compute rsq.
  rsq = norm2sq(r, size);

  // iterate till convergence
  for(k = 0; k< max_iter; k++) {
    
    // 2. beta_i = - rsq / pAp. Which is a weird switch from the normal notation, but whatever.
    beta_prev = beta; 
    beta = -rsq/dot(p, Ap, size);
    //cout << "beta = " << beta << "\n";
    
    for (n = 0; n < n_shift_rem; n++)
    {
      // 3. Calculate beta_i^sigma, zeta_i+1^sigma according to 2.42 to 2.44.
      // zeta_{i+1}^sigma = complicated...
      tmp = zeta_s[n]; // Save zeta_i to pop into the prev zeta.
      zeta_s[n] = (zeta_s[n]*zeta_s_prev[n]*beta_prev)/(beta*alpha*(zeta_s_prev[n]-zeta_s[n]) + zeta_s_prev[n]*beta_prev*(1.0-shifts[n]*beta));
      zeta_s_prev[n] = tmp; 
      
      //cout << "zeta_n = " << zeta_s[n] << ", zeta_{n-1} = " << zeta_s_prev[n];
      
      // beta_i^sigma = beta_i zeta_{n+1}^sigma / zeta_n^sigma
      beta_s[n] = beta*zeta_s[n]/zeta_s_prev[n];
      
      // 4. x_s = x_s - beta_s p_s
      caxpy(-beta_s[n], p_s[n], phi[n], size); 
      
      //cout << ", beta_n = " << beta_s[n] << "\n"; 
    }

    // 5. r = r + beta Ap
    caxpy(beta, Ap, r, size);
    
    // Exit if new residual is small enough
    rsqNew = norm2sq(r, size);

    // Comment this out to save your sanity.
    std::cout << "[CG-M-STATUS]: Iter " << k+1 << " RelTol " << sqrt(rsqNew)/bsqrt << "\n";

    // The residual of the shifted systems is zeta_s[n]*sqrt(rsqNew). Stop iterating on converged systems.
    if (k % resid_freq_check == 0)
    {
      for (n = 0; n < n_shift_rem; n++)
      {
        if (zeta_s[n]*sqrt(rsqNew) < eps*bsqrt) // if the residual of vector 'n' is sufficiently small...
        {
          n_shift_rem = n;
          break;
        }
      }
    }

    if ((abs(zeta_s[0])*sqrt(rsqNew) < eps*bsqrt) || n_shift_rem == 0 || k == max_iter-1)
    {
      break;
    }
    
    
  
    // 6. alpha = rsqNew / rsq.
    alpha = rsqNew / rsq;
    rsq = rsqNew; 
    
    //cout << "alpha = " << alpha << "\n";  
    
    for (n = 0; n < n_shift_rem; n++)
    {
      // 7. alpha_s = alpha * zeta_s * beta_s / (zeta_s_prev * beta)
      alpha_s[n] = alpha*zeta_s[n]*beta_s[n]/(zeta_s_prev[n] * beta);
      //cout << "alpha_n = " << alpha_s[n] << "\n";
      
      // 8. p_s = zeta_s_prev r + alpha_s p_s
      caxpby(zeta_s[n], r, alpha_s[n], p_s[n], size);
    }
    
    // Compute the new Ap.
    cxpay(r, alpha, p, size);
    Mphi(Ap, p, NodeList, param);
  } 
    
  if(k == max_iter-1) {
    std::cout << "[CG-M-STATUS]: WARNING! Did not converge within maximum interations.\n";
  }

  k++;
  
  // Calculate explicit rsqs.
  for (n = 0; n < n_shift; n++) {
    zero_vector(Ap, size);
    Mphi(Ap, phi[n], NodeList, param);
    caxpy(shifts[n], phi[n], Ap, size);
    std::cout << "[CG-M-STATUS]: Shift " << n << " ShiftVal " << shifts[n] << " RelTol " <<  sqrt(diffnorm2sq(Ap, phi0, size))/bsqrt << "\n";
  }
  
  
  // Free all the things!
  delete[] r;
  delete[] p;
  delete[] Ap;
  
  for (i = 0; i < n_shift; i++)
  {
    delete[] p_s[i];
  }
  delete[] p_s;
  
  delete[] alpha_s;
  delete[] beta_s;
  delete[] zeta_s;
  delete[] zeta_s_prev;

  std::cout << "[CG-M-STATUS]: Complete! Iter " << k << "\n";
}

//Wrapper for AdS code.
void Minv_phi_ms(double **phi, double *phi0, vector<Vertex> NodeList, Param p){

  int n_shift = p.n_shift;
  int size = (endNode(p.Levels,p) + 1) * p.t;
  int resid_freq_check = 10;
  int max_iter = p.MaxIter;
  double msqr = p.msqr;
  double delta_msqr = p.delta_msqr;
  double eps = p.tol;
  double *shifts = (double*)malloc(n_shift*sizeof(double));
  for(int i=0; i<n_shift; i++) shifts[i] = p.C_msqr*msqr + i*delta_msqr;

  cg_multishift(phi, phi0, n_shift, size, resid_freq_check,
		max_iter, eps, shifts, NodeList, p);
  
  free(shifts);  
}

void latticeScaling(vector<Vertex> &NodeList, Param& p){
  /*
  int n_shift = p.n_shift;
  int latVol  = p.latVol;
  
  double* phi0= new double[latVol];
  double* phi = new double[latVol];
  double* b   = new double[latVol];  
  for(int s=0; s<sources; s++) { 
    for(int i=0; i<latVol; i++) {
      b[i] = 0.0;
      phi[i] = 0.0;
    }
    b[endNode(p.Levels-1,p) + 1 + s] = 1.0;    
    Minv_phi(phi, b, NodeList, p);
    cout<<phi[endNode(p.Levels-1,p) + 1 + s]<<endl;
    NodeList[endNode(p.Levels-1,p) + 1 + s].oneLoopCorr = phi[0][endNode(p.Levels-1,p) + 1 + s];
  }
  
  gsl_fit_wlinear(const double * x, const size_t xstride, const double * w, const size_t wstride, const double * y, const size_t ystride, size_t n, double * c0, double * c1, double * cov00, double * cov01, double * cov11, double * chisq);


  int n_shift = p.n_shift;
  int size = (endNode(p.Levels,p) + 1) * p.t;
  int resid_freq_check = 10;
  int max_iter = p.MaxIter;
  double msqr = p.msqr;
  double delta_msqr = p.delta_msqr;
  double eps = p.tol;
  double *shifts = (double*)malloc(n_shift*sizeof(double));
  for(int i=0; i<n_shift; i++) shifts[i] = p.C_msqr*msqr + i*delta_msqr;
  
  cg_multishift(phi, phi0, n_shift, size, resid_freq_check,
		max_iter, eps, shifts, NodeList, p);
  
  free(shifts);  
  */  
}



//This function will calculate the M(x,x)^{-1} elements on the qth
//portion of the graph. 
void oneLoopCorrection(vector<Vertex> &NodeList, Param& p){
  
  int latVol  = p.latVol;
  int sources = (endNode(p.Levels,p) - endNode(p.Levels-1,p)) / p.q ;  

  cout<<"Solving for "<<sources<<" sources."<<endl;
  
  double* phi = new double[latVol];
  double* b   = new double[latVol];  
  for(int s=0; s<sources; s++) { 
    for(int i=0; i<latVol; i++) {
      b[i] = 0.0;
      phi[i] = 0.0;
    }
    b[endNode(p.Levels-1,p) + 1 + s] = 1.0;    
    Minv_phi(phi, b, NodeList, p);
    cout<<phi[endNode(p.Levels-1,p) + 1 + s]<<endl;
    NodeList[endNode(p.Levels-1,p) + 1 + s].oneLoopCorr = phi[endNode(p.Levels-1,p) + 1 + s];
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
  
  delete[] phi;
  delete[] b;
}

// v1 = 0
void zero_vector(double* v1, const int size)
{
  for (int i = 0; i < size; i++)
    v1[i] = 0.0;
}

// v1 = v2
void copy_vector(double* v1, double* v2, const int size)
{
  for (int i = 0; i < size; i++)
    v1[i] = v2[i];
}

// v2 += alpha v1
void caxpy(const double alpha, double* v1, double* v2, const int size)
{
  for (int i = 0; i < size; i++)
    v2[i] += alpha*v1[i];
}

// v2 = v1 + alpha * v2
void cxpay(double* v1, const double alpha, double* v2, const int size)
{
  for (int i = 0; i < size; i++)
    v2[i] = v1[i] + alpha*v2[i];
}

// v2 = alpha v1 + beta v2
void caxpby(const double alpha, double* v1, const double beta, double* v2, const int size)
{
  for (int i = 0; i < size; i++)
    v2[i] = alpha*v1[i] + beta*v2[i];
}

// v1 dot v2
double dot(double* v1, double* v2, const int size)
{
  double dv = 0.0;
  for (int i = 0; i < size; i++)
    dv += v1[i]*v2[i];

  return dv;
}

// ||v1||^2
double norm2sq(double* v1, const int size)
{
  double dv = 0.0;
  for (int i = 0; i < size; i++)
    dv += v1[i]*v1[i];

  return dv;
}

// ||v1 - v2||^2
double diffnorm2sq(double* v1, double* v2, const int size)
{
  double dv = 0.0;
  for (int i = 0; i < size; i++)
    dv += (v1[i]-v2[i])*(v1[i]-v2[i]);

  return dv;
}
