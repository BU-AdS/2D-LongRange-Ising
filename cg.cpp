#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>

#include "cg.h" 
#include "util.h"
#include "hyp_util.h"
#include "fit/gsl_fit.h"

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
  
  //printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,(double)rsq,(double)truersq);
  
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
  
  // There can't be an initial guess... though it is sort of possible,
  // in reference to:
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
void Minv_phi_ms(double **phi, double *phi0,
		 vector<Vertex> NodeList, Param p){
  
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

  int latVol  = p.latVol;
  int lower = endNode(p.Levels-1,p)+1;

  int t1, t2, delta_t;
  double delta = 1.0 + sqrt(1 + p.msqr);
  double theta, r, r_p;
  complex<double> ratio;
  complex<double> snk;
  double* analytic_prop = new double[p.S1*p.Lt/2];
  double* xi_invariant  = new double[p.S1*p.Lt/2];

  // Construct the xi invariant and the analytic propagator
  
  int j = lower;
  complex<double> src = NodeList[j].z;  
  //Loop over timeslices
  for(int t=0; t<p.Lt/2; t++) {
    int T_offset = (endNode(p.Levels,p) + 1) * t;

    //Loop over outer circle of the Poincare disk
    for(long unsigned int k = 0; k < p.S1; k++) {

      //Construct i (full lattice index)
      int i = T_offset + lower + k;
      //Construct s (surface lattice index)
      int s = t*p.S1 + k;
      
      ratio = NodeList[i].z/NodeList[j].z;
      theta = atan2( ratio.imag() , ratio.real() );
      //index divided by disk size, using the int floor feature/bug,
      //gives the timeslice for each index.
      t1 = j / p.AdSVol;
      t2 = i / p.AdSVol;

      //Assume PBC.
      delta_t = (t2-t1);
      snk = NodeList[i].z;
      r   = abs(NodeList[i].z);
      r_p = abs(NodeList[j].z);
      
      xi_invariant[s]  = log( ((1-r)*(1-r_p))/(cosh(delta_t)*(1+r)*(1+r_p)
					       - 4*r*r_p*cos(theta)) );
      
      analytic_prop[s] = log( exp(-delta*sigma(src,snk,delta_t)) /
			      (1 - exp(-2*sigma(src,snk,delta_t))));
      
      //if(s<10) cout<<"s="<<s<<" xi="<<xi_invariant[s]<<" ap="<<analytic_prop[s]<<endl;
    }
  }

  double* c = new double[2];
  double* cov_ssq = new double[4];
  gsl_fit_linear(xi_invariant, 1, analytic_prop, 1, p.S1*p.Lt/2, &c[0], &c[1],
		 &cov_ssq[0], &cov_ssq[1], &cov_ssq[2], &cov_ssq[3]);

  cout<<"Target data"<<endl;
  cout<<"GSL data: C="<<c[0]<<" M="<<c[1]<<endl;
  cout<<"          covar00 = "<<cov_ssq[0]<<endl;
  cout<<"          covar01 = "<<cov_ssq[1]<<endl;
  cout<<"          covar11 = "<<cov_ssq[2]<<endl;
  cout<<"          sum_sq  = "<<cov_ssq[3]<<endl;

  double grad = c[1];
  double inter= c[0];
  
  double d_grad = 100;
  double d_inter = 100;

  double grad_tol = 1e-4;
  double inter_tol = 1e-3;

  //Search in wisdom file for a shorcut
  bool preTuned = false;
  bool tuneWin = false;
  ifstream fileInput;
  string line;
  char params[256];
  sprintf(params, "%d %d %d %.4f", p.q, p.Levels, p.Lt, p.msqr);
  char* search = params;
  unsigned int curLine = 0;
  // open file to search
  fileInput.open("ads_wisdom");
  if(fileInput.is_open()) {
    while(getline(fileInput, line)) { 
      curLine++;
      if (line.find(search, 0) != string::npos) {
	cout<<"Found data for your problem! Lucky you..."<<endl;
	preTuned = true;
	getline(fileInput, line);
	p.N_latt = stof(line);
	getline(fileInput, line);
	p.C_msqr = stof(line);
	cout<<"Using N_latt="<<p.N_latt<<endl;
	cout<<"Using C_msqr="<<p.C_msqr<<endl;
	fileInput.close();
      }
    }
    if(!preTuned)
      cout<<endl<<"AdS wisdom data not found. Strap in for some tuning..."<<endl<<endl;
  }
  else cout<<endl<<"AdS wisdom file not found. Strap in for some tuning..."<<endl<<endl;

  //Begin search for correct scaling factors.
  int iter = 0;
  //while(abs(d_grad) > grad_tol || abs(d_inter) > inter_tol) {
  while(abs(d_grad) > grad_tol) {
    
    double* phi_ave = new double[latVol];  
    double* phi = new double[latVol];
    double* b   = new double[latVol];

    //Take the average data from sources in the the qth sector
    //of the outer level. 
    int sources = (endNode(p.Levels,p) - endNode(p.Levels-1,p)) / p.q ;

    //initialise, invert, average.
    for(int i=0; i<latVol; i++) phi_ave[i] = 0.0;    
    for(int s=0; s<sources; s++) {
      for(int i=0; i<latVol; i++) {
	b[i] = 0.0;
	phi[i] = 0.0;
      }
      b[lower + s] = 1.0;    
      Minv_phi(phi, b, NodeList, p);
      for(int i=0; i<latVol; i++) phi_ave[i] += phi[i]/sources;
    }
    
    //Use current lattice normalisation.
    for(int i=0; i<latVol; i++) {
      phi_ave[i] = log(p.N_latt*phi_ave[i]);
    }
    
    //phi_ave now contains an averaged solution vector. We now
    //perform linear regression on this vector and the analytic
    //prop (w.r.t the xi invariant) to compute the relevant
    //scaling and normalsation.
    
    double* latt_prop  = new double[p.S1*p.Lt/2];
    //Loop over timeslices
    for(int t=0; t<p.Lt/2; t++) {
      int T_offset = (endNode(p.Levels,p) + 1) * t;
      //Loop over H2 disk
      for(long unsigned int k = 0; k < p.S1; k++) {
	
	//Construct i (full lattice AdS2p1 index)
	int i = T_offset + lower + k;
	//Construct s (surface lattice 2D index)
	int s = t*p.S1 + k;
	
	latt_prop[s] = phi_ave[i];
      }
    }
  
    //Extract linear fit data from log-log plot.
    gsl_fit_linear(xi_invariant, 1, latt_prop, 1, p.S1*p.Lt/2, &c[0], &c[1],
		   &cov_ssq[0], &cov_ssq[1], &cov_ssq[2], &cov_ssq[3]);
    
    cout<<"At iteration "<<iter<<endl;
    cout<<"GSL data: C="<<c[0]<<" M="<<c[1]<<endl;
    cout<<"          covar00 = "<<cov_ssq[0]<<endl;
    cout<<"          covar01 = "<<cov_ssq[1]<<endl;
    cout<<"          covar11 = "<<cov_ssq[2]<<endl;
    cout<<"          sum_sq  = "<<cov_ssq[3]<<endl;
    
    //Adjust the parameters and start over.
    d_grad = (c[1] - grad)/abs(grad);
    //d_inter= (c[0] - inter)/abs(inter);
    cout<<"D inter = "<<d_inter<<" D grad = "<<d_grad<<endl;
    cout<<"N_latt  = "<<p.N_latt<<" C_msqr = "<<p.C_msqr<<endl<<endl;
    
    double fac1 = abs(p.C_msqr);
    
    if( abs(d_grad) > grad_tol ) {
      if(c[1] > grad) p.msqr < 0 ? p.C_msqr += (fac1*abs(d_grad) + 1*grad_tol) : p.C_msqr -= (fac1*abs(d_grad) + 1*grad_tol);
      if(c[1] < grad) p.msqr < 0 ? p.C_msqr -= (fac1*abs(d_grad) + 1*grad_tol) : p.C_msqr += (fac1*abs(d_grad) + 1*grad_tol);
    }

    if(cov_ssq[0] != cov_ssq[0]) {
      cout<<"GSL failed!"<<endl;
      exit(0);
    }
    
    //double fac2 = 0.01;    
    //if( abs(d_grad) < grad_tol ) {
    //if(c[0] > inter) p.N_latt -= (fac2*abs(d_inter) + grad_tol/10);
    //if(c[0] < inter) p.N_latt += (fac2*abs(d_inter) + grad_tol/10);
    //}
    
    //Did it tune properly?
    //if(d_inter < inter_tol && d_grad < grad_tol) tuneWin = true;
    if(d_grad < grad_tol) tuneWin = true;
    
    delete[] phi_ave;
    delete[] phi;
    delete[] b;
    delete[] latt_prop;

    iter++;
    
  }

  
  
  //If this is a new problem, and it tuned properly, save the data.
  if(!preTuned && tuneWin) {
    FILE *file = fopen("ads_wisdom", "a");
    fprintf(file,"%d %d %d %.4f\n%f\n%f\n",
	    p.q, p.Levels, p.Lt, p.msqr,
	    p.N_latt,
	    p.C_msqr);
    fclose(file);
  }

  delete[] c;
  delete[] cov_ssq;

  delete[] analytic_prop;
  delete[] xi_invariant;
  
}

//This function will calculate the M(x,x)^{-1} elements on the qth
//sector of the Poincare disk.. 
void oneLoopCorrection(double *LR_coupling,
		       vector<Vertex> &NodeList,
		       Param& p){
  
  int sources = (endNode(p.Levels,p) - endNode(p.Levels-1,p)) / p.q ;
  cout<<"Solving for "<<sources<<" sources."<<endl;

  int pos = endNode(p.Levels-1,p) + 1;
  double* phi_ave = new double[p.latVol];
  for(int i=0; i<p.latVol; i++) phi_ave[i] = 0.0;
  double* phi = new double[p.latVol];
  double* b   = new double[p.latVol];  
  for(int s=0; s<sources; s++) {
    //initialise
    for(int i=0; i<p.latVol; i++) {
      b[i] = 0.0;
      phi[i] = 0.0;
    }
    //set source
    b[endNode(p.Levels-1,p) + 1 + s] = 1.0;
    //invert
    Minv_phi(phi, b, NodeList, p);
    //normalise 
    for(int i=0; i<p.latVol; i++) {
      phi[i] *= p.N_latt;
      phi_ave[i] += phi[i]/sources;
    }

    //----------------------//
    //Here is where the propagator has been normalised, so
    //we can use it as our long range coupling. We popuate
    //only the unique part of the i index, then copy the data
    //afterwards.
    for(int j=0; j<p.surfaceVol; j++){      
      //LR_coupling[s + j*p.surfaceVol] = phi[pos + j%p.S1 + (j/p.S1)*p.AdSVol];
      //cout<<s<<" "<<j<<" "<<j%p.S1<<" "<<(j/p.S1)*p.AdSVol<<" "<<LR_coupling[s + j*p.surfaceVol]<<endl;
    }
    
    //record one loop data
    for(int s=0; s<sources; s++) {
      NodeList[endNode(p.Levels-1,p) + 1 + s].oneLoopCorr = phi[pos + s];
    }
  }
  
  //Propagate the one loop correction data...
  for(int i=0; i<p.t; i++) {
    for(int q=0; q<p.q; q++) {    
      for(int s=0; s<sources; s++) {
	NodeList[pos + s + q*sources + i*p.AdSVol].oneLoopCorr = NodeList[pos + s].oneLoopCorr;
	//...and LR coupling data too.
	for(int j=0; j<p.surfaceVol; j++) {
	  //LR_coupling[s + q*sources + i*p.S1 + j*p.surfaceVol] = LR_coupling[s + ((p.surfaceVol - q*sources - i*p.S1 + j)%p.surfaceVol)*p.surfaceVol];
	  //cout<<s<<" "<<q<<" "<<i<<" "<<j<<" "<<LR_coupling[s + q*sources + i*p.S1 + j*p.surfaceVol]<<endl;
	}
      }
    }
  }

  double rad = 0.0;
  for(int i=0; i<p.S1; i++)
    rad += abs(NodeList[pos + i%p.S1 + (i/p.S1)*p.AdSVol].z)/p.S1;

  cout<<rad<<endl;
  //exit(0);
  vector<complex<double>> circum(p.S1);
  for(int i=0; i<p.S1; i++) {
    circum[i].real(rad*cos(2*i*(M_PI/p.S1)));
    circum[i].imag(rad*sin(2*i*(M_PI/p.S1)));
  }

  double delta = 1.0 + sqrt(1 + p.msqr);
  
  for(int i=0; i<p.surfaceVol; i++){
    for(int j=0; j<p.surfaceVol; j++){
      //index divided by disk size, using the int floor feature/bug,
      //gives the timeslice for each index.
      int t1 = i / p.S1;
      int t2 = j / p.S1;
      double delta_t = abs(t2-t1) > p.Lt/2 ? abs(p.Lt - abs(t2-t1)) : abs(t2-t1);
      delta_t = delta_t*(p.t_weight_scale);
      if(i != j) LR_coupling[i+j*p.surfaceVol] = exp(-delta*sigma(circum[i%p.S1], circum[j%p.S1], delta_t)) / (1 - exp(-2*sigma(circum[i%p.S1], circum[j%p.S1], delta_t)));
      //if(i != j) LR_coupling[i+j*p.surfaceVol] = sigma(circum[i%p.S1], circum[j%p.S1], delta_t);
      //else LR_coupling[i+j*p.surfaceVol] = 0.0;
      //if ( j == i+32 ) cout<<i<<" "<<j<<" "<<d12(circum[i%p.S1], circum[j%p.S1])<<" "<<LR_coupling[i + j*p.surfaceVol]<<endl;
    }
  }
  
  //exit(0);
  delete[] phi_ave;
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
