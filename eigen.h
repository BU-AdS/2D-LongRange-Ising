#ifndef EIGEN_H
#define EIGEN_H
// Copyright (c) 2017 Evan S Weinberg
// A reference piece of code which computes matrix elements
// of a reference sparse Laplace operator to fill a dense
// Eigen matrix, then computes the spectrum and prints
// the Eigenvalues.

// This code is for real, symmetric matrices.
// This code lives on github at github.com/weinbe2/eigens-with-eigen/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>

// Borrow dense matrix eigenvalue routines.
#include <Eigen/Dense>

using namespace std; 
using namespace Eigen;

// This is just for convenience. By default, all matrices
// in Eigen are column major. You can replace "Dynamic"
// with a specific number to template just one size.
// You can also ust use "MatrixXd".
typedef Matrix<Float, Dynamic, Dynamic, ColMajor> dMatrix;

int Mphi_ev(vector<Vertex> NodeList, Param p) {
  Float *in_real;
  Float *out_real;
  int Levels = p.Levels;
  int N = p.t*(endNode(Levels,p)+1);
  Float *phi  = (Float*)malloc(N*sizeof(Float));
  Float *phi0 = (Float*)malloc(N*sizeof(Float));
  // Set output precision to be long.
  cout << setprecision(10);

  // Allocate.
  in_real = new Float[N];
  out_real = new Float[N];

  // Zero out.
  for (int i = 0; i < N; i++)
  {
    in_real[i] = out_real[i] = 0.0;
  }

  //////////////////////////
  // REAL, SYMMETRIC CASE //
  //////////////////////////

  std::cout << "Real, Symmetric case.\n";
  std::cout << "NeV Requested = " << N << "\n";
  p.print();
  
  // Allocate a sufficiently gigantic matrix.
  dMatrix mat_real = dMatrix::Zero(N, N);

  // Form matrix elements. This is where it's important that
  // dMatrix is column major.
  for (int i = 0; i < N; i++) {
    // Set a point on the rhs for a matrix element.
    // If appropriate, zero out the previous point.
    if (i > 0) in_real[i-1] = 0.0;
    in_real[i] = 1.0;

    // Zero out the "out" vector. I defined "laplace_1d" to
    // not require this, but I put this here for generality.
    for (int j = 0; j < N; j++) out_real[j] = 0.0;
    
    // PUT YOUR MAT-VEC HERE.
    // laplace_1d(out_real, in_real, N, m_sq);
    for(int i = 0; i< N; i++) phi0[i] = in_real[i];
    Mphi(phi, phi0, NodeList, p);
    for(int j = 0; j < N; j++) out_real[j] = phi[j];
    
    // Copy your output into the right memory location.
    // If your data layout supports it, you can also pass
    // "mptr" directly as your "output vector" when you call
    // your mat-vec.
    Float* mptr = &(mat_real(i*N));
    
    for (int j = 0; j < N; j++) mptr[j] = out_real[j];
    
  }

  // We've now formed the dense matrix. We print it here
  // as a sanity check if it's small enough.
  if (N <= 16)
  {
    std::cout << mat_real << "\n";
  }

  // Get the eigenvalues and eigenvectors.
  SelfAdjointEigenSolver<dMatrix> eigsolve_real(N);
  eigsolve_real.compute(mat_real, EigenvaluesOnly);
  
  // Remark: if you only want the eigenvalues, you can call
  // eigsolve_real.compute(mat_real, EigenvaluesOnly);

  // Print the eigenvalues.
  //std::cout << "The " << N <<" eigenvalues are:\n" << eigsolve_real.eigenvalues() << "\n\n";
  

  dMatrix evals = eigsolve_real.eigenvalues();
  
  FILE *fp;
  char efname[256];
  sprintf(efname, "EV%d_q%d_Lev%d_T%d_msqr%.3e_src%d_%s_%s.dat",
	  N, 
	  p.q,
	  p.Levels, 
	  p.t, 
	  (double)p.msqr,
	  p.src_pos,
	  p.bc == true ? "Dirichlet" : "Neumann",
	  p.Vcentre == true ? "Vertex" : "Circum");
  fp=fopen(efname, "w");
  for(int i=0; i<N; i++) {
    fprintf(fp, "%d %.16e\n", i, (double)evals(i) );
  }
  fclose(fp);
  
  // Print the eigenvectors if the matrix is small enough.
  // As a remark, this also shows you how to access the eigenvectors. 
  if (N <= 16)
  {
    for (int i = 0; i < N; i++)
    {
      // You can use "VectorXd" as the type instead.
      dMatrix evec = eigsolve_real.eigenvectors().col(i);
      
      // Print the eigenvector.
      std::cout << "Eigenvector " << i << " equals:\n" << evec << "\n\n";

      // You can also copy the eigenvector into another array as such:
      for (int j = 0; j < N; j++)
      {
        out_real[j] = evec(j);
      }
    }
  }

  // Clean up.
  delete[] in_real;
  delete[] out_real;

  return 0;
}

void corr_eigs(double **corr, Param p) {

  MatrixXd A(p.S1/2,p.S1/2);;
  for(int i=0; i<p.S1/2; i++)
    for(int j=0; j<p.S1/2; j++)
      A(i,j) = corr[i][j];

  EigenSolver<MatrixXd> solve(A);
  cout << "The eigenvalues are:" << endl << solve.eigenvalues() << endl;
}



#endif
