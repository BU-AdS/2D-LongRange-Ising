#ifndef HYP_UTIL_H
#define HYP_UTIL_H
#include <complex>
#include <cstring>
#include <vector>

#include "util.h"

#define I std::complex<double>(0.0,1.0)

/*
  Basic Hyperbolic Algebra. 

  Moebius Transforms 

  Hyperbolic reflection for geodesic in UHP

  Given Line: (x-a)^2 + y^2 = r^2 
  Determined by x = a \pm r at y = 0
  x = a, y = r


  RL(a,r, u) = a + r^2/(conj(u) - a)

  Preserves the line and swap a  to infty.

  Map to Disc: z = D(u) =  (u -I)/(1 -I * u) with u = x + i y
  Map to UHP   u = U(z) = (z + I)/(1 + I * z);

  Find UHP circle that hits  +/- theta_0 on  Disc

  |z - A|^2 = R^2  
  |z|^2 - |z| |A| cos(theta) + |A|^2 = R^2
  boundary 1 - |A| cos (theta) + |A|^2 = R^2
  pick A = real. when a = 0 with map

  Need 3 point to define the mobius. Circle to Circle. 
*/
  
std::complex<double> T(std::complex<double> z,  std::complex<double> w);
std::complex<double> R(std::complex<double> z, std::complex<double> omega);
std::complex<double> flip(std::complex<double> z, std::complex<double> z1, std::complex<double> z2);
std::complex<double> DisktoUHP(std::complex<double> z);
std::complex<double> UHPtoDisk(std::complex<double> u);
std::complex<double> inversion(std::complex<double> z0, double r);
std::complex<double> squareInversion(std::complex<double>z0, double r1, double r2 );
std::complex<double> newVertex(std::complex<double> z,std::complex<double> z0,int k, int q);

double s(std::complex<double> z);
double r(double s );
double d12(std::complex<double> z1, std::complex<double> z2);
double sigma(std::complex<double> z1, std::complex<double> z2, int t);
double s3p(int q);
double area3q(int q);
double areaGeneral(Param P, double A, double B, double C);
double centralRad(double s);
double greens2D(std::complex<double> z, std::complex<double> w);
double greensM2D(std::complex<double> z, std::complex<double> w, Param p);


#endif
