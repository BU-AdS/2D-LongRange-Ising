#ifndef UTIL_H
#define UTIL_H
#include <complex>
#include <cstring>
#include <vector>

#define I std::complex<double>(0.0,1.0)

class Param{

 public:

  int q = 8;
  
  bool bc        = true;  //if true, use Dirichlet. If false, use Neumann
  bool Vcentre   = true;  //if true, place vertex at centre. If false, use circumcentre.
  bool verbosity = false; //if true, print all data. If false, print summary.
  bool lattice   = false;  //if true, use a square lattice. If false, use the AdS space too.
  int MaxIter = 100000;
  double tol = pow(10,-6);
  int t = 32;
  double msqr = 0.059;
  double C_msqr = 1.0;
  
  double N_latt = 1.0;
  int n_shift = 1;
  double delta_msqr = 0.01;
  int Levels = 2;
  int src_pos = -1;
  double hyp_rad = 5.0;
  int r_min_pos = 0;
  int r_max_pos = 0;
  
  char fname[256];

  int S1 = 0;
  int Lt = 0;
  int AdSVol = 0;
  int R = 9;
  int surfaceVol = 0;
  int latVol = 0;
  double lambda = 1.0;
  double musqr  = -1.275;

  int n_therm=100000;
  int n_meas=1000;
  int n_skip=1000;
  int n_wolff=8;
  double delta_phi = 1.5;

  void print();
  void init(int argc, char **argv);
};

class Vertex{
 public:
  //If the pos value is -1, it is not connected
  //to the graph.
  int pos = -1;
  
  //Neighbours for up to q=9 and 2 temporal directions.
  int nn[11] = {0,0,0,0,0,0,0,0,0,0,0};

  //How many forward links (important in the
  //buildGraph() function.
  int fwdLinks;

  //Positon on the Poincare disk.
  std::complex<double> z;

  //The temporal weight (z dependent.)
  double temporal_weight = 1.0;

  //M(x,x)^{-1} value at ths vertex.
  double oneLoopCorr = 0.0;
  
  //Field value at ths vertex.
  double phi = 0.0;
  
};


//Using the formula c(n) = (q-4)*c(n-1) - c(n-2) where c is the
//number of nodes on circumference at level n, we can construct
//the address of the end node on a given level for triangulation q:
//EN(lev,q) = SUM c(n) n=0..level
long unsigned int endNode(int lev, Param &P);

//Print the hyperbolic and poincare radii
void radiusCheck(std::vector<Vertex> &NodeList, Param P);

//Checks the area of each hyperbolc trangle in the graph
void checkArea(const std::vector<Vertex> NodeList, Param P);

//Checks that every edge length is the same(ish)
void checkEdgeLength(const std::vector<Vertex> NodeList, Param P);

//- Edge length from center z = 0
double edgeLength(int q); 


#endif
