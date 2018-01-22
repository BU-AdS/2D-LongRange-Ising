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
  bool lattice   = false; //if true, use a square lattice. If false, use the AdS space too.
  int MaxIter = 100000;
  double tol = pow(10,-6);
  int t = 32;
  double msqr = 4.0;
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

  int S1 = 32;
  int Lt = 32;
  int AdSVol = 32;
  int R  = 9;
  int surfaceVol = 0;
  int latVol = 0;
  double lambda = 1.0;
  double musqr  = -1.2;

  int n_therm=50000;
  int n_meas=100;
  int n_skip=1000;
  int n_wolff=20;
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

  //Positon on the poincare disk.
  std::complex<double> z;

  //The temporal weight (z dependent.)
  double temporal_weight = 1.0;

  //Field value at ths vertex.
  double phi = 0.0;
  
};


// Functions to calculate node/graph properties.
//----------------------------------------------

//Using the formula c(n) = (q-4)*c(n-1) - c(n-2) where c is the
//number of nodes on circumference at level n, we can construct
//the address of the end node on a given level for triangulation q:
//EN(lev,q) = SUM c(n) n=0..level
long unsigned int endNode(int lev, Param P);

//Get the z coordinates of every node on the Poincare disk 
void GetComplexPositions(std::vector<Vertex> &NodeList, Param& P);

//For each node n, with a link to another node,
//it checks that the neighbour table on the linked
//node contains the original node n as a neighbour.
void connectivityCheck(std::vector<Vertex> &NodeList, Param P);

//Truncate the graph according to the hyperbolic radius
//condition |z| < s.
void hypRadGraph(std::vector<Vertex> &NodeList, Param &P);

//Print the hyperbolic and poincare radii
void radiusCheck(std::vector<Vertex> &NodeList, Param P);

//Checks the area of each hyperbolc trangle in the graph
void CheckArea(const std::vector<Vertex> NodeList, Param P);

//Checks that every edge length is the same(ish)
void CheckEdgeLength(const std::vector<Vertex> NodeList, Param P);

// Functions to print, dump, or visualise data.
//----------------------------------------------
void PrintNodeTables(const std::vector<Vertex> NodeList, Param P);
void PrintComplexPositions(const std::vector<Vertex> NodeList, Param P);

//One-Size-Fits-All data file for lattice/analytical propagator data.
void DataDump(std::vector<Vertex> NodeList, double *phi, Param p, int level,
	      int t_range, int shift);

void visualiserSqr(std::vector<double> phi_cyl, double barr, Param p);
void visualiserAdS(std::vector<Vertex> NodeList, double barr, Param p);
void visualiserPhi2(double **phi_cyl, Param p, int iter);

// Functions to process data in the fly
//---------------------------------------
void correlatorsAdS(double **corr, double **corr_ave, int corr_norm,
		    std::vector<Vertex> NodeList, double avePhi, Param p);
void correlatorsSqr(double **corr, double **corr_ave, int corr_norm,
		    std::vector<double> phi, double avePhi, Param p);

void autocorrelation(double *PhiAb_arr, double avePhiAbs, int meas);




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
double s(std::complex<double> z);
double r(double s );
double d12(std::complex<double> z1, std::complex<double> z2);
double sigma(std::complex<double> z1, std::complex<double> z2, int t);
double s3p(int q);
double area3q(int q);
double areaGeneral(Param P, double A, double B, double C);
double centralRad(double s);
std::complex<double> DisktoUHP(std::complex<double> z);
std::complex<double> UHPtoDisk(std::complex<double> u);
std::complex<double> inversion(std::complex<double> z0, double r);
std::complex<double> squareInversion(std::complex<double>z0, double r1, double r2 );
double greens2D(std::complex<double> z, std::complex<double> w);
double greensM2D(std::complex<double> z, std::complex<double> w, Param p);
std::complex<double> newVertex(std::complex<double> z,std::complex<double> z0,int k, int q);

void radiusCheck(std::vector<Vertex> &NodeList, Param P);
void PrintNodeTables(const std::vector<Vertex> NodeList, Param P);

//- Edge length from center z = 0
double edgeLength(int q); 
#endif
