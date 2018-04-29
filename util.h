#ifndef UTIL_H
#define UTIL_H
#include <complex>
#include <cstring>
#include <vector>

#define I std::complex<double>(0.0,1.0)

typedef enum theoryType_s {
  PHI4,
  ISING
} theoryType;


typedef enum latType_s {
  TWO_D,
  ADS
} latType;

typedef enum couplingType_s {
  SR,
  POW,
  RAD
} couplingType;

class Param{

 public:

  int q = 8;
  
  bool bc        = true;  //If true, use Dirichlet. If false, use Neumann
  bool Vcentre   = true;  //If true, place vertex at centre. 
                          //If false, use circumcentre.
  bool verbosity = false; //If true, print all data. If false, print summary.
  bool useWolff  = true;  //If true, use Wolff. Else, use SW.
  bool usePowLaw = true;  //If true, use the 1/r^a LR coupling, else use the
                          //user defined LR.

  int nThreads = 0;
  int nBlocks = 0;
  
  bool useGPUMetro = false;
  bool useGPUCluster = false;
  bool doMetroCheck = false;

  theoryType theory_type = ISING;
  latType lat_type = TWO_D;
  couplingType coupling_type = RAD;
  
  int MaxIter = 100000;
  double tol = pow(10,-6);
  double msqr = 1.0;
  double C_msqr = 10.0;
  
  double N_latt = 0.01;
  int n_shift = 1;
  double delta_msqr = 0.01;
  int Levels = 2;
  int src_pos = -1;
  double hyp_rad = 5.0;
  int r_min_pos = 0;
  int r_max_pos = 0;
  double t_scale = 1.0;
  
  char fname[256];

  int S1 = 16;
  int Lt = 64;
  int AdSVol = 0;
  int R = 9;
  int surfaceVol = 0;
  int latVol = 0;
  double musqr  = -1.2725;
  double J = 1.0;
  double h = 0.0;
  double lambda = 1.0;
  double sigma = 10.0;

  int n_metro_cool = 0;
  int n_therm = 1000;
  int n_meas = 1000;
  int n_write = 100;
  int n_skip = 100;
  int n_cluster = 8;
  double delta_phi = 1.5;

  void usage(char **argv);
  void print();
  int init(int argc, char **argv, int *idx);
  
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

class observables {
  
 private:
  int n_meas = 1;
  
 public:
  //Arrays holding measurements for error analysis
  double *E_arr;
  double *E2_arr;
  double *PhiAb_arr;
  double *Phi_arr;
  double *Phi2_arr;
  double *Phi4_arr;
  double *Suscep;
  double *SpecHeat;
  double *Binder;
  
  //Running averages
  double tmpE     = 0.0;  
  double aveE     = 0.0;
  double aveKE    = 0.0;
  double avePE    = 0.0;
  double aveE2    = 0.0;
  double avePhiAb = 0.0;
  double avePhi   = 0.0;
  double avePhi2  = 0.0;
  double avePhi4  = 0.0;
  double MagPhi   = 0.0;

  //temps
  double KE = 0.0;
  double PE = 0.0;
  
  //constructor
  observables(int meas);
  
  //destructor
  //~observable();
  
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
inline double edgeLength(int q){
  return sqrt( 1 - 4*sin(M_PI/q)*sin(M_PI/q) );
}


#endif
