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


typedef enum couplingType_s {
  SR,
  POW,
  RAD
} couplingType;

class Param{

 public:

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
  couplingType coupling_type = RAD;  
  double t_scale = 1.0;
  
  char fname[256];

  int S1 = 16;
  int Lt = 64;
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
  int n_jkblock = 10;
  double delta_phi = 1.5;

  void usage(char **argv);
  void print();
  int init(int argc, char **argv, int *idx);
  
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
  double *Phi3_arr;
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
  double phi3Ave  = 0.0;
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

#endif
