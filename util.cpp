#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <string>
#include <random>
#include <unistd.h>

#include "util.h"
#include "hyp_util.h"

using namespace std;

void Param::usage(char **argv) {

  printf("\n\n");
  printf("This is an exhaustive list of run-time options. Those not set at\n");
  printf("run time will be given a default value. Important input options will\n");
  printf("dumped to stdout at execution. Please ensure that they are sensible!\n\n");

  // The user is advised to study these parameters and set them accordingly.
  //-----------------------------------------------------------------------
  //printf("--latType <2D, AdS>              A 2D lattice is a periodic, Euclidean lattice. Use the --Lt <n> and --S1 <m> options to set the extents.\n");
  printf("--theory <phi4, ising>           Perform calculations for either the Ising model, or phi4.\n");  
  printf("--couplingType <SR, POW, RAD>    Set the coupling between the lattice sites. For SR the coupling is nearest neighbour.\n");
  printf("                                 For POW the coupling is of the form |x - y|^{-(d+sigma)}. For RAD the coupling is of the form (cosh(dt) - cos(dtheta))^{-(d+sigma)/2}\n");
  printf("--metroArch <GPU/CPU>            Use GPU or CPU for metropolis step\n");
  printf("--metroCheck <GPU/CPU>           Use GPU for metropolis step, check against the CPU result\n");
  printf("--clusterArch <GPU/CPU>          Use GPU or CPU for (GPU Wolff only, FIXME) cluster step\n");
  printf("--clusterAlg <WOLFF, SW>         Use either the Wolff or Swendsen-Wang cluster algorithm.\n");
  printf("--nBlocks <n>                    The number of GPU blocks.\n");
  printf("--nThreads <n>                   The number of GPU threads per block.\n");
  printf("--Lt <n>                         The temporal extent of the 2D lattice.\n");
  printf("--S1 <n>                         The spatial extent of the 2D lattice..\n");
  printf("--muSqr <Float>                  The phi^4 mass^2 term. This is usually negative...\n");
  printf("--lambda <Float>                 The phi^4 interaction term.\n");
  printf("--sigma <Float>                  The exponent in the LR power law.\n");
  printf("--J <Float>                      The (constant) coupling in the Ising model.\n");
  printf("--h <Float>                      The external field strength in the Ising model.\n");
  printf("--tScale <Float>                 The value of alpha in cosh(alpha * delta_t) in radial LR coupling.\n");
  //MC params
  printf("--nMetroCool <n>                 The number of pure Metropolis steps to perfom in the initial cooldown\n");
  printf("--nTherm <n>                     The number of thermalisation steps (1 x metro + nCluster x Cluster)\n");
  printf("--nSkip  <n>                     The number of samples to skip between measuremnets.\n");
  printf("--nMeas  <n>                     The number of measurements to make.\n");
  printf("--nWrite <n>                     Jackknife and dump the data every nth measurement.\n");
  printf("--nCluster <n>                   The number of cluster sweeps to make per Metropolis step.\n");
  printf("--nJkBlock <n>                   The number jackknife blocks.\n");
  printf("--deltaPhi <float>               The parameter that governs how much an individual phi value is changed in the Metropolis test.\n");
  
  printf("--verbosity <V/v, Q/q>           Sets the verbosity as either Verbose or Quiet.\n");
  
}

int Param::init(int argc, char **argv, int *idx) {

  int ret = -1;
  int i = *idx;

  if( strcmp(argv[i], "--help")== 0){
    usage(argv);
  }

  //Theory
  if( strcmp(argv[i], "--theory") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string theory_in(argv[i+1]);
    if (theory_in == "phi4" ||
	theory_in == "PHI4" ||
	theory_in == "Phi4") {
      theory_type = PHI4;
      
    } else if(theory_in == "ising" ||
	      theory_in == "Ising" ||
	      theory_in == "ISING") {
      theory_type = ISING;
    } else {
      cout<<"Invalid theory type ("<<theory_in<<") given. Options are "<<endl;
      cout<<"phi4: Phi fourth theory"<<endl;
      cout<<"ising: The Ising model."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }
  
  //Coupling Type
  if( strcmp(argv[i], "--couplingType") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    std::string couplingType_in(argv[i+1]);    
    if (couplingType_in == "sr"  ||
	couplingType_in == "SR"  ||
	couplingType_in == "Sr") {
      coupling_type = SR;      
    } else if(couplingType_in == "pow" ||
	      couplingType_in == "POW" ||
	      couplingType_in == "Pow") {
      coupling_type = POW;
      usePowLaw = true;
    } else if(couplingType_in == "rad" ||
	      couplingType_in == "RAD" ||
	      couplingType_in == "Rad") {
      coupling_type = RAD;
      usePowLaw = false;      
    } else {
      cout<<"Invalid coupling type "<<couplingType_in<<" given. Options are "<<endl;
      cout<<"SR: Short range (local) interacton"<<endl;
      cout<<"Pow: Long range (non-local) |x-y|^{d+s}"<<endl;
      cout<<"Rad: Long range (non-local) (cosh(dt) - cos(dtheta))^{d+s}"<<endl;    
      exit(0);
    }

    i++;
    ret = 0;
    goto out;
  } 

  //Metro Arch
  if( strcmp(argv[i], "--metroArch") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    std::string ARCH(argv[i+1]);
    if (ARCH == "GPU" || 
	ARCH == "gpu" ||
        ARCH == "device" ) {
      useGPUMetro = true;
    } else if (ARCH == "CPU" || 
	       ARCH == "cpu" ||
	       ARCH == "host" ) {
      useGPUMetro = false;
    } else {
      cout<<"Invalid Metropolis architecture ("<<ARCH<<") given. Use GPU or CPU (...abacus not implemented)."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //Metro Check
  if( strcmp(argv[i], "--metroCheck") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    std::string CHECK(argv[i+1]);
    if (CHECK == "true" || 
	CHECK == "TRUE" ||
        CHECK == "yes" ) {
      doMetroCheck = true;
    } else if (CHECK == "false" || 
	       CHECK == "FALSE" ||
	       CHECK == "no" ) {
      doMetroCheck = false;
    } else {
      cout<<"Invalid Metropolis check ("<<CHECK<<") given. Use TRUE or FALSE."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }
  
  //Cluster Arch
  if( strcmp(argv[i], "--clusterArch") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    std::string ARCH(argv[i+1]);
    if (ARCH == "GPU" || 
	ARCH == "gpu" ||
        ARCH == "device" ) {
      useGPUCluster = true;
    } else if (ARCH == "CPU" || 
	       ARCH == "cpu" ||
	       ARCH == "host" ) {
      useGPUCluster = false;
    } else {
      cout<<"Invalid cluster architecture ("<<ARCH<<") given. Use GPU or CPU (...abacus not implemented)."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }
  
  //Cluster Algorithm
  if( strcmp(argv[i], "--clusterAlg") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    std::string CLUSTER(argv[i+1]);
    if (CLUSTER == "wolff" || 
	CLUSTER == "Wolff" ||
	CLUSTER == "WOLFF" ) {
      useWolff = true;
    } else if (CLUSTER == "SW" || 
	       CLUSTER == "sw" ||
	       CLUSTER == "Sw" ) {
      useWolff = false;
    } else {
      cout<<"Invalid cluster algorithm ("<<CLUSTER<<") given. Use WOLFF for Wolff or SW for Swendsen Wang."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //Blocks
  if( strcmp(argv[i], "--nBlocks") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    nBlocks = atoi(argv[i+1]);
    if(nBlocks < 0) {
      cout<<"Invalid number of blocks ("<<nBlocks<<") given. Please ensure that nBlocks > 0, or nBlocks = 0 for automatic assignment."<<endl;
      exit(0);
    }    
    i++;
    ret = 0;
    goto out;
  }

  //Blocks
  if( strcmp(argv[i], "--nBlocks") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    nBlocks = atoi(argv[i+1]);
    if(nBlocks < 0) {
      cout<<"Invalid number of blocks ("<<nBlocks<<") given. Please ensure that nBlocks > 0, or nBlocks = 0 for automatic assignment."<<endl;
      exit(0);
    }    
    i++;
    ret = 0;
    goto out;
  }
  
  //Threads
  if( strcmp(argv[i], "--nThreads") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    nThreads = atoi(argv[i+1]);
    if(nThreads < 0) {
      cout<<"Invalid number of threads ("<<nThreads<<") given. Please ensure that nThreads > 0, or nThreads = 0 for automatic assignment."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //Temporal sites on the boundary
  if( strcmp(argv[i], "--Lt") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    Lt = atoi(argv[i+1]);
    if(Lt<2) {
      cout<<"Invalid number of temporal sites ("<<Lt<<") given. Please ensure that Lt >= 2."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  
  //Spatial sites on the boundary
  if( strcmp(argv[i], "--S1") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    S1 = atoi(argv[i+1]);
    if(S1<2) {
      cout<<"Invalid number of spatial sites ("<<S1<<") given. Please ensure that S1 >= 2."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //Boundary (mu) mass^2
  if( strcmp(argv[i], "--muSqr") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    musqr = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  //Boundary (lambda) interaction
  if( strcmp(argv[i], "--lambda") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    lambda = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  //Long-range coupling exponent
  if( strcmp(argv[i], "--sigma") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    sigma = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  //Ising (constant) coupling
  if( strcmp(argv[i], "--J") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    J = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  //Ising external field
  if( strcmp(argv[i], "--h") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    h = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  //Temporal scaling
  if( strcmp(argv[i], "--tScale") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    t_scale = atof(argv[i+1]);
    if(t_scale < 0) {
      cout<<"Invalid t_scale ("<<t_scale<<") given. Please ensure that tscale > 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }
  
  //--------------------//
  //     MC params      //
  //--------------------//
  
  //Metro cool-down steps.
  if( strcmp(argv[i], "--nMetroCool") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    n_metro_cool = atoi(argv[i+1]);
    if(n_metro_cool < 0) {
      cout<<"Invalid number of metro cool down steps ("<<n_metro_cool<<") given. Please ensure that n_metro_cool >= 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //Thermalisation sweeps
  if( strcmp(argv[i], "--nTherm") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    n_therm = atoi(argv[i+1]);
    if(n_therm < 0) {
      cout<<"Invalid number of thermalisation steps ("<<n_therm<<") given. Please ensure that n_therm >= 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //Measurements
  if( strcmp(argv[i], "--nMeas") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    n_meas = atoi(argv[i+1]);
    if(n_meas < 0) {
      cout<<"Invalid number of measurements ("<<n_meas<<") given. Please ensure that n_meas >= 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //Measurements
  if( strcmp(argv[i], "--nWrite") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    n_write = atoi(argv[i+1]);
    if(n_write <= 0) {
      cout<<"Invalid number of measurements ("<<n_write<<") given. Please ensure that n_write > 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

   //Skip steps
  if( strcmp(argv[i], "--nSkip") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    n_skip = atoi(argv[i+1]);
    if(n_skip < 0) {
      cout<<"Invalid number of skip steps ("<<n_skip<<") given. Please ensure that n_skip >= 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //Cluster steps
  if( strcmp(argv[i], "--nCluster") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    n_cluster = atoi(argv[i+1]);
    if(n_cluster < 0) {
      cout<<"Invalid number of cluster steps ("<<n_cluster<<") given. Please ensure that n_cluster >= 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

    //Jackknife blocks
  if( strcmp(argv[i], "--nJkBlock") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    n_jkblock = atoi(argv[i+1]);
    if(n_jkblock < 0) {
      cout<<"Invalid number of jackknife blocks ("<<n_jkblock<<") given. Please ensure that n_jkblock > 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }
  
  //Delta phi
  if( strcmp(argv[i], "--deltaPhi") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    delta_phi = atof(argv[i+1]);
    if(t_scale == 0) {
      cout<<"Invalid t_scale ("<<t_scale<<") given. Please ensure that tscale != 0.0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //Run-time Verbosity
  if( strcmp(argv[i], "--verbosity") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string verbose(argv[i+1]);
    if (verbose == "V" || verbose == "v") {
      verbosity = true;
    } else if(verbose == "Q" || verbose == "q") {
      verbosity = false;
    } else {
      cout<<"Invalid Verbosity condition ("<<verbose<<") given. Use v/q for verbose/quiet"<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  } 


  //Sanity checks...
  // if(nBlocks*nThreads != S1*Lt) {
  //   cout<<"nBlocks * nThreads is not the same at S1 * Lt."<<endl;
  //   cout<<"nBlocks = "<<nBlocks<<endl;
  //   cout<<"ntThreads = "<<nThreads<<endl;
  //   cout<<"S1 = "<<S1<<endl;
  //   cout<<"Lt = "<<Lt<<endl;
  // }

   
 out:
  *idx = i;
  return ret ;
  
}

void Param::print() {

  bool useGPU = false;
  bool useOMP = false;

#ifdef USE_GPU
  useGPU = true;
#endif
#ifdef USE_OMP
  useOMP = true;
#endif
  
  cout<<endl;
  cout<<"********************************"<<endl;
  cout<<"*      Parameter status        *"<<endl;
  cout<<"********************************"<<endl;
  cout<<endl;
  cout<<"Architectures compiled: "<<(useGPU == true ? "GPU and " : "")<<(useOMP == true ? "CPU (OMP)" : "CPU (serial)")<<endl;
  cout<<"Architecture for Metropolis: "<<(useGPUMetro == true ? "GPU" : (useOMP == true ? "CPU (OMP)" : "CPU (serial)"))<<endl;
  cout<<"Architecture for Cluster: "<<(useGPUCluster == true ? "GPU" : (useOMP == true ? "CPU (OMP)" : "CPU (serial)"))<<endl;
  cout<<"Theory type = "<<(theory_type == PHI4 ? "Phi Fourth" : "Ising")<<endl;
  cout<<"Coupling type = "<<(coupling_type == SR ? "Short-Range" :
			     coupling_type == POW ? "Power Law" : "Radial")<<endl;
  cout<<"Cluster Algorithm = "<<(useWolff == true ? "Wolff" : "Swendsen-Wang")<<endl;
  cout<<"TimeSlices = "<<Lt<<endl;
  cout<<"Circumference = "<<S1<<endl;
  if(theory_type == PHI4) {
    cout<<"(Phi4) Boundary Mass squared (mu^2) = "<<musqr<<endl;  
    cout<<"(Phi4) Boundary coupling (Lambda) = "<<lambda<<endl;
  } else {
    cout<<"(Ising) External field strength (h) = "<<h<<endl;
    cout<<"(Ising) Constant coupling (J) = "<<J<<endl;
  }   
  if(coupling_type != SR) cout<<"Long-Range coupling sigma = "<<sigma<<endl;
  cout<<"Temporal Scaling = "<<t_scale<<endl;
  
  cout<<endl<<"* MC Params *"<<endl;
  cout<<"Metro Cool-Down = "<<n_metro_cool<<endl;
  cout<<"Thermalisations = "<<n_therm<<endl;
  cout<<"Skip steps = "<<n_skip<<endl;
  cout<<"Measurements = "<<n_meas<<endl;
  cout<<"Cluster steps per Metro = "<<n_cluster<<endl;

  cout<<endl;
}

// Container class for observable quantities.
observables::observables(int meas) {

  n_meas = meas;
  
  E_arr = new double[n_meas];
  E2_arr = new double[n_meas];
  PhiAb_arr = new double[n_meas];
  Phi_arr = new double[n_meas];
  Phi2_arr = new double[n_meas];
  Phi3_arr = new double[n_meas];
  Phi4_arr = new double[n_meas];
  Suscep = new double[n_meas];
  SpecHeat = new double[n_meas];
  Binder = new double[n_meas];
  for(int i=0; i<n_meas; i++) {
    E_arr[i]     = 0.0;
    E2_arr[i]    = 0.0;
    PhiAb_arr[i] = 0.0;
    Phi_arr[i]   = 0.0;
    Phi2_arr[i]  = 0.0;
    Phi3_arr[i]  = 0.0;
    Phi4_arr[i]  = 0.0;
    Suscep[i]    = 0.0;
    SpecHeat[i]  = 0.0;
    Binder[i]    = 0.0;
  }

  tmpE     = 0.0;  
  aveE     = 0.0;
  aveKE    = 0.0;
  avePE    = 0.0;
  aveE2    = 0.0;
  avePhiAb = 0.0;
  avePhi   = 0.0;
  avePhi2  = 0.0;
  phi3Ave  = 0.0;
  avePhi4  = 0.0;
  MagPhi   = 0.0;
   
}


