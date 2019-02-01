#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>

using namespace std;

//int seed = clock();
int seed = 1234;
mt19937 rng(seed);
uniform_real_distribution<double> unif(0.0,1.0);
int CACHE_LINE_SIZE = 64;
int sze = 512;

#include "util.h"
#include "mcPhiFourth2D.h"

int main(int argc, char **argv) {

  cout<<setprecision(10);
  
  Param p;
  //Process Command line arguments
  for (int i=1; i<argc; i++){
    if(p.init(argc, argv, &i) == 0){
      continue;
    }
    printf("ERROR: Invalid option: %s\n", argv[i-1]);
    p.usage(argv);
    exit(0);
  }

  int k=0;
  if(argc > 1) p.init(argc, argv, &k);

  //Sanity checks
  //---------------------------
#ifndef USE_GPU
  if(p.useGPUMetro || p.useGPUCluster) {
    printf("ERROR: Cannot use GPU routines if they haven't been built. ");
    printf("Please revise your options or rebuild to enable GPU support\n");
    exit(0);
  }
#endif

  if(p.doMetroCheck && !p.useGPUMetro) {
    printf("WARNING: Cannot perform CPU metro check when not performing\n");
    printf("Metropolis on the GPU. Overriding doMetroCheck to false.\n");
    p.doMetroCheck = false;
  }

    
  //Points on the 2D surface
  p.surfaceVol = p.S1*p.Lt;
  
  //Print paramters
  p.print();
  
  if(p.theory_type == PHI4) {
    PhiFourth2D Sim(p);
    Sim.runSimulation(p);
  }
  else {
    Ising2D Sim(p);
    Sim.runSimulation(p);
  }
  
  
  return 0;  
}
