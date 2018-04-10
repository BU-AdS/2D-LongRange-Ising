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
//int CACHE_LINE_SIZE = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);
int CACHE_LINE_SIZE = 64;
int sze = 512;

#include "util.h"
#include "graph.h"
#include "cg.h"
#include "mcPhiFourth2D.h"
#include "monte_carlo_ads_local.h"

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
#ifndef USE_GPU
  if(p.useGPUMetro || p.useGPUCluster) {
    printf("ERROR: Cannot use GPU routines if they haven't been built. ");
    printf("Please revise your options or rebuild to enable GPU support\n");
    exit(0);
  }
#endif
  
  if(p.lat_type == ADS) {

#if 0
    //Nodes on the 2D surface of the AdS space.
    p.surfaceVol = p.S1*p.Lt;
    //Nodes on the Poincare disk
    p.AdSVol = endNode(p.Levels,p) + 1;
    //Nodes in the entire lattice Volume.
    p.latVol = p.AdSVol * p.Lt;

    //Print paramters
    p.print();

    //Object to hold graph data
    vector<Vertex> NodeList(p.latVol);
    
    //-1 in NodeList indicates that node n has no connections.
    //This is used during construction to indicate if the node is yet
    //to be populated. During truncation, nodes to be removed are
    //assigned a position of -1, and all connections to that node are
    //removed.
    for(int n = 0; n < p.latVol; n++)
      for(int mu = 0; mu < p.q+2; mu++) 
	NodeList[n].nn[mu] = -1;
    
    //Construct neighbour table.
    buildGraph(NodeList, p);
    //Get the z-coords and temporal weighting
    getComplexPositions(NodeList, p);
    //Get lattice/analytic scaling law
    latticeScaling(NodeList, p);
    //Calculate the one loop corrections, store in NodeList,
    //populate LR AdS Couplings
    //oneLoopCorrection(LR_couplings, NodeList, p);

    //Print endnode data for reference.
    for(int i=1;i<10;i++) cout<<"Endnode("<<i<<") = "<<endNode(i,p)<<endl;
    cout<<"Total Surface nodes = "<<p.surfaceVol<<endl;
    cout<<"Total AdS nodes     = "<<p.latVol<<endl;
    
    runMonteCarloAdSL(NodeList, p);
#endif
  } else {    
    //Nodes on the 2D surface of the AdS space.
    p.surfaceVol = p.S1*p.Lt;
    
    //Print paramters
    p.print();
    
    if(p.theory_type == PHI4) {
      PhiFourth2D Sim(p);
      cout<<"Yo!"<<endl;
      Sim.runSimulation(p);
    }
    else {
      Ising2D Sim(p);
      Sim.runSimulation(p);
    }
  }
  
  return 0;  
}
