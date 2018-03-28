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

int seed = clock();
//int seed = 1234;
mt19937 rng(seed);
uniform_real_distribution<double> unif(0.0,1.0);
int CACHE_LINE_SIZE = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);
int sze = 512;

#include "util.h"
#include "data_io.h"
#include "graph.h"
#include "cg.h"
#include "eigen.h"
#include "mc_update_sr.h"
#include "mc_update_lr.h"
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
      
  if(p.lat_type == ADS) {
   
    
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
  } else {

    //Nodes on the 2D surface of the AdS space.
    p.surfaceVol = p.S1*p.Lt;

    //Print paramters
    p.print();
    
    MonteCarlo2DIsing Sim(p);
    Sim.runSimulation(p);
  }
  
  return 0;  
}
