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
  if(argc > 1) p.init(argc, argv);
  
  //-- Populate problem dependent data. --//
  //Unless a specific positive input for the circumference is given,
  //it will be the nodes on the outer AdS perimeter.
  if(p.S1 <= 0)  p.S1 = endNode(p.Levels,p) - endNode(p.Levels-1,p);
  //Nodes on the 2D surface of the AdS space.
  p.surfaceVol = p.S1*p.Lt;
  //Nodes on the Poincare disk
  p.AdSVol = endNode(p.Levels,p) + 1;
  //Nodes in the entire lattice Volume.
  p.latVol = p.AdSVol * p.Lt;

  //Print endnode data for reference.
  for(int i=1;i<10;i++) cout<<"Endnode("<<i<<") = "<<endNode(i,p)<<endl;
  cout<<"Total Surface nodes = "<<p.surfaceVol<<endl;
  cout<<"Total AdS nodes     = "<<p.latVol<<endl;
  
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
  
  if(p.lat_type == ADS) {
    
    //Construct neighbour table.
    buildGraph(NodeList, p);
    //Get the z-coords and temporal weighting
    getComplexPositions(NodeList, p);
    //Get lattice/analytic scaling law
    latticeScaling(NodeList, p);
    //Calculate the one loop corrections, store in NodeList,
    //populate LR AdS Couplings
    //oneLoopCorrection(LR_couplings, NodeList, p);
    
    runMonteCarloAdSL(NodeList, p);    
  } else {
    MonteCarlo2DIsing Sim(p);
    Sim.runSimulation(p);
  }
  
  return 0;  
}
