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
mt19937_64 rng(seed);
uniform_real_distribution<double> unif(0.0,1.0);

#include "util.h"
#include "data_io.h"
#include "graph.h"
#include "cg.h"
#include "eigen.h"
#include "update.h"

int main(int argc, char **argv) {

  cout<<setprecision(10);
  
  Param p;
  //Process Command line arguments
  if(argc > 1) p.init(argc, argv);
  
  //Populate problem dependent data
  p.surfaceVol = (endNode(p.Levels,p) - endNode(p.Levels-1,p))*p.t;
  p.S1 = endNode(p.Levels,p) - endNode(p.Levels-1,p);
  p.Lt = p.t;
  p.AdSVol = endNode(p.Levels,p) + 1;
  if(p.lattice == true) {    
    p.latVol = p.surfaceVol;
  }
  else p.latVol = p.AdSVol * p.Lt;

  for(int i=1;i<10;i++) cout<<"Endnode("<<i<<") = "<<endNode(i,p)<<endl;
  if(p.lattice == true) cout<<"Total Square nodes = "<<p.surfaceVol<<endl;
  else cout<<"Total AdS nodes = "<<p.latVol<<endl;
  
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
  
  if(!p.lattice) {
    //Construct neighbour table.
    buildGraph(NodeList, p);
    //Get the z-coords and temporal weighting
    getComplexPositions(NodeList, p);
    //Get lattice/analytic scaling law
    latticeScaling(NodeList, p);
    //Calculate the one loop corrections, store in NodeList.
    oneLoopCorrection(NodeList, p);
    //PrintNodeTables(NodeList, p);
    //connectivityCheck(NodeList, p);
  }
  
  runMonteCarlo(NodeList, p);

  return 0;
  
}
