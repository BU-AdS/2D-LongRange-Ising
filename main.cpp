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
#include "graph.h"
#include "cg.h"
#include "cg_multishift.h"
#include "eigen.h"
#include "update.h"

int main(int argc, char **argv) {

  cout<<setprecision(10);
  
  Param p;
  //Process Command line arguments
  if(argc > 1) p.init(argc, argv);
    
  //Print paramters
  p.print();

  //Print graph endnode info
  for(int i=1; i<20; i++) cout<<"Endnode("<<i<<") = "<<endNode(i,p)<<endl;

  //Total number of nodes in the graph.
  int TotNumber = (endNode(p.Levels,p) + 1) * p.t;

  cout<<"Total AdS nodes = "<<TotNumber<<endl;

  //Object to hold index positions of vertices
  vector<Vertex> NodeList(TotNumber);

  //-1 in NodeList indicates that node n has no connections.
  //This is used during construction to indicate if the node is yet
  //to be populated. During truncation, nodes to be removed are
  //assigned a position of -1, and all connections to that node are
  //removed.
  for(int n = 0; n <TotNumber;n++)
    for(int mu = 0; mu < p.q+2; mu++) 
      NodeList[n].nn[mu] = -1;
  
  //Construct neighbour table.
  BuildGraph(NodeList, p);

  //Get the z-coords and temporal weighting
  GetComplexPositions(NodeList, p);
    
  if(p.lattice == true) runMonteCarloSqr(NodeList, p);
  else runMonteCarloAdS(NodeList, p);

  /*
  //If the specified source position is < 0, place the point source
  //on the outer circumference.
  if(p.src_pos < 0) p.src_pos = (endNode(p.Levels,p) + 1) - 1;
  //if(p.src_pos < 0) p.src_pos = 0;
  
  //Debug tools
  //CheckEdgeLength(NodeList, p);
  //CheckArea(NodeList, p);  
  //PrintNodeTables(NodeList, p);
  //PrintComplexPositions(NodeList, p);
  //radiusCheck(NodeList, p);
  
  if(p.src_pos < 0 || p.src_pos > TotNumber) {
    cout<<"ERROR: Source Position must be g.e. 0 and l.e. "<<TotNumber-1;
    cout<<endl;
    exit(0);
  }
  if(NodeList[p.src_pos].pos == -1) {
    cout<<"ERROR: Source Position must be within the hyperbolic ";
    cout<<"radius, i.e., a point on the truncated graph."<<endl;
    exit(0);
  }
  
  //---------------//
  // Multishift CG //
  //---------------//

  cout<<"SOURCE = "<<p.src_pos<<endl;  
  int n_shift = p.n_shift;
  Float** phi = new Float*[n_shift];
  for(int i=0; i<n_shift; i++) {
    phi[i] = new Float[TotNumber];
    for(int j=0; j<TotNumber; j++) phi[i][j] = 0.0;
  }
  Float* b = new Float[TotNumber];
  for(int i=0; i<TotNumber; i++) b[i] = 0.0;
  
  b[p.src_pos] = 1.0;
  Minv_phi(phi[0], b, NodeList, p);
  //Minv_phi_ms(phi, b, NodeList, p);
    
  for(int i=0; i<n_shift; i++) {
    DataDump(NodeList, phi[i], p, p.Levels, p.t/2, i);
  }

  //Mphi_ev(NodeList, p);
  
  for(int i=0; i<n_shift; i++) delete[] phi[i];
  delete[] phi;
  delete[] b;
  */


  
  return 0;
}
