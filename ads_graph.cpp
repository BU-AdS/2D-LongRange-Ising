#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>

#define Float double

#include <util.h>
#include <graph.h>
#include <cg.h>
#include <eigen.h>

using namespace std;

int main(int argc, char **argv) {

  Param p;
  if(argc > 1) p.init(argc, argv);

  //If the specified source position is < 0, place the point source
  //on the outer circumference.
  if(p.src_pos < 0) p.src_pos = endNode(p.Levels-1,p) + (endNode(p.Levels,p) - endNode(p.Levels-1,p) )/2;  

  //Print paramters
  p.print();

  //Print graph endnode info
  for(int i=1; i<20; i++) cout<<"Endnode("<<i<<") = "<<endNode(i,p)<<endl;

  //Total number of nodes in the graph.
  int TotNumber = (endNode(p.Levels,p) + 1) * p.t;

  //Object to hold index positions of vertices
  vector<Vertex> NodeList(TotNumber);

  //Initialise. -1 in NodeList indicates that node
  //is not yet populated.
  for(int n = 0; n <TotNumber;n++)
    for(int mu = 0; mu < p.q+2; mu++) 
      NodeList[n].nn[mu] = -1;
  
  //Construct neighbour table and z-coords
  BuildGraph(NodeList, p);
  GetComplexPositions(NodeList, p);
  
  //Debug tools
#if 1
  
  ConnectivityCheck(NodeList, p);
  CheckEdgeLength(NodeList, p);
  CheckArea(NodeList, p);
  if(p.verbosity) {
    PrintNodeTables(NodeList, p);  
    PrintComplexPositions(NodeList, p);
  }
  
#endif

  //------------//
  // CG routine //
  //------------//
  
  vector<Float> b(TotNumber,0.0);
  vector<Float> phi(TotNumber,0.0);
  vector<Float> phi0(TotNumber,0.0);  
  
  if(p.src_pos < 0 || p.src_pos > TotNumber) {
    cout<<"ERROR: Source Position must be g.e. 0 and l.e. "<<TotNumber-1<<endl;
    exit(0);
  }
  
  b[p.src_pos] = 1.0;
  
  Float truesq = 0.0;
  truesq = Minv_phi(phi, phi0, b, NodeList, p);
  cout<<"Tolerance = "<<p.tol<<" True Residual = "<<sqrt(truesq)<<endl;  

  DataDump(NodeList, phi, p);  
  //Mphi_ev(NodeList, p);
  
  return 0;
}

