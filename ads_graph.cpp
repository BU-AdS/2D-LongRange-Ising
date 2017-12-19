#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>

#define Float long double

#include <util.h>
#include <graph.h>
#include <cg.h>
#include <cg_multishift.h>
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
  if(p.verbosity) {
    ConnectivityCheck(NodeList, p);
    CheckEdgeLength(NodeList, p);
    CheckArea(NodeList, p);
    PrintNodeTables(NodeList, p);  
    PrintComplexPositions(NodeList, p);
  }
  
  if(p.src_pos < 0 || p.src_pos > TotNumber) {
    cout<<"ERROR: Source Position must be g.e. 0 and l.e. "<<TotNumber-1;
    cout<<endl;
    exit(0);
  }
  
  //-------------//
  // CG routines //
  //-------------//

  Float *phi = (Float*)malloc(TotNumber*sizeof(Float));
  Float *b   = (Float*)malloc(TotNumber*sizeof(Float));
  for(int i=0; i<TotNumber; i++) {
    phi[i] = 0.0;
    b[i]   = 0.0;
  }
  
  b[p.src_pos] = 1.0;
  
  Float truesq = 0.0;
  truesq = Minv_phi(phi, b, NodeList, p);
  cout<<"Tolerance = "<<p.tol<<" True Residual = "<<sqrt(truesq)<<endl;
  //DataDump(NodeList, phi, p);  

  int n_shift = p.n_shift;
  Float **phi_ms = (Float**)malloc(n_shift*sizeof(Float*));
  for(int i=0; i<n_shift; i++) {
    phi_ms[i] = (Float*)malloc(TotNumber*sizeof(Float));
    for(int j=0; j<TotNumber; j++) phi_ms[i][j] = 0.0;
  }
  
  Minv_phi_ms(phi_ms, b, NodeList, p);
  
  //Mphi_ev(NodeList, p);
  
  return 0;
}
