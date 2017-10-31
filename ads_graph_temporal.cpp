#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <string>
#include "graph.h"
#include "eigen.h"
using namespace std;

int main(int argc, char *argv[]) {

  Param p;
  if(argc > 1) p.Levels = atoi(argv[1]);
  int q = p.q;
  p.bc = true;  // true is Dirichlet, false Neumann.
  
  //Checking...
  //for(int i=1; i<20; i++) cout<<"ENDNODE("<<i<<")="<<endNode(i,q)<<"   RATIO="<<setprecision(10)<<1.0*endNode(i,q)/endNode(i-1,q)<<endl;
  
  int TotNumber = (endNode(p.Levels,q) + 1) * p.t;
  vector<Vertex> NodeList(TotNumber);
  vector<Vertex> AuxNodeList(TotNumber);
    
  //set default null. Better with constructor
  for(int n = 0; n <TotNumber;n++)
    for(int mu = 0; mu < q+2; mu++) {
      NodeList[n].nn[mu] = -1;
      AuxNodeList[n].nn[mu] = 0;
    }

  //Construct neighbour table and z-coords
  BuildGraph(NodeList, p);
  GetComplexPositions(NodeList, p);

  //Debug tools
  BooleanCheck(NodeList, AuxNodeList, p);
  PrintNodeTables(AuxNodeList, p);
  PrintNodeTables(NodeList, p);
  //PrintComplexPositions(NodeList, p);
  //CheckEdgeLength(NodeList, p);
  
  
#if 1
  //CG routine
  vector<double> b(TotNumber,0.0);
  vector<double> phi(TotNumber,0.0);
  vector<double> phi0(TotNumber,0.0);  

  //int j = endNode(p.Levels-1,q);
  int j = 0;
  b[j] = 1.0; //point source at zero
  
  //cout<< endl << "=========  Before CG ========== " << endl;
  // for(int i = 0;i < TotNumber; i++) cout << phi[i]<< "  ";
  
  double truesq = 0.0;
  if(p.t == 1) truesq = Minv_phi(phi, phi0, b, NodeList, p);
  else truesq = Minv_phi_t(phi, phi0, b, NodeList, p);
  
  cout << endl <<  "True Error Squared at End of CG  = " << truesq  << endl;
  
  cout<< endl << "=========  After CG (Normalised) ========== " << endl;
  for(int i = 0;i < TotNumber; i++) truesq += phi[i]*phi[i];
  //  for(int i = 0;i < TotNumber; i++) cout << phi[i]/sqrt(truesq)<< "  ";

  /*  
  for(int i = 0;i < TotNumber; i++)
    if( i != j) {
      cout << phi[i]/sqrt(truesq)<< "  " << greens2D(NodeList[j].z, NodeList[i].z) << " " << (phi[i]/sqrt(truesq)) / greens2D(NodeList[j].z, NodeList[i].z)<< endl;
    }
  */
  cout << endl;
  
#endif
  
  
  //Mphi_ev(NodeList, p);
  Mphi_ev_t(NodeList, p);
  //eigenLaplace();
   
  return 0;
}

