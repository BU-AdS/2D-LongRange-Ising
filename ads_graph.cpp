#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include "graph.h"
#include "eigen.h"
using namespace std;

void DataDump(vector<Vertex> NodeList, vector<double> phi, Param p);

int main(int argc, char **argv) {

  Param p;
  if(argc > 1) p.init(argc, argv);
  p.print();
   
  int TotNumber = (endNode(p.Levels,p) + 1) * p.t;  
  vector<Vertex> NodeList(TotNumber);
  vector<Vertex> AuxNodeList(TotNumber);

  //set default null. Better with constructor
  for(int n = 0; n <TotNumber;n++)
    for(int mu = 0; mu < p.q+2; mu++) {
      NodeList[n].nn[mu] = -1;
      AuxNodeList[n].nn[mu] = 0;
    }

  //Construct neighbour table and z-coords
  BuildGraph(NodeList, p);
  GetComplexPositions(NodeList, p);
  
  //Debug tools
#if 0
  BooleanCheck(NodeList, AuxNodeList, p);
  PrintNodeTables(AuxNodeList, p);
  PrintNodeTables(NodeList, p);
  PrintComplexPositions(NodeList, p);
  CheckEdgeLength(NodeList, p);
  CheckArea(NodeList, p);
#endif

  //------------//
  // CG routine //
  //------------//
  
  vector<double> b(TotNumber,0.0);
  vector<double> phi(TotNumber,0.0);
  vector<double> phi0(TotNumber,0.0);  
  
  if(p.src_pos < 0 || p.src_pos > TotNumber) {
    cout<<"ERROR: Source Position must be g.e. 0 and l.e. "<<TotNumber-1<<endl;
    exit(0);
  }
  
  b[p.src_pos] = 1.0;
  
  double truesq = 0.0;
  truesq = Minv_phi(phi, phi0, b, NodeList, p);
  cout<<"Tolerance = "<<p.tol<<" True Residual = "<<sqrt(truesq)<<endl;  
  DataDump(NodeList, phi, p);
  
  //Mphi_ev(NodeList, p);
  
  return 0;
}

void DataDump(vector<Vertex> NodeList, vector<double> phi, Param p) {

  long unsigned int TotNumber = (endNode(p.Levels,p) + 1) * p.t;
  long unsigned int j = p.src_pos;
  
  //Data file for lattice/analytical propagator data,
  //Complex positions (Poincare, UHP), etc.
  if(strncmp(p.fname, "", 1) == 0) 
    sprintf(p.fname, "q%d_Lev%d_T%d_msqr%.3f_src%d_%s_%s.dat",
	    p.q,
	    p.Levels, 
	    p.t, 
	    p.msqr,
	    p.src_pos,
	    p.bc == true ? "Dirichlet" : "Neumann",
	    p.Vcentre == true ? "Vertex" : "Circum");
  FILE *fp1;
  fp1=fopen(p.fname, "w");  
  
  double norm = 0.0;
  for(long unsigned int i = 0;i < TotNumber; i++) norm += phi[i]*phi[i];
  for(long unsigned int i = 0;i < TotNumber; i++) phi[i] /= sqrt(norm);  
  for(long unsigned int i = endNode(0,p)+1; i < endNode(p.Levels,p)+1; i++) {    
    if(i != j && (i > endNode(p.Levels-1,p) &&  i < endNode(p.Levels,p) ) ) {
      fprintf(fp1, "%e %e %e %e %e %e %e %e %e %e %e %e %s", 
	      atan( ( NodeList[i].z / NodeList[j].z).real() / (NodeList[i].z / NodeList[j].z).imag() ), //1 source/sink angle
	      phi[i], //2 lattice prop
	      greens2D(NodeList[j].z, NodeList[i].z), //3 real prop
	      (greens2D(NodeList[j].z, NodeList[i].z) - phi[i])/greens2D(NodeList[j].z, NodeList[i].z), //4 rel error
	      phi[i]/greens2D(NodeList[j].z, NodeList[i].z), //5 ratio
	      abs(NodeList[i].z) - abs(NodeList[j].z), //6 delta r (euclidean)
	      d12( (abs(NodeList[i].z) / abs(NodeList[j].z))*NodeList[j].z , NodeList[j].z), //7 delta r (hyperbolic)
	      d12(NodeList[i].z , NodeList[j].z), //8 delta s (hyperbolic)
	      DisktoUHP(NodeList[i].z).real(), //9 real UHP position
	      DisktoUHP(NodeList[i].z).imag(), //10 imag UHP position
	      abs(DisktoUHP(NodeList[i].z)), //11 abs UHP position
	      atan( DisktoUHP(NodeList[i].z).imag() / DisktoUHP(NodeList[i].z).real() ), //12 angle UHP position
	      //
	      abs((greens2D(NodeList[j].z, NodeList[i].z) - phi[i])/greens2D(NodeList[j].z, NodeList[i].z)) < 0.01 ? "<---\n" : "\n");
    }
  }
  
  fclose(fp1);
}
