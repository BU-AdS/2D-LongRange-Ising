#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>

using namespace std;
mt19937_64 rng(137);
uniform_real_distribution<double> unif;
#define Float long double

#include <util.h>
#include <graph.h>
#include <cg.h>
#include <cg_multishift.h>
#include <eigen.h>
#include <update.h>

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

  //Object to hold index positions of vertices
  vector<Vertex> NodeList(TotNumber);

  //-1 in NodeList indicates that node n has no connections.
  //This is used during construction to indicate if the node is yet
  //to be populated. During truncation, nodes to be removed are
  //assined a position of -1, and all connections to that node are
  //removed.
  for(int n = 0; n <TotNumber;n++)
    for(int mu = 0; mu < p.q+2; mu++) 
      NodeList[n].nn[mu] = -1;
  
  //Construct neighbour table and z-coords
  BuildGraph(NodeList, p);
  GetComplexPositions(NodeList, p);
  //Truncate the graph to within a hyperbolic radius.
  hypRadGraph(NodeList, p);
  //PrintNodeTables(NodeList, p);
  radiusCheck(NodeList, p);
  
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
  Float* b   = new Float[TotNumber];
  for(int i=0; i<TotNumber; i++) b[i] = 0.0;
  
  b[p.src_pos] = 1.0;  

<<<<<<< HEAD
  Minv_phi_ms(phi, b, NodeList, p);
=======
  Minv_phi(phi[0], b, NodeList, p);
  //Minv_phi_ms(phi, b, NodeList, p);
>>>>>>> 19b88bac0aa9505717051b36cb9807727bcc07f2
    
  for(int i=0; i<n_shift; i++) {
    DataDump(NodeList, phi[i], p, p.Levels, p.t/2, i);
  }

  Mphi_ev(NodeList, p);
  
  for(int i=0; i<n_shift; i++) delete[] phi[i];
  delete[] phi;
  delete[] b;

  /*
    Jobs to do
    incorporate the the internal DoF.

    check that the cylynder is the same as DS code.
    use the tuned ADS mass to match the CFT

    find the factor that brings the kinetic term to the surface term.
    
    Kinetic term = Cosh(t) - cos(phi)

   */

  
  /*****************  Monte Carlo Update *********************/
  //p.S1 =  endNode( p.Levels,q) - endNode( p.Levels-1,q);
  //Test against    ./2d_phi4 -0.7 0.5 32 16384 16384
  //compare with Dave Schaich's code
  //p.S1 = endNode(p.Levels,p) - endNode(p.Levels-1,p);
  //p.Lt = p.S1;
  p.S1 = 32;
  p.Lt = 32;
  p.SurfaceVol = p.S1 * p.Lt;
  
  //Debug against Shaich code  mu^2 = -0.7 lambda = 0.5 Size = 32  "^-0.7,0.5" 32-50.csv
  //Table I  in  paper by David Schaich, Will Loinaz https://arxiv.org/abs/0902.0045
  
  vector<double> phi_cyl(p.SurfaceVol,0.0);
  vector<int> s(p.SurfaceVol,0.0);
  p.latVol = p.SurfaceVol;
  
  // p.cluster = (int *) malloc( p.latVol * sizeof (int));  //Swendsen Wang Data Structure
  // p.stack = (int *) malloc(p.latVol * sizeof (int));    //Wolff Data Structure
    
  double KE = 0.0, PE = 0.0;
  double Etot = 0.0;
  double mag_phi = 0.0;
  double delta_mag_phi = 0.0;
  
  for(int i = 0;i < p.SurfaceVol; i++) {
    phi_cyl[i] = 2.0*unif(rng) - 1.0;
    s[i] = (phi_cyl[i] > 0) ? 1:-1;
    mag_phi +=  phi_cyl[i];
  }
  
  cout<<"p.Levels = "<<p.Levels<<" Lattice Size "<<p.S1<<" x "<<p.Lt<<" = ";
  cout<<p.Lt*p.S1<<" Initial Mag is"<<mag_phi<<endl;

  Etot = action_phi(phi_cyl, s, p, KE, PE);
  
  cout<<" K = "<<KE<<" U = "<<PE<<" K + U = "<<KE+PE;
  cout<<" (Check) E = "<<Etot<<endl;

  
  Etot = action_phi(phi_cyl, s, p, KE, PE);

  cout<<"K = "<<KE<<" U = "<<PE<<" K + U = "<<KE+PE<<endl;
  cout<<"Etot = "<<Etot<<endl;
  cout<<"Etot - (K+U) = "<<Etot-(KE+PE)<<endl;
  
  for(int iter = 0;iter < p.n_therm; iter++) {
    metropolis_update_phi(phi_cyl, s, p, delta_mag_phi, iter);
    if((iter+1)%5000 == 0) cout<<"Therm sweep "<<iter+1<<endl;
  }
  
  for(int i = 0;i < p.SurfaceVol; i++) {
    mag_phi += phi_cyl[i];
  }

  double tmp     = 0.0;
  
  double AverE   = 0.0;
  double AverE2  = 0.0;
  double AvePhiAb= 0.0;
  double AvePhi  = 0.0;
  double AvePhi2 = 0.0;
  double AvePhi4 = 0.0;
  double MagPhi  = 0.0;
  double rhoVol  = 1.0/(double)p.SurfaceVol;
  double rhoVol2 = rhoVol*rhoVol;
  double rhoVol4 = rhoVol2*rhoVol2;
  
  for(int iter = 0;iter < p.n_skip*p.n_meas; iter++) {
    
    metropolis_update_phi(phi_cyl, s, p, delta_mag_phi, iter);
    
    mag_phi += rhoVol*delta_mag_phi;
    tmp      = action_phi(phi_cyl,s, p, KE, PE);
    AverE   += rhoVol*tmp;
    AverE2  += rhoVol2*tmp*tmp;
    
    MagPhi   = 0.0;
    for(int i = 0;i < p.SurfaceVol; i++) {
      MagPhi += phi_cyl[i];
    }
    MagPhi    *= rhoVol;
    
    AvePhiAb  += abs(MagPhi);
    AvePhi    += MagPhi;
    AvePhi2   += MagPhi*MagPhi;
    AvePhi4   += MagPhi*MagPhi*MagPhi*MagPhi;
    
    if((iter+1) % p.n_skip == 0) {      
      cout<<"Mearurement "<<(iter+1)/p.n_skip<<" Running averages."<<endl;
      cout<<setprecision(8);
      double norm = 1.0/(iter+1);
      
      cout<<"Average Mag. Density = "<<norm*mag_phi<<endl;
      cout<<"Average Eng. Density = "<<AverE*norm<<endl;
      cout<<"Ave |phi| = "<<AvePhiAb*norm<<endl;
      cout<<"Ave phi   = "<<AvePhi*norm<<endl;
      cout<<"Ave phi^2 = "<<AvePhi2*norm<<endl;
      cout<<"Ave phi^4 = "<<AvePhi4*norm<<endl;
      cout<<"Suscep    = "<<(AvePhi2*norm - pow(AvePhiAb*norm,2))/rhoVol<<endl;
      cout<<"Spec Heat = "<<(AverE2*norm - pow(AverE*norm,2))/rhoVol<<endl;
      cout<<"Binder    = "<<1.0 - AvePhi4/(3.0*AvePhi2*AvePhi2*norm)<<endl<<endl;
					
    }
  }
  
  return 0;
}
