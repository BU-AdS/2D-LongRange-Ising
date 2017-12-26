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
  //PrintNodeTables(NodeList, p);
  
  //Check that hyperbolic radius and levels requested
  //are compatible. Shift radius to maximize connectivity.
  int r_min_pos = 0;
  Float r_min = 1.0;
  for(int n=endNode(p.Levels-1,p)+1; n<endNode(p.Levels,p)+1; n++){
    //if(s(NodeList[n].z) < p.hyp_rad) {
    //cout<<"ERROR: Please request either more levels or a smaller ";
    //cout<<"hypergeometric radius"<<endl;
    //exit(0);
    //}
    if(abs(NodeList[n].z) < r_min) {
      r_min = abs(NodeList[n].z);
      r_min_pos = n;
    }
  }
  p.hyp_rad = s(NodeList[r_min_pos].z);
  cout<<"HYP_RAD = "<<p.hyp_rad<<endl;
  
  hypRadGraph(NodeList, p);
  //PrintNodeTables(NodeList, p);
  //radiusCheck(NodeList, p);
  
  //If the specified source position is < 0, place the point source
  //on the outer circumference.
  int graph_size = 0;
  for(int n=0; n<endNode(p.Levels-1,p); n++)
    if(NodeList[n].pos != -1) {
      graph_size++;
      cout<<n<<endl;
    }
  
  if(p.src_pos < 0) p.src_pos = r_min_pos;

  
  
  //Debug tools
  //CheckEdgeLength(NodeList, p);
  //CheckArea(NodeList, p);  

  if(p.src_pos < 0 || p.src_pos > TotNumber) {
    cout<<"ERROR: Source Position must be g.e. 0 and l.e. "<<TotNumber-1;
    cout<<endl;
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
  for(int i=0; i<TotNumber; i++) b[i]   = 0.0;
  
  b[p.src_pos] = 1.0;  

  Minv_phi(phi[0], b, NodeList, p);
  //Minv_phi_ms(phi, b, NodeList, p);
    
  for(int i=0; i<n_shift; i++) {
    //DataDump(NodeList, phi[i], p, p.Levels, p.t/2, i);
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
  //p.S1 = endNode(p.Levels, p); // compare with Dave Schaich's code
  p.Lt = 32;
  p.S1 = 32;
  p.SurfaceVol = p.S1 * p.Lt;
  p.lambda = 0.5;
  p.msqr = 0.7;
  
  //Debug against Shaich code  mu^2 = -0.7 lambda = 0.5 Size = 32  "^-0.7,0.5" 32-50.csv
  //Table I  in  paper by David Schaich, Will Loinaz https://arxiv.org/abs/0902.0045
  
  vector<double> phi_cyl( p.SurfaceVol,0.0);
  vector<int> s( p.SurfaceVol,0.0);
  p.latVol =  p.SurfaceVol;
  
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
  
  cout <<"p.Levels =  "<< p.Levels  << "  Lattice Size is p.S1 x p.Lt =  "<< p.S1
       << " x "<<p.Lt<< "  Initial Mag is " << mag_phi << endl;
  Etot =   action_phi(phi_cyl,s, p,  KE, PE);
  
  cout << " KE = " << KE <<" PE " << PE<< " KE + PE "<< KE + PE << "   " << Etot << endl;
  
  int Thermalize = 100000;
  
  for(int iter = 0;iter < Thermalize;iter++) {
    metropolis_update_phi(phi_cyl, s, p, delta_mag_phi);
    // cout <<"   " <<  metropolis_update_phi(phi_cyl, s, p, delta_mag_phi) << "   " << delta_mag_phi<< endl;
  }
  
  for(int i = 0;i < p.SurfaceVol; i++) {
    mag_phi +=  phi_cyl[i];
  }
    
  int AverIter = 4000000;
  double AverE = 0.0;
  double tmpPhi = 0.0;
  double AveAbsPhi = 0.0;
  double AvePhi2 =0.0;
  double AvePhi4 = 0.0;
  double MagPhi = 0.0;
  double rhoVol = 1.0/(double) p.SurfaceVol;
  for(int iter = 0;iter < AverIter;iter++) {
    metropolis_update_phi(phi_cyl, s, p, delta_mag_phi);
    mag_phi += delta_mag_phi;
    AverE +=   action_phi(phi_cyl,s, p,  KE, PE);
    MagPhi = 0.0;
    for(int i = 0;i < p.SurfaceVol; i++) {
      tmpPhi = phi_cyl[i];
      MagPhi+= phi_cyl[i];
      //       AveAbsPhi +=  abs(tmpPhi);
      //tmpPhi *= tmpPhi;
      //AvePhi2 += tmpPhi;
      // tmpPhi *= tmpPhi;
      // AvePhi4 += tmpPhi;
    };
      AveAbsPhi +=  abs(MagPhi);
      MagPhi *= MagPhi;
      AvePhi2 += MagPhi;
      MagPhi *= MagPhi;
      AvePhi4 += MagPhi;
      
      // cout <<"   " <<  metropolis_update_phi(phi_cyl, s, p, delta_mag_phi) << "   " << delta_mag_phi<< endl;
  }
  
  cout<< " average Mag density =  " << rhoVol*mag_phi/(double)AverIter << " averge Energy density " << rhoVol*AverE/(double)AverIter << endl;
  
  // for(int i = 0;i < p.SurfaceVol; i++)
  //{
  //  AveAbsphi +=  abs(phi_cyl[i]);
  //};
  
  double quartPhi = AvePhi4/(double)AverIter ;
  double squaredPhi = AvePhi2/(double)AverIter;
  double cumulant = 1.0 - quartPhi / (3.0 * squaredPhi * squaredPhi);
  
  cout<< "  AveAbsphi  =  " <<  rhoVol*AveAbsPhi/(double)AverIter << endl;
  cout << "  Avephi2 = " <<squaredPhi<<"  Avephi4 = " << quartPhi<< "  Binder= " << cumulant << endl;
  
  return 0;
}
