#include <iostream>
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
uniform_real_distribution<double> unif;
#define Float long double

#include <util.h>
#include <graph.h>
#include <cg.h>
#include <cg_multishift.h>
#include <eigen.h>
#include <update_2dAdS.h>
//#include <update_2dSqr.h>


//-------------  Monte Carlo Update  ----------------//
void runMonteCarloAdS(vector<Vertex> NodeList, Param p) {
  
  p.AdSVol = endNode(p.Levels,p)+1;
  p.Lt = p.t;
  p.latVol = p.AdSVol * p.Lt;
  p.surfaceVol = endNode(p.Levels,p) - endNode(p.Levels-1,p);
  
  double KE = 0.0, PE = 0.0;
  double Etot = 0.0;
  double mag_phi = 0.0;
  double delta_mag_phi = 0.0;

  vector<int> s(p.latVol,0.0);
  
  for(int i = 0;i < p.latVol; i++) {
    NodeList[i].phi = 2.0*unif(rng) - 1.0;
    s[i] = (NodeList[i].phi > 0) ? 1:-1;
    mag_phi += NodeList[i].phi;
  }
  
  cout<<"Levels = "<<p.Levels<<" Lattice Size "<<p.AdSVol<<" x "<<p.Lt<<" = ";
  cout<<p.latVol<<" Initial Mag = "<<mag_phi<<endl;

  Etot = action_phi_AdS(NodeList, s, p, KE, PE);

  cout<<"K = "<<KE<<" U = "<<PE<<" K + U = "<<KE+PE<<endl;
  cout<<"Etot = "<<Etot<<endl;
  cout<<"Etot - (K+U) = "<<Etot-(KE+PE)<<endl;

  //Thermalisation
  for(int iter = 0;iter < p.n_therm; iter++) {
    //if p.n_wolff is set to zero, no Wolff steps occur.
    if((iter+1)%p.n_wolff == 0 && p.n_wolff != 0) {
      wolff_update_phi_AdS(NodeList, s, p, delta_mag_phi, iter);
    }

    metropolis_update_phi_AdS(NodeList, s, p, delta_mag_phi, iter);

    if((iter+1)%p.n_skip == 0) cout<<"Therm sweep "<<iter+1<<endl;
  }

  //Arrays holding measurements for error analysis
  double E_arr[p.n_meas];
  double E2_arr[p.n_meas];
  double PhiAb_arr[p.n_meas];
  double Phi_arr[p.n_meas];
  double Phi2_arr[p.n_meas];
  double Phi4_arr[p.n_meas];
  for(int i=0; i<p.n_meas; i++) {
    E_arr[i]     = 0.0;
    E2_arr[i]    = 0.0;
    PhiAb_arr[i] = 0.0;
    Phi_arr[i]   = 0.0;
    Phi2_arr[i]  = 0.0;
    Phi4_arr[i]  = 0.0;
  }

  //Running averages
  double tmpE     = 0.0;  
  double aveE     = 0.0;
  double aveE2    = 0.0;
  double avePhiAb = 0.0;
  double avePhi   = 0.0;
  double avePhi2  = 0.0;
  double avePhi4  = 0.0;
  double MagPhi   = 0.0;

  //Density of lattice sites
  double rhoVol  = 1.0/(double)p.latVol;
  double rhoVol2 = rhoVol*rhoVol;

  //correlation function arrays. One keeps the running average,
  //one is temporary to pass around single measurent values.
  //int s_max = p.surfaceVol/2;
  //if(p.surfaceVol/2 % 2 != 0) s_max =  
 
  double **corr_tmp = (double**)malloc((p.surfaceVol/2)*sizeof(double*));
  double **corr_ave = (double**)malloc((p.surfaceVol/2)*sizeof(double*));
  for(int i=0; i<p.surfaceVol/2; i++) {
    corr_tmp[i] = (double*)malloc((p.Lt/2)*sizeof(double));
    corr_ave[i] = (double*)malloc((p.Lt/2)*sizeof(double));
  }
  for(int i=0; i<p.surfaceVol/2; i++)
    for(int j=0; j<p.Lt/2; j++) 
      corr_ave[i][j] = 0.0;
  
  double **phi_sq_arr = (double**)malloc(p.surfaceVol*sizeof(double*));
  for(int i=0; i<p.surfaceVol; i++) {
    phi_sq_arr[i] = (double*)malloc(p.Lt*sizeof(double));
    for(int j=0; j<p.Lt; j++) phi_sq_arr[i][j] = 0.0;
  }

  int idx = 0;
  double norm;
  
  for(int iter = 0;iter < p.n_skip*p.n_meas; iter++) {
    
    if((iter+1)%p.n_wolff == 0 && p.n_wolff != 0) {
      wolff_update_phi_AdS(NodeList, s, p, delta_mag_phi, iter);
    }    
    metropolis_update_phi_AdS(NodeList, s, p, delta_mag_phi, iter);

    //Take measurements.
    if((iter+1) % p.n_skip == 0) {

      int offset = endNode(p.Levels-1,p) + 1;
      for(int i = 0;i < p.surfaceVol; i++) {
	for(int j=0; j<p.Lt; j++) 
	  phi_sq_arr[i][j] += pow(NodeList[i + offset + p.AdSVol*j].phi,2);
      }
      
      tmpE     = action_phi_AdS(NodeList, s, p, KE, PE);
      aveE    += rhoVol*tmpE;
      aveE2   += rhoVol2*tmpE*tmpE;

      MagPhi = 0.0;
      for(int i = 0;i < p.latVol; i++) MagPhi += NodeList[i].phi;
      MagPhi *= rhoVol;
      
      avePhiAb  += abs(MagPhi);
      avePhi    += MagPhi;
      avePhi2   += MagPhi*MagPhi;
      avePhi4   += MagPhi*MagPhi*MagPhi*MagPhi;
      
      E_arr[idx]     = rhoVol*tmpE;
      E2_arr[idx]    = rhoVol*tmpE*tmpE;
      PhiAb_arr[idx] = abs(MagPhi);
      Phi_arr[idx]   = MagPhi;
      Phi2_arr[idx]  = MagPhi*MagPhi;
      Phi4_arr[idx]  = MagPhi*MagPhi*MagPhi*MagPhi;
      
      idx++;
      
      cout<<setprecision(8);
      norm = 1.0/(idx);

      //Dump to stdout
      cout<<"Measurement "<<(iter+1)/p.n_skip<<" Sweep "<<iter+1<<endl;
      cout<<"Ave Energy= "<<aveE*norm<<endl;
      cout<<"Ave |phi| = "<<avePhiAb*norm<<endl;
      cout<<"Ave phi   = "<<avePhi*norm<<endl;
      cout<<"Ave phi^2 = "<<avePhi2*norm<<endl;
      cout<<"Ave phi^4 = "<<avePhi4*norm<<endl;
      cout<<"Suscep    = "<<(avePhi2*norm-pow(avePhiAb*norm,2))/rhoVol<<endl;
      cout<<"Spec Heat = "<<(aveE2*norm-pow(aveE*norm,2))/rhoVol<<endl;
      cout<<"Binder    = "<<1.0-avePhi4/(3.0*avePhi2*avePhi2*norm)<<endl;

      //Visualisation tools
      visualiser_AdS(NodeList, avePhiAb*norm, p);
      visualiser_phi2_AdS(phi_sq_arr, p, idx);      

      //Calculate correlaton functions and update the average.
      correlators(corr_tmp, corr_ave, idx, NodeList, avePhi*norm, p);
      sprintf(p.fname, "./data_dump/correlators.dat");
      FILE *fp1;
      fp1=fopen(p.fname, "w");
      for(int i=0; i<p.Lt/2; i++) fprintf(fp1, "%d %.8e\n", i, corr_tmp[0][i]);
      fclose(fp1);
      
      
      //if(idx%50 == 0) corr_eigs(corr_run, p);
      
    }
  }

  //autocorrelation(PhiAb_arr, avePhiAb, p.n_meas);
  
  //correlators(corr_run, corr_ave, idx, phi_cyl, avePhi*norm, p);
  //corr_eigs(corr_run, p);
  
  //free(corr_run);
  //free(corr_ave);
}


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

  cout<<"Total nodes = "<<TotNumber<<endl;

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
  
  runMonteCarloAdS(NodeList, p);
   
  return 0;
}
