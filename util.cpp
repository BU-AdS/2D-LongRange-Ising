#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>

#include "util.h"
#include "hyp_util.h"

using namespace std;

void Param::init(int argc, char **argv) {
  
  std::string BC(argv[1]);
  if (BC == "D" || BC == "d") {
    bc = true;
  } else if (BC == "N" || BC == "n") {
    bc = false;
    } else {
    cout<<"Invalid boundary condition given. Use D/d for Dirichlet or N/n for Neumann."<<endl;
    exit(0);
  }
  
  std::string Centre(argv[2]);
  if (Centre == "V" || Centre == "v") {
    Vcentre = true;
  } else if (Centre == "C" || Centre == "c") {
    Vcentre = false;
  } else {
    cout<<"Invalid centre condition given. Use V/v for Vertexcentred or C/c for Circumcentred."<<endl;
    exit(0);
  }
  
  std::string verbose(argv[3]);
  if (verbose == "V" || verbose == "v") {
    verbosity = true;
  } else if(verbose == "Q" || verbose == "q") {
    verbosity = false;
  } else {
    cout<<"Invalid Verbosity conditions given. Use v/q for verbose/quiet"<<endl;
    exit(0);
  }

  std::string latType_in(argv[4]);
  
  if (latType_in == "sq_local"  ||
      latType_in == "SQ_LOCAL"  ||
      latType_in == "Sq_Local") {
    latType = SQ_LOCAL;
    
  } else if(latType_in == "sq_nonlocal" ||
	    latType_in == "SQ_NONLOCAL" ||
	    latType_in == "Sq_NonLocal") {
    latType = SQ_NONLOCAL;
    
  } else if(latType_in == "ads" ||
	    latType_in == "AdS") {
    latType = ADS;
    
  } else {
    cout<<"Invalid Lattice type given. Options are "<<endl;
    cout<<"SQ_LOCAL: square lattice with local interacton"<<endl;
    cout<<"SQ_NONLOCAL: square lattice with non-local interacton"<<endl;
    cout<<"ADS: lattice with non-local interacton via AdS space"<<endl;
    exit(0);
  }
  
  MaxIter    = atoi(argv[5]);
  tol        = atof(argv[6]);
  t          = atoi(argv[7]);
  msqr       = atof(argv[8]);
  delta_msqr = atof(argv[9]);
  Levels     = atoi(argv[10]);
  src_pos    = atoi(argv[11]);
  
  if(atof(argv[12]) == 0) {
    if(t > 1) C_msqr = (1.57557326 + 1.56565549/msqr);
    else C_msqr = -0.0126762/msqr + 0.0689398*msqr + 2.02509;
  }
  else C_msqr = atof(argv[12]);
  
  if(atof(argv[13]) == 0) N_latt = 0.294452/(msqr + 0.766901) + 0.0788137;
  else N_latt = atof(argv[13]);
  
  q = atoi(argv[14]);
  n_shift = atoi(argv[15]);
  
  //MC params  
  n_therm = atoi(argv[16]);
  n_meas  = atoi(argv[17]);
  n_skip  = atoi(argv[18]);
  n_wolff = atoi(argv[19]);
  musqr   = atof(argv[20]);
  lambda  = atof(argv[21]);
  sigma   = atof(argv[22]);

  t_weight_scale  = atof(argv[23]);  
  
}

void Param::print() {
  cout<<"Parameter status:"<<endl;
  cout<<"Triangulation = "<<q<<endl;
  cout<<"B.C. = "<< (bc ? ("Dirichlet") : ("Neumann") ) << endl;
  cout<<"Centre = "<< (Vcentre ? ("Vertex") : ("Circum") ) << endl;
  cout<<"Source Position = "<<src_pos<<endl;
  cout<<"CG MaxIter = "<<MaxIter<<endl;
  cout<<"CG Tol = "<<tol<<endl;
  cout<<"MCG Number of Shifts = "<<n_shift<<endl;
  cout<<"MCG Msqr increment = "<<delta_msqr<<endl;
  cout<<"TimeSlices = "<<t<<endl;
  cout<<"Mass squared = "<<msqr<<endl;
  cout<<"Levels = "<<Levels<<endl;
  cout<<"Mass squared scaling = "<<C_msqr<<endl;
  cout<<"Lattice normalisation = "<<N_latt<<endl;

}

//Using the formula c(n) = (q-4)*c(n-1) - c(n-2) where c is the
//number of nodes on circumference at level n, we can construct
//the address of the end node on a given level for triangulation q:
//EN(lev,q) = SUM c(n) n=0..level
long unsigned int endNode(int lev, Param &P) { 
  
  int q = P.q;
  
  //This is the case where a vertex is at the centre of the disk
  if(P.Vcentre == true) {
    //Explicit results for level <= 2.
    if(lev==0) return 0;
    if(lev==1) return q;
    if(lev==2) return q*(q-4) + q;
    
    //level >= 3
    int p=q-4; //A convenience  
    long unsigned int c[20]; //Array to hold the circumference info
    
    //initialise array
    for(int i=0; i<20; i++) c[i] = 0;
    //manually set 0,1,2 entries:
    c[0] = 0;
    c[1] = q;
    c[2] = p*q;
    
    long unsigned int EN = 0; //address of node at end of level;  
    //Loop up to lev
    for(int n=3; n<lev+1; n++) {
      c[n] = p*c[n-1] - c[n-2];
      EN += c[n];
    }
    EN += (c[1] + c[2]);
    return EN;
  }

  //This is the case where a circumcentre is at the centre of the disk.
  else {
    //Explicit results for level <= 1.
    if(lev==0) return 2;
    if(lev==1) return (q-3)*3 + 2;
    
    //level >= 2
    int p=q-4; //A convenience  
    long unsigned int c[20]; //Array to hold the circumference info
    
    //initialise array
    for(int i=0; i<20; i++) c[i] = 0;
    //manually set 0,1,2 entries:
    //NB!!! There are three nodes on level 0, but the END NODE is addressed
    //as 2. Therefore, when correcting the number of nodes on a
    //circumference, we use 3 (there are three nodes)
    //but when giving the endNode count, we use 2 (we count from 0)
    c[0] = 3;       
    c[1] = (q-3)*3;
    
    long unsigned int EN = 0; //address of node at end of level;  
    //Loop up to lev
    for(int n=2; n<lev+1; n++) {
      c[n] = p*c[n-1] - c[n-2];
      EN += c[n];
    }
    EN += (2 + c[1]); //0,1,2 nodes to add from circumference 0
    return EN;
  }
}


//Print the hyperbolic and poincare radii
void checkRadius(vector<Vertex> &NodeList, Param P){

  double hyp_rad_ave = 0.0;
  double poi_rad_ave = 0.0;
  double hyp_rad_sig = 0.0;
  double poi_rad_sig = 0.0;

  int q = P.q;
  int count = 0;
  int TotNumber = (endNode(P.Levels,P)+1);
  
  for(long unsigned int n=0; n<TotNumber; n++) {
    //loop over neighbours. Break if an outer node is found.
    for(int m=0; m<q; m++) {
      if(NodeList[n].nn[m] == -1 && NodeList[n].pos != -1) {
	//cout<<"n = "<<n<<" |z| = "<<abs(NodeList[n].z)<<" s = ";
	//cout<<s(NodeList[n].z)<<endl;
	
	poi_rad_ave += abs(NodeList[n].z);
	hyp_rad_ave += s(NodeList[n].z);
	count++;
	m=q;
      }
    }
  }
  
  hyp_rad_ave /= count;
  poi_rad_ave /= count;
  
  for(long unsigned int n=0; n<TotNumber; n++) {
    for(int m=0; m<q; m++) {
      if(NodeList[n].nn[m] == -1 && NodeList[n].pos != -1) {
	hyp_rad_sig += pow(hyp_rad_ave - s(NodeList[n].z),2);
	poi_rad_sig += pow(poi_rad_ave - abs(NodeList[n].z),2);
	m=q;
      }
    }
  }

  hyp_rad_sig /= (count-1);
  poi_rad_sig /= (count-1);
  hyp_rad_sig = sqrt(hyp_rad_sig);
  poi_rad_sig = sqrt(poi_rad_sig);
  
  cout<<"HYP RAD AVE = "<<hyp_rad_ave<<endl;
  cout<<"HYP RAD SIG = "<<hyp_rad_sig<<endl;
  cout<<"POI RAD AVE = "<<poi_rad_ave<<endl;
  cout<<"POI RAD SIG = "<<poi_rad_sig<<endl;

  
  //Eyeball the output. Something out of place will
  //stick out like a sore thumb.
  //if(P.verbosity) PrintNodeTables(NodeList, P);  
}


void checkArea(const vector<Vertex> NodeList, Param P) {

  double length_01 = 0.0;
  double length_02 = 0.0;
  double length_12 = 0.0;
  double equi_area = area3q(P.q);
  double ave       = 0.0;

  double sig1 = 0.0;
  double sig2 = 0.0;
  int count = 0;

  if(P.verbosity) cout<<endl<<"Checking boundary areas"<<endl;
  for(long unsigned int n=endNode(P.Levels-2,P) + 1; n<endNode(P.Levels-1,P)+1; n++) {
    for(int k=0; k<NodeList[n].fwdLinks; k++) {
      length_01 = d12(NodeList[n].z, NodeList[NodeList[n].nn[k+1]].z);
      length_02 = d12(NodeList[n].z, NodeList[NodeList[n].nn[k+2]].z);
      length_12 = d12(NodeList[NodeList[n].nn[k+1]].z, NodeList[NodeList[n].nn[k+2]].z);      
      ave += areaGeneral(P, length_01, length_02, length_12);
      count++;
    }
  }
  
  ave /= count;
  
  for(long unsigned int n=endNode(P.Levels-2,P) + 1; n<endNode(P.Levels-1,P)+1; n++) {
    for(int k=0; k<NodeList[n].fwdLinks; k++) {
      length_01 = d12(NodeList[n].z, NodeList[NodeList[n].nn[k+1]].z);
      length_02 = d12(NodeList[n].z, NodeList[NodeList[n].nn[k+2]].z);
      length_12 = d12(NodeList[NodeList[n].nn[k+1]].z, NodeList[NodeList[n].nn[k+2]].z);      
      if(P.verbosity) cout<<"n="<<n<<" area "<<k+1<<" = "<<areaGeneral(P, length_01, length_02, length_12)<<endl;
      sig1 += pow(equi_area - areaGeneral(P, length_01, length_02, length_12),2);
      sig2 += pow(ave       - areaGeneral(P, length_01, length_02, length_12),2);
    }
  }

  sig1 /= count - 1;
  sig2 /= count - 1;
  sig1 = sqrt(sig1);
  sig2 = sqrt(sig2);

  cout<<"Boundary areas"<<endl;
  cout<<"AREA EQUI = "<<equi_area<<endl;
  cout<<"AREA STD DEV W.R.T. EQUI = "<<sig1<<endl;
  cout<<"AREA AVE = "<<ave<<endl;  
  cout<<"AREA STD DEV W.R.T AVE = "<<sig2<<endl;  

}

void checkEdgeLength(const vector<Vertex> NodeList, Param P) {
  
  int q = P.q;
  int Levels = P.Levels;
  double length = 0.0;
  double sig = 0.0;
  int  nn_node;
  bool Vcentre = P.Vcentre;
  double length_0 = d12(NodeList[0].z, NodeList[1].z);
  double tol = 1e-2;

  //Level 0 is specific to how the graph is centred.
  if(Vcentre) {
    if(P.verbosity) cout<<" lev =  " << 0 << endl;
    if(P.verbosity) cout<<endl<<" Node number = "<<0<<" : "<<endl;
    for(int i = 0; i < q; i++){
      nn_node = NodeList[0].nn[i];
      length = d12(NodeList[0].z, NodeList[nn_node].z);
      if(P.verbosity) cout<<" "<<NodeList[0].nn[i]<<" > "<<length<<"  ";
    }
  }
  else {
    if(P.verbosity) cout<<" lev = "<<0<<endl;
    if(P.verbosity) cout<<endl<<" Node number = "<<0<<" : "<<endl;
    for(int i = 0; i < 2; i++){
      nn_node = NodeList[0].nn[i];
      length = d12(NodeList[0].z, NodeList[nn_node].z);
      if(P.verbosity) cout << NodeList[0].nn[i] << " >  " << length<< "  ";
    }
  }
  
  for(int lev = 1; lev < Levels+1; lev++)  {
    if(P.verbosity) cout<<endl<<endl<<" lev = "<<lev<<endl;      
    for(long unsigned int n = endNode(lev-1,P) + 1;n < endNode(lev,P) + 1 ;n++) {
      if(P.verbosity) cout<<endl<<" Node number = "<<n<<":"<<endl;
      sig += pow( length_0 - d12(NodeList[n].z, NodeList[NodeList[n].nn[q-1]].z), 2);
      
      for(int i = 0; i <q; i++){
	nn_node = NodeList[n].nn[i];
	if(NodeList[n].nn[i] != -1 ) {
	  length = d12(NodeList[n].z, NodeList[nn_node].z);
	  if(P.verbosity) {
	    cout<<" to "<<NodeList[n].nn[i]<<" = "<<length<<" ";
	    if(abs(length - length_0)/length_0 > tol) cout<<"<-! "<<endl;
	    else cout<<"    "<<endl;
	  }
	}
      }
    }
  }
  sig /= endNode(Levels,P);
  sig = sqrt(sig);
  cout<<endl<<"LENGTH STD DEV = "<<sig<<endl;
  if(sig>tol) {
    cout<<"WARNING: Hypergeometric length STD_DEV has diverged over "<<tol<<endl;
    //exit(0);
  }
}

//inline double edgeLength(int q) {
//return sqrt( 1 - 4*sin(M_PI/q)*sin(M_PI/q) );
//}

