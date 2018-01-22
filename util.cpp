#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>

#include "util.h"

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

  std::string latType(argv[4]);
  if ( (latType == "sq" || latType == "SQ") || latType == "Sq"){
    lattice = true;
  } else if(latType == "ads" || latType == "AdS") {
    lattice = false;
  } else {
    cout<<"Invalid Lattice type given. Please use sq/AdS for square/AdS lattices."<<endl;
    exit(0);
  }
  
  MaxIter    = atoi(argv[5]);
  tol        = atof(argv[6]);
  t          = atoi(argv[7]);
  msqr       = atof(argv[8]);
  delta_msqr = atof(argv[9]);
  Levels     = atoi(argv[10]);
  src_pos    = atoi(argv[11]);
  
  //if(atof(argv[11]) == 0) C_msqr = -0.0126762/msqr + 0.0689398*msqr + 2.02509;
  if(atof(argv[12]) == 0) {
    if(t > 1) C_msqr = (1.57557326 + 1.56565549/msqr);
    else C_msqr = -0.0126762/msqr + 0.0689398*msqr + 2.02509;
  }
  else C_msqr = atof(argv[12]);
  
  //if(atof(argv[12]) == 0) N_latt = 0.294452/(msqr + 0.766901) + 0.0788137;
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
long unsigned int endNode(int lev, Param P) { 
  
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

//- Get the z coordinates of every node on the Poincare disk 
void GetComplexPositions(std::vector<Vertex> &NodeList, Param& P){

  int q = P.q;
  int Levels = P.Levels;
  int T_offset = endNode(P.Levels,P)+1;
  
  if(P.Vcentre == true) {
    //Assume for now that the origin (level 0) is a vertex
    NodeList[0].z = 0.0;
    //Assert that node 1 is on the real axis
    complex<double> init(edgeLength(q),0.0);
    NodeList[1].z = init;
    //Rotate to create level level 1
    for(int k=1; k<q+1; k++) {
      NodeList[k].z = newVertex(init, 0.0, k-1, q);
    }
    //For every node on level >=1, the zeroth
    //nn is the (n-1)th link on the same level. If
    //n=1, the node address is the endnode value.
    for(int l=1; l<Levels+1; l++) {
      for(long unsigned int n=endNode(l-1,P)+1; n<endNode(l,P)+1; n++) {
	for(int k=0; k<q; k++) {
	  if(NodeList[n].nn[k] != -1) {
	    NodeList[NodeList[n].nn[k]].z = newVertex(NodeList[NodeList[n].nn[0]].z, NodeList[n].z, k, q);
	    NodeList[NodeList[n].nn[k]].temporal_weight = 1.0/ ((1+pow(abs(NodeList[NodeList[n].nn[k]].z),2))/(1-pow(abs(NodeList[NodeList[n].nn[k]].z),2)));
	  }
	}
      }
    }    
  }
  else {
    
    double numer = sqrt(cos(M_PI*(q+6)/(6*q)) - sin(M_PI/q));
    double denom = sqrt(sin(M_PI/q) + sin(M_PI*(q+3)/(3*q)));    
    double init_mod = sqrt(norm(numer/denom));
    
    //Assume that node 0 lies on +ve real axis
    complex<double> init_0(init_mod,0.0);
    NodeList[0].z = init_0;
    
    //Assert that node 1 is node 0 rotated by 2*PI/3
    complex<double>init_1(init_mod*cos(2.0*M_PI/3.0),
			 init_mod*sin(2.0*M_PI/3.0));
    
    NodeList[1].z = init_1;
    //Rotate node 1 about node 0 to create level 0 (the equilateral triangle)
    NodeList[2].z = newVertex(init_1, init_0, 1, q);

    //For every node on level >=1, the zeroth
    //nn is the (n-1)th link on the same level. If
    //n=1, the node address is the endnode value.
    for(long unsigned int n=0; n<endNode(1,P)+1; n++) {
      for(int k=0; k<q; k++) {
	if(NodeList[n].nn[k] != -1) {
	  NodeList[NodeList[n].nn[k]].z = newVertex(NodeList[NodeList[n].nn[0]].z, NodeList[n].z, k, q);
	}
      }
    }
    for(int l=1; l<Levels+1; l++) {
      for(long unsigned int n=endNode(l-1,P)+1; n<endNode(l,P)+1; n++) {
	for(int k=0; k<q; k++) {
	  if(NodeList[n].nn[k] != -1) {
	    NodeList[NodeList[n].nn[k]].z = newVertex(NodeList[NodeList[n].nn[0]].z, NodeList[n].z, k, q);
	    NodeList[NodeList[n].nn[k]].temporal_weight = 1.0/ ((1+pow(abs(NodeList[NodeList[n].nn[k]].z),2)) / (1-pow(abs(NodeList[NodeList[n].nn[k]].z),2)));
	  }
	}
      }
    }
  }

  if(P.t > 1) {
    //Copy all 2D complex positions and weights along the cylinder
    for(long unsigned int n=0; n<endNode(P.Levels,P)+1; n++) 
      for(int t=1; t<P.t; t++) {
	NodeList[n + T_offset*t].z = NodeList[n].z;
	NodeList[n + T_offset*t].temporal_weight = NodeList[n].temporal_weight;
      }
  }
}

//- For each node n, with a link to another node,
//  it checks that the neighbour table on the linked
//  node contains the original node n as a neighbour.
void connectivityCheck(vector<Vertex> &NodeList, Param P){

  int q = P.q;
  int Levels = P.Levels;
  int T = P.t;
  int TotNumber = T*(endNode(Levels,P)+1);
  int t_offset  = 0;
  T == 1 ? t_offset = 0 : t_offset = 2;
  
  //Object to hold boolean values of graph connectivity.
  vector<Vertex> AuxNodeList(TotNumber);
  //Initialise to 0.
  for(int n = 0; n <TotNumber;n++)
    for(int mu = 0; mu < q+t_offset; mu++) {
      AuxNodeList[n].nn[mu] = 0;
    }
  
  for(long unsigned int n=0; n<TotNumber; n++) {
    //Check that the node is valid
    if(NodeList[n].pos != -1) {      
      for(int m=0; m<q+t_offset; m++) {
	//Check that the link is valid
	if(NodeList[n].nn[m] != -1) {
	  for(int p=0; p<q+t_offset; p++) {
	    //Loop over all links on the linked node,
	    //check if original node exists in neighbour
	    //table.
	    if( n == NodeList[ NodeList[n].nn[m] ].nn[p] ) {
	      AuxNodeList[n].nn[m] = 1;
	    }
	  }
	}
      }
    }
  }

  //Eyeball the output. something out of place will
  //stick out like a sore thumb.
  PrintNodeTables(AuxNodeList, P);
}

//Truncate the graph according to the hyperbolic radius
//condition |z| < s.
void hypRadGraph(vector<Vertex> &NodeList, Param &P){

  int q = P.q;
  int Levels = P.Levels;
  int T = P.t;
  int TotNumber = T*(endNode(Levels,P)+1);
  int t_offset  = 0;
  T == 1 ? t_offset = 0 : t_offset = 2;

  //Find radius to maximize connectivity.
  int r_min_pos = 0;
  int r_max_pos = 0;
  double r_min = 1.0;
  double r_max = 0.0;

  //Locate the node on the outer circumference with the smallest
  //radius
  for(int n=endNode(P.Levels-1,P)+1; n<endNode(P.Levels,P)+1; n++){
    if(abs(NodeList[n].z) < r_min) {
      r_min = abs(NodeList[n].z);
      r_min_pos = n;
    }
  }
  P.hyp_rad = s(NodeList[r_min_pos].z);
  cout<<"HYP_RAD = "<<P.hyp_rad<<endl;

  double hyp_rad = P.hyp_rad;
  
  for(long unsigned int n=0; n<TotNumber; n++) {
    if(s(NodeList[n].z) >= hyp_rad + 0.0) {     
      //This node must be removed. loop over its neighbours
      //and remove this specific connection.
      if(P.verbosity) cout<<"Deletng node: "<<n<<" connections: ";
      NodeList[n].pos = -1;
      for(int m=0; m<q+t_offset; m++)
	for(int k=0; k<q+t_offset; k++) {
	  if(NodeList[NodeList[n].nn[m]].nn[k] == n) {
	    NodeList[NodeList[n].nn[m]].nn[k] = -1;
	    NodeList[n].nn[m] = -1;
	    if(P.verbosity) cout<<NodeList[n].nn[m]<<" ";
	  }
	}
      if(P.verbosity) cout<<endl;
    }
  }
  
  for(int n=endNode(P.Levels-2,P)+1; n<endNode(P.Levels-1,P)+1; n++){
    if(abs(NodeList[n].z) < r_min && NodeList[n].pos != -1) {
      r_min = abs(NodeList[n].z);
      r_min_pos = n;
    }
    if(abs(NodeList[n].z) > r_max && NodeList[n].pos != -1) {
      r_max = abs(NodeList[n].z);
      r_max_pos = n;
    }
  }
  P.r_max_pos = r_max_pos;
  cout<<"R MAX POS = "<<P.r_max_pos<<endl;
  P.r_min_pos = r_min_pos;
  cout<<"R MIN POS = "<<P.r_min_pos<<endl;
  
  //Eyeball the output. Something out of place will
  //stick out like a sore thumb.
  //if(P.verbosity == "d") radiusCheck(NodeList, P); 
}

//Print the hyperbolic and poincare radii
void radiusCheck(vector<Vertex> &NodeList, Param P){

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

void PrintNodeTables(const vector<Vertex> NodeList, Param P) {

  int q = P.q;
  int Levels = P.Levels;  
  int T = P.t;
  int t_offset  = 0;
  T == 1 ? t_offset = 0 : t_offset = 2;
  
  for(int t=0; t<1; t++) {
    int offset = t*( endNode(Levels,P) + 1 );    
    cout<<endl<<"lev = "<<0<<" T = "<<t<<endl;
    
    if(P.Vcentre) {
      cout<<endl<<" Node number = "<<0 + offset<<" : ";
      for(int i = 0; i < q+t_offset; i++) cout<<NodeList[offset + 0].nn[i]<<"  ";
    }      
    else {
      for(long unsigned int n = 0; n < endNode(0,P)+1; n++) {
	cout<<endl<<"n= "<<n<<" Node number = "<<n + offset<<" FL="<<NodeList[n].fwdLinks<<" |z|= "<<abs(NodeList[n].z)<<" s= "<<s(NodeList[n].z)<<" : ";
	for(int i = 0; i < q+t_offset; i++) cout<<NodeList[offset + n].nn[i]<< "  ";
      }
    }
    for(int lev = Levels-1; lev < Levels+1; lev++)  {
      cout<<endl<<"lev = "<<lev<<" T = "<<t<<endl;
      long unsigned int low = endNode(lev-1,P)+1;
      long unsigned int high= endNode(lev,P)+1;
      
      for(long unsigned int n = low; n < high; n++) {
	cout<<endl<<"n= "<<n<<" Node number = "<<NodeList[n + offset].pos<<" FL="<<NodeList[n].fwdLinks<<" |z|= "<<abs(NodeList[n].z)<<" s= "<<s(NodeList[n].z)<<" : ";
	for(int i = 0; i < q+t_offset; i++)
	  cout << NodeList[offset + n].nn[i] << "  ";
      }
    }      
  }  
  cout<<endl;
}

void PrintComplexPositions(const vector<Vertex> NodeList, Param P) {

  int Levels = P.Levels;
  
  cout<<endl<<"#Printing for Level 0"<<endl;
  for(long unsigned int n=0; n<endNode(0,P)+1; n++) {
    //cout<<"n="<<n<<" z="<<NodeList[n].z.real()<<","<<NodeList[n].z.imag();
    cout<<"n= "<<n<<" |z|="<<abs(NodeList[n].z)<<" s="<<s(NodeList[n].z);    
  }
  for(int l=1; l<Levels+1; l++) {
    cout<<endl<<"Printing for Level "<<l<<endl;
    for(long unsigned int n=endNode(l-1,P)+1; n<endNode(l,P)+1; n++) {
      //cout<<"n="<<n<<" z="<<NodeList[n].z.real()<<","<<NodeList[n].z.imag();
      cout<<"n= "<<n<<" |z|="<<abs(NodeList[n].z)<<" s="<<s(NodeList[n].z);
      cout<<endl;
    }
  }
}


void CheckArea(const vector<Vertex> NodeList, Param P) {

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

void CheckEdgeLength(const vector<Vertex> NodeList, Param P) {
  
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

//Data file for lattice/analytical propagator data,
void DataDump(vector<Vertex> NodeList, double *phi, Param p, int level, int t_range, int shift) {

  long unsigned int TotNumber = (endNode(p.Levels,p) + 1) * p.t;
  double norm = 0.0;
  for(long unsigned int i = 0;i < TotNumber; i++) norm += phi[i]*phi[i];
  for(long unsigned int i = 0;i < TotNumber; i++) phi[i] /= sqrt(norm); 
  
  long unsigned int j = p.src_pos;

  int T = p.t;
  int T_offset = 0;
  double theta = 0.0;
  //double msqr = (double)p.msqr + (double)p.delta_msqr*shift;
  double delta = p.t < 2 ? 0.5 + sqrt(0.25 + p.msqr) : 1.0 + sqrt(1 + p.msqr);
  complex<double> ratio;
  complex<double> src = NodeList[j].z;
  
  //Loop over timeslices
  for(int t=0; t<T/2; t++) {
    T_offset = (endNode(p.Levels,p) + 1) * t;

    //Loop over circumference levels
    for(int lev=0; lev<p.Levels; lev++) {

      if(lev+1 == level && t<t_range) {
	sprintf(p.fname, "./data_dump/q%d_Lev%d_T%d_BASEmsqr%.5e_LATTmsqr%.5e_srct0_srcpos%d_sinkt%dLev%d_%s_%s.dat",
		p.q,
		p.Levels,
		p.t,
		(double)p.msqr,
		(double)(p.C_msqr*p.msqr + shift*p.delta_msqr),
		p.src_pos,
		0,
		lev+1,
		p.bc == true ? "Dirichlet" : "Neumann",
		p.Vcentre == true ? "Vertex" : "Circum");
	FILE *fp1;
	fp1=fopen(p.fname, "a");
	
	//Loop over H2 disk
	for(long unsigned int k = endNode(lev,p)+1; k < endNode(lev+1,p)+1; k++) {
	  
	  //Construct i
	  int i = k + T_offset;
	  
	  ratio = NodeList[i].z/NodeList[j].z;
	  theta = atan2( ratio.imag() , ratio.real() );
	  complex<double> snk = NodeList[i].z;
	  
	  //index divided by disk size, using the int floor feature/bug,
	  //gives the timeslice for each index.
	  int t1 = j / (TotNumber/p.t);
	  int t2 = i / (TotNumber/p.t);
	  //Assume PBC.
	  int delta_t = (t2-t1);// > p.t/2 ? (t2-t1) - p.t : (t2-t1);
	  
	  double r = abs(NodeList[i].z);
	  double r_p = abs(NodeList[j].z);
	  double xi = (cosh(delta_t)*(1+r)*(1+r_p) - 4*r*r_p*cos(theta)) / ((1-r)*(1-r_p));
	  
	  if( i != j && (abs(theta) > M_PI/2) ) {
	    fprintf(fp1, "%d %d %.8e %.8e %.8e %.8e %.8e\n",
		    //1 Timeslice, 2 H2 pos
		    t, (int)k,
		    
		    //3 source/sink angle
		    (double)theta,
		    
		    //4 lattice prop
		    (double)p.N_latt*(double)phi[i],
		    
		    //5 invariant
		    (double)1.0/((double)xi),
		    //(double)( ((double)1.0-abs(snk))*((double)1.0-abs(src))/(pow(abs(snk - src),2))),
		    
		    //6 AdS2p1 formula
		    (double)(exp(-delta*sigma(src,snk,delta_t)) / (1 - exp(-2*sigma(src,snk,delta_t)))),
		    
		    //7 geodesic
		    (double)sigma(src,snk,delta_t)
		    );
	  }
	}
	fclose(fp1);
      }
    }
  }
}

void visualiserSqr(vector<double> phi_cyl, double barr, Param p) {  
  
  for(int i=0; i<p.S1; i++) {
    for(int j=0; j<p.Lt; j++) {
      if(phi_cyl[i + p.AdSVol*j] < -1.5*barr) cout<<"\033[1;41m \033[0m";
      if(-1.5*barr < phi_cyl[i + p.AdSVol*j] && phi_cyl[i + p.AdSVol*j] < -0.5*barr) cout<<"\033[1;43m \033[0m";
      if(-0.5*barr < phi_cyl[i + p.AdSVol*j] && phi_cyl[i + p.AdSVol*j] <  0.5*barr) cout<<"\033[1;42m \033[0m";
      if( 0.5*barr < phi_cyl[i + p.AdSVol*j] && phi_cyl[i + p.AdSVol*j] <  1.5*barr) cout<<"\033[1;46m \033[0m";
      if( 1.5*barr < phi_cyl[i + p.AdSVol*j]) cout<<"\033[1;44m \033[0m";
    }
    cout<<endl;
  }
  //usleep(250000);
}

void visualiserAdS(vector<Vertex> NodeList, double barr, Param p) {  

  cout<<endl<<"Level "<<0<<endl;
  for(int j=0; j<p.Lt; j++) {
    int i = 0;
    if(NodeList[i + (p.AdSVol)*j].phi <= -1.5*barr) cout<<"\033[1;41m \033[0m";
    if(-1.5*barr < NodeList[i + p.AdSVol*j].phi && NodeList[i + p.AdSVol*j].phi <= -0.5*barr) cout<<"\033[1;43m \033[0m";
    if(-0.5*barr < NodeList[i + p.AdSVol*j].phi && NodeList[i + p.AdSVol*j].phi <=  0.5*barr) cout<<"\033[1;42m \033[0m";
    if( 0.5*barr < NodeList[i + p.AdSVol*j].phi && NodeList[i + p.AdSVol*j].phi <=  1.5*barr) cout<<"\033[1;46m \033[0m";
    if( 1.5*barr < NodeList[i + p.AdSVol*j].phi) cout<<"\033[1;44m \033[0m";
  }
  cout<<endl;
  for(int l=1; l<p.Levels+1; l++) {
    cout<<"Level "<<l<<endl;
    for(int i=endNode(l-1,p)+1; i<endNode(l,p)+1; i++) {
      for(int j=0; j<p.Lt; j++) {
	if(NodeList[i + (p.AdSVol)*j].phi <= -1.5*barr) cout<<"\033[1;41m \033[0m";
	if(-1.5*barr < NodeList[i + p.AdSVol*j].phi && NodeList[i + p.AdSVol*j].phi <= -0.5*barr) cout<<"\033[1;43m \033[0m";
	if(-0.5*barr < NodeList[i + p.AdSVol*j].phi && NodeList[i + p.AdSVol*j].phi <=  0.5*barr) cout<<"\033[1;42m \033[0m";
	if( 0.5*barr < NodeList[i + p.AdSVol*j].phi && NodeList[i + p.AdSVol*j].phi <=  1.5*barr) cout<<"\033[1;46m \033[0m";
	if( 1.5*barr < NodeList[i + p.AdSVol*j].phi) cout<<"\033[1;44m \033[0m";
      }
      cout<<endl;
    }
  }
}

void visualiserPhi2(double **phi_sq, Param p, int iter) {
  
  double barr_t[p.S1];
  double barr = 0.0;
  
  for(int i=0; i<p.S1; i++) {
    barr_t[i] = 0.0;
    for(int j=0; j<p.Lt; j++) {
      barr_t[i] += phi_sq[i][j];
    }
    barr_t[i] /= p.Lt;
    barr      += barr_t[i];
  }
  barr /= p.S1;
  
  for(int i=0; i<p.S1; i++) {
    for(int j=0; j<p.Lt; j++) {
      if( phi_sq[i][j] < barr ) cout<<"\033[1;41m \033[0m";
      if( barr < phi_sq[i][j] ) cout<<"\033[1;44m \033[0m";
    }
    cout<<" phi_sq_ave[s="<<i<<"] = "<<barr_t[i]/(iter+1)<<endl;
  }
}


void correlatorsAdS(double **corr, double **corr_ave, int corr_norm,
		    vector<Vertex> NodeList, double avePhi, Param p) {
  int s_idx = 0;
  int t_idx = 0;

  int s_size = p.S1;
  int t_size = p.Lt;

  int offset = endNode(p.Levels-1,p)+1;
  int disk   = p.AdSVol;
  
  //Each value of \delta t and \delta \theta
  //must be normalised seperately.
  int norm[s_size/2][t_size/2];
  for(int i=0; i<s_size/2; i++)
    for(int j=0; j<t_size/2; j++) {
      norm[i][j] = 0;
      corr[i][j] = 0.0;
    }
  
  //loop over sink/source *theta*
  for(int is=0; is<s_size; is++)
    for(int js=0; js<s_size; js++) {
      
      s_idx = abs(is-js);
      if(s_idx != s_size/2) {
	if(s_idx >= s_size/2) s_idx = s_size - s_idx - 1;
      
	//loop over sink/source *temporal*
	for(int il=0; il<t_size; il++) 
	  for(int jl=0; jl<t_size; jl++) {
	    
	    t_idx = abs(il-jl);
	    if(t_idx != t_size/2) {
	      if(t_idx >= t_size/2) t_idx = t_size - t_idx - 1;
	      
	      corr[s_idx][t_idx] += ((NodeList[offset + is + il*disk].phi - avePhi) *
				     (NodeList[offset + js + jl*disk].phi - avePhi));
	      norm[s_idx][t_idx]++;
	    }
	  }
      }
    }
  
  //Normalise, add to running average
  for(int i=0; i<s_size/2; i++)
    for(int j=0; j<t_size/2; j++) {
      corr[i][j] /= norm[i][j];
      corr_ave[i][j] += corr[i][j];
    }

  //Corr dump to stdout and collect running average.
  if(p.verbosity) cout<<setprecision(4);      
  if(p.verbosity) cout<<"Corr Sample: "<<endl;
  for(int i=0; i<s_size/2; i++) {
    if(p.verbosity) cout<<"theta "<<2*M_PI*i/(double)(s_size)<<":";
    for(int j=0; j<t_size/2; j++) {
      if(p.verbosity) cout<<" "<<corr_ave[i][j]/corr_norm;
      corr[i][j] = corr_ave[i][j]/corr_norm;
    }
    if(p.verbosity) cout<<endl;
  }

  if(!p.verbosity) {
    cout<<setprecision(4);      
    cout<<"Corr Sample:";
    for(int j=0; j<t_size/2; j++) {
      cout<<" "<<corr[0][j];
    }
    cout<<endl;
  }
  
  //**corr now contains the current, normalised, running average correlation
  //matrix.  
}


void correlatorsSqr(double **corr, double **corr_ave, int corr_norm,
		    vector<double> phi, double avePhi, Param p) {
  
  int s_idx = 0;
  int t_idx = 0;

  int s_size = p.S1;
  int t_size = p.Lt;

  //Each value of \delta t and \delta \theta
  //must be normalised seperately.
  int norm[s_size/2][t_size/2];
  for(int i=0; i<s_size/2; i++)
    for(int j=0; j<t_size/2; j++) {
      norm[i][j] = 0;
      corr[i][j] = 0.0;
    }
  
  //loop over sink/source *theta*
  for(int is=0; is<s_size; is++)
    for(int js=0; js<s_size; js++) {
      
      s_idx = abs(is-js);
      if(s_idx != s_size/2) {
	if(s_idx >= s_size/2) s_idx = s_size - s_idx - 1;
      
	//loop over sink/source *temporal*
	for(int il=0; il<t_size; il++) 
	  for(int jl=0; jl<t_size; jl++) {
	    
	    t_idx = abs(il-jl);
	    if(t_idx != t_size/2) {
	      if(t_idx >= t_size/2) t_idx = t_size - t_idx - 1;
	      
	      corr[s_idx][t_idx] += ((phi[is + p.S1*il] - avePhi) *
				     (phi[js + p.S1*jl] - avePhi));
	      norm[s_idx][t_idx]++;
	    }
	  }
      }
    }
  
  //Normalise, add to running average
  for(int i=0; i<s_size/2; i++)
    for(int j=0; j<t_size/2; j++) {
      corr[i][j] /= norm[i][j];
      corr_ave[i][j] += corr[i][j];
    }

  //Corr dump to stdout and collect running average.
  if(p.verbosity) cout<<setprecision(4);      
  if(p.verbosity) cout<<"Corr Sample: "<<endl;
  for(int i=0; i<s_size/2; i++) {
    if(p.verbosity) cout<<"theta "<<2*M_PI*i/(double)(s_size)<<":";
    for(int j=0; j<t_size/2; j++) {
      if(p.verbosity) cout<<" "<<corr_ave[i][j]/corr_norm;
      corr[i][j] = corr_ave[i][j]/corr_norm;
    }
    if(p.verbosity) cout<<endl;
  }

  if(!p.verbosity) {
    cout<<setprecision(4);      
    cout<<"Corr Sample:";
    for(int j=0; j<t_size/2; j++) {
      cout<<" "<<corr[0][j];
    }
    cout<<endl;
  }
  
  //**corr now contains the current, normalised, running average correlation
  //matrix.  
}

void autocorrelation(double *PhiAb_arr, double avePhiAbs, int meas) {
  
}


// Translate w to 0 
complex<double> T(complex<double> z,  complex<double> w){ //translate w to 0
  return (z - w)/(z*conj(w) + (double)1.0);
}

// Rotate at z = 0
complex<double> R(complex<double> z, complex<double> omega){
  // rotate by omega = exp [i theta] about z= 0
  return omega*z;
 }

//Reflection z accross the z1 -- z2 line
complex<double> flip(complex<double> z, complex<double> z1, complex<double> z2){
  // reflection (or flip)  z across (z1,z2)
  return T( R( conj( T(z,z1) ), z2/conj(z2) ), -z1 );
}

//Hyperbolic distance s, from origin to z
double  s(complex<double> z){
   return log(((double)1.0+abs(z))/((double)1.0-abs(z)));
}

//Poincare distance |z| from origin to s
double r(double s){
  return tanh(s/2);
}

//Geodesic distance from z1 to z2
double d12(complex<double> z, complex<double> w) {
  return log ( (abs((double)1.0-conj(z)*w) + abs(z-w))/(abs((double)1.0-conj(z)*w) - abs(z-w)));
  
}

//Geodesic distance from z1,t1 to z2,t2
double sigma(complex<double> z, complex<double> w, int delta_t) {

  double theta = atan2( (w/z).imag() , (w/z).real() );
  double r = abs(z);
  double r_p = abs(w);  
  double xi = (cosh(delta_t)*(1+r)*(1+r_p) - 4*r*r_p*cos(theta)) / ((1-r)*(1-r_p)); 
  
  return acosh(xi);
    
}

// length of arc q fold triangle to origin.
double s3p(int q){  //vertex centered Arc lengeth
  return (double)2.0*acosh((double)1.0/sin(M_PI/(double)q));
}

// Area equilateral triangle with angles 2 pi/q
double area3q(int q){
  //pi - (3 * hyp_angle) = defect
  return M_PI - (double)3.0*(2.0*M_PI/(double)q);
}

// Area non-equilateral triangle with side length a,b,c
double areaGeneral(Param P, double a, double b, double c) {
  //pi - (A+B+C) = defect
  
  // use general cosine law:
  // cos(C) = (cosh(c) - cosh(a)cosh(b)) / sinh(a)sinh(b)
  double C = acos( -(cosh(c) - cosh(a)*cosh(b)) / (sinh(a)*sinh(b)) );
  
  // use general sine law:
  // sin(A)/sinh(a) = sin(B)/sinh(b) = ...
  double B = asin( sinh(b)*sin(C)/sinh(c) );
  double A = asin( sinh(a)*sin(C)/sinh(c) );

  return M_PI - (A+B+C);
}

//
double centralRad(double s){
  return (sqrt( cosh(s/2.0) - (double)1.0/4.0) - 0.5*sqrt(3.0))/sqrt(cosh(s/2.0) -(double)1.0);
}

//
complex<double> DisktoUHP(complex<double> z) {
  // z = -1 -i, 1, i maps to u =  -1, 0 1, infty
  return (z + I)/((double)1.0 + I * z);
}
//
complex<double> UHPtoDisk(complex<double> u) {
  // u = 0, 1, infty  maps to -1 -i , 1, i  
  return (u - I)/((double)1.0 - I*u); 
}

//- Rotate z about z0 by 2*k*pi/q 
complex<double> newVertex(complex<double> z,complex<double> z0, int k, int q) {

  complex<double> w( 0.0, 2.0 * sin(k * M_PI/q) );
  complex<double> a( cos(k*M_PI/q)*((double)1.0 - norm(z0)), sin(k*M_PI/q)*((double)1.0 + norm(z0)) ); 
  w = w*z0;
  
  //cout<<"New z = "<<-(a*z - w)/(conj(w)*z - conj(a))<<endl;
  return - (a*z - w)/(conj(w)*z - conj(a)); 
}


complex<double> inversion(complex<double> z0, double r) {
  // z_image conj(z0) = r^2
  return r*2/conj(z0);
}

complex<double> squareInversion(complex<double>z0,double r1,double r2 ) {
  return inversion(inversion(z0, r1),r2);
}

double greens2D(complex<double> z, complex<double> w) {
  return -log( tanh ( log ( (abs((double)1.0-conj(z)*w) + abs(z-w))/(abs((double)1.0-conj(z)*w) - abs(z-w)) )/2 ) );    
}

double greensM2D(complex<double> z, complex<double> w, Param p) {
  
  //First compute 2F1  

  double delta = p.t > 1 ? 1.0+sqrt(1.0+p.msqr) : 0.5+sqrt(0.25+p.msqr);
  double h = 1;
  double result = 0.0;
  double result_0 = 0.0;
  double geo = exp(-2*d12(z,w));
  double a,b,c;
  double tol = 1e-10;
  int n=0;
  bool conv = false;

  while( !conv && n < 10000 ) {    
    result_0 = result;
    a = tgamma(delta + n)/tgamma(delta);
    b = tgamma(h + n)/tgamma(h);
    c = tgamma(delta+1-h + n)/tgamma(delta+1-h);
    result += ( (a*b) / (c*tgamma(n+1)) ) * pow(geo,n);
    if( abs(result_0 - result)/result_0 < tol ) conv = true;
    n++;
    if(n%10000 == 0) cout<<n<<" 2F1 iters "<<geo<<" "<<abs(result_0 - result)/result_0<<endl; 
  }
  
  //Then compute coefficient. 
  result *= pow(geo,delta/2) * tgamma(delta) / (2*pow(M_PI,h)*tgamma(delta+1-h));

  return result;
}

double edgeLength(int q) {
  return sqrt( 1 - 4*sin(M_PI/q)*sin(M_PI/q) );
}

