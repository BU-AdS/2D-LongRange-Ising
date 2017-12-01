#ifndef UTILH
#define UTIL_H
#include <complex>
#include <cstring>

using namespace std;

#define Nlinks 7
#define I complex<double>(0.0,1.0)

class Param{

 public:

  int q = Nlinks;
  
  bool bc = false;       //if true, use Dirichlet. If false, use Neumann
  bool Vcentre = true;  //if true, place vertex at centre. If false, use circumcentre.
  bool verbosity = false;  //if true, print all data. If false, print summary.
  int MaxIter = 100000;
  double tol = pow(10,-6);
  int t = 1;
  double msqr = 0.1;
  int Levels = 3;
  int src_pos = 0;
  char fname[256];

  
  void print(){
    cout<<"Parameter status:"<<endl;
    cout<<"Triangulation = "<<q<<endl;
    cout<<"B.C. = "<< (bc ? ("Dirichlet") : ("Neumann") ) << endl;
    cout<<"Centre = "<< (Vcentre ? ("Vertex") : ("Circum") ) << endl;
    cout<<"MaxIter = "<<MaxIter<<endl;
    cout<<"Tol = "<<tol<<endl;
    cout<<"TimeSlices = "<<t<<endl;   
    cout<<"Mass squared = "<<msqr<<endl;
    cout<<"Levels = "<<Levels<<endl;
    cout<<"Source Position = "<<src_pos<<endl;
  }
  
  void init(int argc, char **argv) {
    
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
      cout<<"Invalid Verbosity conditions given. Use verbose/quiet"<<endl;
      exit(0);
    }

    MaxIter = atoi(argv[4]);
    tol     = atof(argv[5]);
    t       = atoi(argv[6]);    
    msqr    = atof(argv[7]);
    Levels  = atoi(argv[8]);
    src_pos = atoi(argv[9]);
    
  }    
};

struct Vertex{
  int nn[Nlinks+2];
  int fwdLinks;
  complex<double> z;
};

typedef vector<Vertex> Graph;


complex<double> T(complex<double> z,  complex<double> w);
complex<double> R(complex<double> z, complex<double> omega);
complex<double> flip(complex<double> z, complex<double> z1, complex<double> z2);
double s(complex<double> z);
double r(double s );
double d12(complex<double> z1, complex<double> z2);
double s3p(int q);
double area3q(int q);
double areaGeneral(Param P, double A, double B, double C);
double centralRad(double s);
complex<double> DisktoUHP(complex<double> z);
complex<double> UHPtoDisk(complex<double> u);
complex<double> inversion(complex<double> z0, double r);
complex<double> squareInversion(complex<double>z0,double r1,double r2 );
double greens2D(complex<double> z, complex<double> w);
complex<double> newVertex(complex<double> z,complex<double> z0,int k, int q);

void PrintNodeTables(const vector<Vertex> NodeList, Param P);

//- Edge length from center z = 0
double edgeLength(int q) {
  return sqrt( 1 - 4*sinl(M_PI/q)*sinl(M_PI/q) );
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
void GetComplexPositions(Graph &NodeList, Param& P){

  int q = P.q;
  int Levels = P.Levels;
  double lower = 1.0;
  double upper = 0.0;  

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
	  }
	}
	if(l==Levels) {
	  if( abs(NodeList[n].z) < lower ) lower = abs(NodeList[n].z);
	  if( abs(NodeList[n].z) > upper ) upper = abs(NodeList[n].z);
	}
      }
    }    
  }
  else {
    
    double numer = sqrt(cosl(M_PI*(q+6)/(6*q)) - sinl(M_PI/q));
    double denom = sqrt(sinl(M_PI/q) + sinl(M_PI*(q+3)/(3*q)));    
    double init_mod = sqrt(norm(numer/denom));
    
    //Assume that node 0 lies on +ve real axis
    complex<double> init_0(init_mod,0.0);
    NodeList[0].z = init_0;
    
    //Assert that node 1 is node 0 rotated by 2*PI/3
    complex<double>init_1(init_mod*cosl(2.0*M_PI/3.0),
			  init_mod*sinl(2.0*M_PI/3.0));
    
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
      //int counter = 0;
      for(long unsigned int n=endNode(l-1,P)+1; n<endNode(l,P)+1; n++) {
	for(int k=0; k<q; k++) {
	  if(NodeList[n].nn[k] != -1) {
	    NodeList[NodeList[n].nn[k]].z = newVertex(NodeList[NodeList[n].nn[0]].z, NodeList[n].z, k, q);
	  }
	}
	if(l==Levels) {
	  if( abs(NodeList[n].z) < lower ) lower = abs(NodeList[n].z);
	  if( abs(NodeList[n].z) > upper ) upper = abs(NodeList[n].z);
	}
      }
    }
  }
  
  /*
  //Histogram
  int N = 5000;
  double arr[N];
  double fac = 0.0001;
  double width = ((1+fac)*upper - (1-fac)*lower)/N;
  int count = 0;
  int distinct = 0;
  cout<<endl<<"Upper = "<<upper<<" Lower = "<<lower<<" circum = "<<endNode(Levels,P) - endNode(Levels-1,P)<<endl;
  
  for(long unsigned int n=0; n<N; n++) arr[n] = 0;    
  for(long unsigned int n=endNode(Levels-1,P)+1; n<endNode(Levels,P)+1; n++) {
    for(int i=0; i<N; i++)
      if(abs(NodeList[n].z) > (1-fac)*lower + i*width && abs(NodeList[n].z) < (1-fac)*lower + (i+1)*width ) {
	arr[i]++;
	count++;
	if(P.verbosity) cout<<abs(NodeList[n].z)<<endl;
      }
  }
  cout<<endl;
  for(int i=0; i<N; i++) {
    if(P.verbosity) cout<<arr[i]<<" ";
    if(arr[i] != 0) distinct++;
  }
  cout<<endl<<"Count = "<<count<<endl;
  cout<<endl<<"Distinct = "<<distinct<<endl;
  */

}

//- For each node n, with a link to another node,
//  it checks that the neighbour table on the linked
//  node contains the original node n as a neighbour.
void ConnectivityCheck(Graph &NodeList, Param P){

  int q = P.q;
  int Levels = P.Levels;
  int T = P.t;
  int TotNumber = T*(endNode(Levels,P)+1);
  
  //Object to hold boolean values of graph connectivity.
  vector<Vertex> AuxNodeList(TotNumber);
  //Initialise to 0.
  for(int n = 0; n <TotNumber;n++)
    for(int mu = 0; mu < q+2; mu++) {
      AuxNodeList[n].nn[mu] = 0;
    }
  
  for(long unsigned int n=0; n<TotNumber; n++) {
    for(int m=0; m<q+2; m++) {
      //Check that the link is valid
      if(NodeList[n].nn[m] != -1) {
	for(int p=0; p<q+2; p++) {
	  //Loop over all links on the linked node,
	  //check if original node exists. 
	  if( n == NodeList[ NodeList[n].nn[m] ].nn[p] ) {
	    AuxNodeList[n].nn[m] = 1;
	  }
	}
      }
    }
  }
  if(P.verbosity) PrintNodeTables(AuxNodeList, P);
}

void PrintNodeTables(const vector<Vertex> NodeList, Param P) {

  int q = P.q;
  int Levels = P.Levels;  
  int T = P.t;
  
  for(int t=0; t<T; t++) {

    int offset = t*( endNode(Levels,P) + 1 );
    
    cout << endl << "lev = " << 0 << "  T = " << t << endl;
    
    if(P.Vcentre) {
      cout << endl<< " Node number = " << 0 + offset << " : ";
      for(int i = 0; i < q+2; i++) cout << NodeList[offset + 0].nn[i] << "  ";
    }      
    else {
      for(long unsigned int n = 0; n < endNode(0,P)+1; n++) {
	cout << endl<< " Node number = " << n + offset << " FL="<<NodeList[n].fwdLinks<<" : ";
	for(int i = 0; i < q+2; i++) cout << NodeList[offset + n].nn[i] << "  ";
      } 
    }
    for(int lev = 1; lev < Levels+1; lev++)  {
      cout << endl << "lev = " << lev << "  T = " << t << endl;
      for(long unsigned int n = endNode(lev-1,P)+1; n < endNode(lev,P)+1; n++) {
	cout << endl<< " Node number = " << n + offset << " FL="<<NodeList[n].fwdLinks<<" : ";
	for(int i = 0; i < q+2; i++) cout << NodeList[offset + n].nn[i] << "  ";
      }
    }      
  }  
  cout<<endl;
}

void PrintComplexPositions(const vector<Vertex> NodeList, Param P) {

  int Levels = P.Levels;
  
  if(P.verbosity) cout<<endl<<"#Printing for Level 0"<<endl;
  for(long unsigned int n=0; n<endNode(0,P)+1; n++) {
    if(P.verbosity) {
      cout<<"n="<<n<<" z="<<NodeList[n].z.real()<<","<<NodeList[n].z.imag();
      cout<<" |z|="<<abs(NodeList[n].z)<<" phi="<<arg(NodeList[n].z);
    }
  }
  for(int l=1; l<Levels+1; l++) {
    if(P.verbosity) cout<<endl<<"Printing for Level "<<l<<endl;
    for(long unsigned int n=endNode(l-1,P)+1; n<endNode(l,P)+1; n++) {
      if(P.verbosity) {
	cout<<"n="<<n<<" z="<<NodeList[n].z.real()<<","<<NodeList[n].z.imag();
	cout<<" |z|="<<abs(NodeList[n].z)<<" phi="<<arg(NodeList[n].z);
	cout<<endl;
      }
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
    cout<<"ERROR: Hypergeometric length STD_DEV has diverged over "<<tol<<endl;
    exit(0);
  }
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


/********************************************
Basic Hyperbolic Algebra. 

Moebius Transforms 

Hyperbolic reflection for geodesic in UHP

Given Line: (x-a)^2 + y^2 = r^2 
Determined by x = a \pm r at y = 0
              x = a, y = r


 RL(a,r, u) = a + r^2/(conj(u) - a)

Preserves the line and swap a  to infty.

Map to Disc: z = D(u) =  (u -I)/(1 -I * u) with u = x + i y
Map to UHP   u = U(z) = (z + I)/(1 + I * z);

Find UHP circle that hits  +/- theta_0 on  Disc

|z - A|^2 = R^2  
|z|^2 - |z| |A| cos(theta) + |A|^2 = R^2
boundary 1 - |A| cos (theta) + |A|^2 = R^2
pick A = real. when a = 0 with map

Need 3 point to define the mobius. Circle to Circle. 

****************************************************************/

// Translate w to 0 
complex<double> T(complex<double> z,  complex<double> w)
{ //translate w to 0
  return (z - w)/(z*conj(w) + 1.0);
}

// Rotate at z = 0
complex<double> R(complex<double> z, complex<double> omega)
{
  // rotate by omega = exp [i theta] about z= 0
  return omega*z;
 }

//Reflection z accross the z1 -- z2 line
complex<double> flip(complex<double> z, complex<double> z1, complex<double> z2)
{
  // reflection (or flip)  z across (z1,z2)
  return T( R( conj( T(z,z1) ), z2/conj(z2) ), -z1 );
}

//Geodesic from z = 0 to z
double  s(complex<double> z)
{
   return logl((1.0+abs(z))/(1.0-abs(z)));
}

//Geodesic distance s from origin
double r(double s)
{
  return tanhl(s/2);
}

//Geodesic distandce from z1 to z2
double d12(complex<double> z1, complex<double> z2) {

  //2asinh( |z1 - z2| / sqrt( (1-|z1|)

  //double denom = sqrt( 1.0 - norm(z1) - norm(z2) + norm(z1)*norm(z2));
  //cout<<"Arg="<<abs(z1 - z2) / sqrt( (1.0 - norm(z1)) * (1.0 - norm(z2)));
  //cout<<" numer="<<abs(z1 - z2);
  //cout<<" denom="<<sqrt( (1.0 - norm(z1)) * (1.0 - norm(z2)))<<endl;

  double length =  2.0 * asinhl(abs(z1 - z2) / sqrt( (1.0 - norm(z1)) * (1.0 - norm(z2)) ) ) ;
  //double length =  2.0 * asinhl( abs(z1 - z2)/denom ) ;
  return length;
}

// length of arc q fold triangle to origin.
double s3p(int q)
{  //vertex centered Arc lengeth
  return 2.0*acoshl(1.0/sinl(M_PI/(double)q));
}

// Area equilateral triangle with angles 2 pi/q
double area3q(int q)
{
  //pi - (3 * hyp_angle) = defect
  return M_PI - 3.0*(2.0*M_PI/(double)q);
}

// Area non-equilateral triangle with side length a,b,c
double areaGeneral(Param P, double a, double b, double c) {
  //pi - (A+B+C) = defect
  
  // use general cosine law:
  // cos(C) = (cosh(c) - cosh(a)cosh(b)) / sinh(a)sinh(b)
  double C = acosl( -(coshl(c) - coshl(a)*coshl(b)) / (sinhl(a)*sinhl(b)) );
  
  // use general sine law:
  // sin(A)/sinh(a) = sin(B)/sinh(b) = ...
  double B = asinl( sinhl(b)*sinl(C)/sinhl(c) );
  double A = asinl( sinhl(a)*sinl(C)/sinhl(c) );

  return M_PI - (A+B+C);
}

//
double centralRad(double s)
{
  return (sqrt( cosh(s/2.0) -1.0/4.0) -0.5*sqrt(3.0))/sqrt(cosh(s/2.0) -1.0);
}

//
complex<double> DisktoUHP(complex<double> z)
{
  // z = -1 -i, 1, i maps to u =  -1, 0 1, infty
  return (z + I)/(1.0 + I * z);
}
//
complex<double> UHPtoDisk(complex<double> u)
{
  // u = 0, 1, infty  maps to -1 -i , 1, i  
  return (u - I)/(1.0 - I*u); 
}

//- Rotate z about z0 by 2*k*pi/q 
complex<double> newVertex(complex<double> z,complex<double> z0,int k, int q) {

  complex<double> w( 0.0, 2.0 * sinl(k * M_PI/q) );
  complex<double> a( cosl(k*M_PI/q)*(1.0 - norm(z0)), sinl(k*M_PI/q)*(1.0 + norm(z0)) ); 
  w = w*z0;
  
  //cout<<"New z = "<<-(a*z - w)/(conj(w)*z - conj(a))<<endl;
  return - (a*z - w)/(conj(w)*z - conj(a)); 
}


complex<double> inversion(complex<double> z0, double r)
{
  // z_image conj(z0) = r^2
  return r*2/conj(z0);
}

complex<double> squareInversion(complex<double>z0,double r1,double r2 )
{
  return inversion(inversion(z0, r1),r2);
}

double greens2D(complex<double> z, complex<double> w)
{
  // cout <<" z="<< z.real()<<","<< z.imag()<<" "<<" w="<< w.real()<<","<< w.imag() <<endl;
  return (-1.0/M_PI)*log(abs(z - w)/abs(1.0 - z * conj(w)));
}

/* (3,p)  packing lenght side = 2 ch^{-1} (1/2\sin (\pi/p)) 

cosh(side/2) = 1/(2 sin(\pi/p)) 

r_0=[ sqrt{ 4 cosh^2(side/2) -1} -(sqrt 3)/2]/ sqrt{cosh^2(side/2) -1}

*/


#endif
