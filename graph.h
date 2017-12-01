#ifndef ALG_H
#define ALG_H
#include <complex>

using namespace std;

#define Nlinks 7

#define I complex<double>(0.0,1.0)

class Param{
 public:  
  bool bc = false;       //if true, use Dirichlet. If false, use Neumann
  bool Vcentre = true;  //if true, place vertex at centre. If false, use circumcentre.
  bool verbose = false;  //if true, print all data. If false, print summary.
  int MaxIter = 100000;
  double tol = pow(10,-6);
  int q = Nlinks;
  int t = 1;
  double msqr = 0.1;
  int Levels = 3;
  int src_pos = 0;
  char fname[256];
  
  void print(){
    cout<<"Parameter status:"<<endl;
    cout<<"B.C. = "<< (bc ? ("Dirichlet") : ("Neumann") ) << endl;
    cout<<"Centre = "<< (Vcentre ? ("Vertex") : ("Circum") ) << endl;
    cout<<"Levels = "<<Levels<<endl;
    cout<<"MaxIter = "<<MaxIter<<endl;
    cout<<"Tol = "<<tol<<endl;
    cout<<"Triangulation = "<<q<<endl;
    cout<<"TimeSlices = "<<t<<endl;
    cout<<"mass squared = "<<msqr<<endl;
  }

  Param(){
    sprintf(fname, "");
  }
};


/*
struct Param{
  bool bc = true;      //if true, use Dirichlet. If false, use Neumann
  bool Vcentre = false;//if true, place vertex at centre. If false, use circumcentre.
  int Levels = 3;      
  int MaxIter = 100;
  double tol = pow(10,-6);
  int q = Nlinks;
  int t = 5;
  double msqr = 0.1;

  void print(){
    cout<<"Parameter status:"<<endl;
    cout<<"B.C. = "<< (bc ? ("Dirichlet") : ("Neumann") ) << endl;
    cout<<"Centre = "<< (Vcentre ? ("Vertex") : ("Circum") ) << endl;
    cout<<"Levels = "<<Levels<<endl;
    cout<<"MaxIter = "<<MaxIter<<endl;
    cout<<"Tol = "<<tol<<endl;
    cout<<"Triangulation = "<<q<<endl;
    cout<<"TimeSlices = "<<t<<endl;
    cout<<"mass squared = "<<msqr<<endl;
  }
};
*/

struct Vertex{
  int Number;
  int nn[Nlinks+2];
  int fwdLinks;
  complex<double> z;
};


typedef vector<Vertex> Graph;

complex<double> newVertex(complex<double> z,complex<double> z0,int k, int q);
void CheckEdgeLength(const vector<Vertex> NodeList, Param P);
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

//- Edge length from center z = 0
inline double edgeLength(int q) {
  return sqrt( 1 - 4*sinl(M_PI/q)*sinl(M_PI/q) );
}

//Using the formula c(n) = (q-4)*c(n-1) - c(n-2) where c is the
//number of nodes on circumference at level n, we can construct
//the address of the end node on a given level for triangulation q:
//EN(lev,q) = SUM c(n) n=0..level
inline long unsigned int endNode(int lev, Param P) { 
  
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
    //NB!!! There are three nodes on level 0, but the END NODE is addressed as 0.
    //Therefore, when correcting the number of nodes on a circumference, we use
    //3, but when giving the endNode count, we use 2.
    c[0] = 3;       //Three nodes on circumference 0
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
    complex<double>init_1(init_mod*cosl(2.0*M_PI/3.0),init_mod*sinl(2.0*M_PI/3.0));
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
	if(P.verbose) cout<<abs(NodeList[n].z)<<endl;
      }
  }
  cout<<endl;
  for(int i=0; i<N; i++) {
    if(P.verbose) cout<<arr[i]<<" ";
    if(arr[i] != 0) distinct++;
  }
  cout<<endl<<"Count = "<<count<<endl;
  cout<<endl<<"Distinct = "<<distinct<<endl;
  */

}

//- Construct the nearest neighbour table
void BuildGraph(Graph &NodeList, Param P){

  int q = P.q;
  int Levels = P.Levels;
  int T = P.t;
  int offset = endNode(Levels,P) + 1;
  
  if(P.Vcentre == true) {

    //Level 0 spatial: trivial
    for(int mu=1; mu<q+1; mu++) {
      NodeList[0].nn[mu-1] = mu;
    }

    // CONVENTION: The link on the same level, going from node n to (n-1)
    //             is the 0th entry in the neighbour table. The next entries
    //             are the links going to the next higher level. The next
    //             entry is same level link going from n to (n+1). The last
    //             entries go to the lower level. This means all links in the 
    //             neighbour table are listed in anti-clockwise order.
    
    //Level 1
    for(long unsigned int n=endNode(0,P)+1; n<endNode(1,P)+1; n++){
      
      //This is the first node, treat it separately.
      n-1 == 0 ? NodeList[n].nn[0] = endNode(1,P) : NodeList[n].nn[0] = n-1;
      
      //Address of first new node on level l+1 from node a
      int x = endNode(1,P)+1 +(n-1)*(q-4); 
      
      //get new nodes on level 2
      for(int i=1; i<q-2; i++) {
	NodeList[n].nn[i] = x+i-1;
	NodeList[x+i-1].nn[q-1] = n;
	NodeList[x+i-1].fwdLinks = q-3;
	//By definition, the first link in this loop has two back links.
	if(i==1) {
	  NodeList[x+i-1].nn[q-2] = n;
	  NodeList[x+i-1].fwdLinks = q-4;
	  n == 1 ? NodeList[x+i-1].nn[q-1] = q : NodeList[x+i-1].nn[q-1] = n-1;
	}
      }
      NodeList[n].nn[q-2] = n%q+1;
      NodeList[n].nn[q-1] = 0;
    }
    //Fix (q-3) link on final node 
    NodeList[q].nn[q-3] = endNode(1,P)+1;
    
    //Level >=2
    for(int l=2; l<Levels+1; l++){
      
      //Loop over all nodes on this level
      //Get first new node on level l+1
      int x = endNode(l,P)+1;
      
      //Loop over all nodes on this level
      for(long unsigned int n=endNode(l-1,P)+1; n<endNode(l,P)+1; n++){      
	
	//Assign links on the same level 
	//Check if first node
	if(n == endNode(l-1,P)+1) {
	  NodeList[n].nn[0] = endNode(l,P);
	} else {
	  NodeList[n].nn[0] = n-1;
	}
	//Check if last node
	if(n == endNode(l,P)) {
	  NodeList[n].nn[NodeList[n].fwdLinks]   = n+1;
	  NodeList[n].nn[NodeList[n].fwdLinks+1] = endNode(l-1,P)+1;
	}
	else NodeList[n].nn[NodeList[n].fwdLinks+1] = n+1;
	
	//Loop over all links on node n to level l+1. If the node
	//has two back links, there are only (q-4) new links
	//to l+1, else there are (q-3).
	//By deiniftion, the first link in this loop has two back links.
	if(l<Levels) {
	  for(int i=1; i<NodeList[n].fwdLinks+1; i++) {
	    NodeList[n].nn[i] = x+i-1;
	    NodeList[x+i-1].nn[q-1] = n;
	    NodeList[x+i-1].fwdLinks = q-3;
	    if(i==1) {
	      NodeList[x+i-1].nn[q-2] = n;
	      NodeList[x+i-1].fwdLinks = q-4;
	      n == endNode(l-1,P)+1 ? NodeList[x+i-1].nn[q-1] = endNode(l-1,P) : NodeList[x+i-1].nn[q-1] = n-1;
	    }
	  }
	}
	x += NodeList[n].fwdLinks-1;
	
	//fix link q-1 on start node
	NodeList[endNode(l-1,P)+1].nn[q-1]=endNode(l-1,P);
	//fix link q-2 on start node
	NodeList[endNode(l-1,P)+1].nn[q-2]=endNode(l-2,P)+1;
	//fix link q-3 on end node
	if(n == endNode(Levels,P)) NodeList[endNode(l,P)].nn[q-3] = -1;
	else NodeList[endNode(l,P)].nn[q-3] = endNode(l,P) + 1;
      }
    }
    
    //Populate temporal links on t=0 disk
    for(int n=0; n<endNode(Levels,P)+1; n++) {
      //Fwd link
      NodeList[n].nn[q  ] = n + offset;
      //Bkd link
      NodeList[n].nn[q+1] = (T-1)*offset + n;
    }
    
    //Construct disks and t links for 0 < t < T
    for(int t=1; t<T; t++)
      for(int n=0; n<endNode(Levels,P)+1; n++) {
	for(int i=0; i<q; i++) {
	  NodeList[(t-1)*offset + n].nn[i] == -1 ?
	    NodeList[t*offset + n].nn[i] = -1 : 
	    NodeList[t*offset + n].nn[i] =
	    NodeList[(t-1)*offset + n].nn[i] + offset;
	}
	//Fwd link
	NodeList[t*offset + n].nn[q  ] = (t+1)*offset + n;    
	//Bkd link
	NodeList[t*offset + n].nn[q+1] = (t-1)*offset + n;      
      }
    
    //Correct forward t links for t = T-1
    int t=T-1;
    for(int n=0; n<endNode(Levels,P)+1; n++) {
      NodeList[t*offset + n].nn[q] = n;
    }
  }
  else {  
    //Level 0
    for(long unsigned int n=0; n<3; n++){
      
      //This is the first node, treat it separately.
      n == 0 ? NodeList[n].nn[0] = 2 : NodeList[n].nn[0] =  n-1;
      
      //Address of first new node on level 1 from node n
      int x = 3 + n*(q-3);
      
      //get new nodes on level 1
      for(int i=1; i<q-1; i++) {
	NodeList[n].nn[i] = x+i-1;
	NodeList[n].fwdLinks = q-2;
	
	NodeList[x+i-1].nn[q-1] = n;
	NodeList[x+i-1].fwdLinks = q-3;
	//Corrections
	if(i==1) {
	  //This node (3,7,11) has two back links to level 0:
	  //3:  0,2
	  //7:  1,0
	  //11: 2,1
	  NodeList[x+i-1].fwdLinks = q-4;
	  NodeList[x+i-1].nn[q-2] = (x+i-1)%3;
	  NodeList[x+i-1].nn[q-1] = ((x+i-1)%3 + 2)%3;	
	}
	n == 2 ? NodeList[n].nn[q-1] = 0 : NodeList[n].nn[q-1] =  n+1;
      }
    }
    //Fix (q-2) link on final node 
    NodeList[2].nn[q-2] = 3;
    
    //Level 1
    
    //Get first new node on level 2
    int x = endNode(1,P)+1;
    
    //Loop over all nodes on level 1.
    for(long unsigned int n=endNode(0,P)+1; n<endNode(1,P)+1; n++){      
      
      //Assign links on the same level 
      //Check if first node
      if(n == endNode(0,P)+1) {
	NodeList[n].nn[0] = endNode(1,P);
      } else {
	NodeList[n].nn[0] = n-1;
      } //OK
      //Check if last node
      if(n == endNode(1,P)) {
	NodeList[n].nn[NodeList[n].fwdLinks]   = n+1;
	NodeList[n].nn[NodeList[n].fwdLinks+1] = endNode(0,P)+1;
      }
      else NodeList[n].nn[NodeList[n].fwdLinks+1] = n+1;
      
      //Loop over new links
      for(int i=1; i<NodeList[n].fwdLinks+1; i++) {
	NodeList[n].nn[i] = x+i-1;
	NodeList[x+i-1].nn[q-1] = n;
	NodeList[x+i-1].fwdLinks = q-3;
	if(i==1) {
	  NodeList[x+i-1].nn[q-2] = n;
	  NodeList[x+i-1].fwdLinks = q-4;
	  n == endNode(0,P)+1 ? NodeList[x+i-1].nn[q-1] = endNode(1,P) : NodeList[x+i-1].nn[q-1] = n-1;
	}
      }
      x += NodeList[n].fwdLinks-1;
    }
    //Fix (q-3) link on final node 
    NodeList[endNode(1,P)].nn[q-3] = endNode(1,P)+1;
    
    //Level >=2
    for(int l=2; l<Levels+1; l++){
      
      //Get first new node on level l+1
      int x = endNode(l,P)+1;    
      //Loop over all nodes on this level
      for(long unsigned int n=endNode(l-1,P)+1; n<endNode(l,P)+1; n++){      
	
	//Assign links on the same level 
	//Check if first node
	if(n == endNode(l-1,P)+1) {
	  NodeList[n].nn[0] = endNode(l,P);
	} else {
	  NodeList[n].nn[0] = n-1;
	}
	//Check if last node
	if(n == endNode(l,P)) {
	  NodeList[n].nn[NodeList[n].fwdLinks]   = n+1;
	  NodeList[n].nn[NodeList[n].fwdLinks+1] = endNode(l-1,P)+1;
	}
	else NodeList[n].nn[NodeList[n].fwdLinks+1] = n+1;
	
	//Loop over all links on node n to level l+1. If the node
	//has two back links, there are only (q-4) new links
	//to l+1, else there are (q-3).
	//By deiniftion, the first link in this loop has two back links.
	if(l<Levels) {
	  for(int i=1; i<NodeList[n].fwdLinks+1; i++) {
	    NodeList[n].nn[i] = x+i-1;
	    NodeList[x+i-1].nn[q-1] = n;
	    NodeList[x+i-1].fwdLinks = q-3;
	    if(i==1) {
	      NodeList[x+i-1].nn[q-2] = n;
	      NodeList[x+i-1].fwdLinks = q-4;
	      n == endNode(l-1,P)+1 ? NodeList[x+i-1].nn[q-1] = endNode(l-1,P) : NodeList[x+i-1].nn[q-1] = n-1;
	    }
	  }
	}
	x += NodeList[n].fwdLinks-1;
	
	//fix link q-1 on start node
	NodeList[endNode(l-1,P)+1].nn[q-1]=endNode(l-1,P);
	//fix link q-2 on start node
	NodeList[endNode(l-1,P)+1].nn[q-2]=endNode(l-2,P)+1;
	//fix link q-3 on end node
	if(n == endNode(Levels,P)) NodeList[endNode(l,P)].nn[q-3] = -1;
	else NodeList[endNode(l,P)].nn[q-3] = endNode(l,P) + 1;
      }
    }
    
    //Populate temporal links on t=0 disk
    for(int n=0; n<endNode(Levels,P)+1; n++) {
      //Fwd link
      NodeList[n].nn[q  ] = n + offset;
      //Bkd link
      NodeList[n].nn[q+1] = (T-1)*offset + n;
    }
    
    //Construct disks and t links for 0 < t < T
    for(int t=1; t<T; t++)
      for(int n=0; n<endNode(Levels,P)+1; n++) {
	for(int i=0; i<q; i++) {
	  NodeList[(t-1)*offset + n].nn[i] == -1 ?
	    NodeList[t*offset + n].nn[i] = -1 : 
	    NodeList[t*offset + n].nn[i] =
	    NodeList[(t-1)*offset + n].nn[i] + offset;
	}
	//Fwd link
	NodeList[t*offset + n].nn[q  ] = (t+1)*offset + n;    
	//Bkd link
	NodeList[t*offset + n].nn[q+1] = (t-1)*offset + n;
      }
    
    //Correct forward t links for t = T-1
    int t=T-1;
    for(int n=0; n<endNode(Levels,P)+1; n++) {
      NodeList[t*offset + n].nn[q] = n;
    }
  }
    
}

int Mphi(vector<double> &phi, const vector<double> phi0,
	 vector<Vertex> NodeList, Param P) {

  int Levels = P.Levels;
  int q = P.q;
  int T = P.t;
  double msqr = P.msqr;
  bool bc = P.bc;
  int InternalNodes = endNode(Levels-1,P)+1;
  int TotNumber = endNode(Levels,P)+1;
  int offset = TotNumber;
  
  for(int t = 0; t<T; t++) {

    //loop over Interior nodes on all disks    
    for(int i = t*offset; i < t*offset + InternalNodes; i++) {
      //cout<<"interior i="<<i<<" t="<<t<<endl;
      //mass term
      phi[i] = msqr * phi0[i];    
      
      for(int mu = 0; mu < q+2; mu++) {
	phi[i] += (phi0[i] - phi0[NodeList[i].nn[mu]]);
      }
    }
  
    // Dirichlet or Neuman at Exterior Nodes.  
    for(int i = t*offset + InternalNodes; i < t*offset + TotNumber; i++){
      
      //cout<<"Exterior i="<<i<<" t="<<t<<endl;
      //mass term
      phi[i] = msqr * phi0[i];
      if(bc == true) phi[i] += q * phi0[i];
      
      //the zeroth link is always on the same level
      if(bc == true) phi[i] +=  - phi0[NodeList[i].nn[0]];
      else           phi[i] += phi0[i]- phi0[NodeList[i].nn[0]];
      
      //The q-1 (and q-2) link(s) go back one level. 
      //The q-3 (or q-2) link is on the same level.
      //We use the member data fwdLinks to sum only links on
      //the same level, or gong back one level. 
      for(int mu = q-1; mu > (NodeList[i].fwdLinks); mu--) {
	if(bc == true) phi[i] +=  - phi0[NodeList[i].nn[mu]];
	else           phi[i] += phi0[i] - phi0[NodeList[i].nn[mu]];
      }

      //Temporal links at exterior
      if(T>1) {
	for(int mu = q; mu < q+2; mu++) {
	  phi[i] += (phi0[i] - phi0[NodeList[i].nn[mu]]);
	}
      }    
    }
  }
  return 0;  
}


double Minv_phi(vector<double> &phi, const vector<double> phi0,
		const vector<double> b, const vector<Vertex> NodeList,
		Param P)
{
  // CG solutions to Mphi = b 
  //  see http://en.wikipedia.org/wiki/Conjugate_gradient_method
  int Levels = P.Levels;
  int diskN = endNode(Levels,P) + 1;
  int N = P.t*diskN;
  
  vector<double> res(N,0.0), resNew(N,0.0),  pvec(N,0.0), Mpvec(N,0.0), pvec_tmp(N,0.0);
  double alpha, beta, denom;
  double rsq = 0, rsqNew = 0, bsqrt = 0, truersq = 0.0;
  int  i;
  
  for(i = 0; i<N; i++){
    res[i] = b[i];
    pvec[i] = res[i];
    bsqrt += b[i]*b[i];
  }
  bsqrt = sqrt(bsqrt);
  
  int maxIter = P.MaxIter;
  double resEpsilon = 0.000001;
  // iterate till convergence
  rsqNew = 100.0;
  int k = 0;
  
  while( (k<maxIter)&&(sqrt(rsqNew) > resEpsilon*bsqrt) ){
    
    k++;
    rsq = 0;
    for (int i = 0; i < N; i++) rsq += res[i]*res[i];
    
    cout << endl << "In CG at iteration = "<<k <<" Residue Squared  = " << rsq << endl;

    //Mat-Vec operation
    Mphi(Mpvec, pvec, NodeList, P);  

    denom = 0;
    for(i=0; i< N; i++) denom += pvec[i]*Mpvec[i];
    alpha = rsq/denom;
    
    for(i=0; i < N; i++) phi[i] +=  alpha * pvec[i];
    for(i=0; i < N; i++) resNew[i] = res[i]- alpha*Mpvec[i];
    
    // Exit if new residual is small enough
    rsqNew = 0;
    for (i = 0; i < N; i++) rsqNew += resNew[i]*resNew[i];
    
    // Update vec using new residual
    beta = rsqNew / rsq;
    for (i = 0; i < N; i++) {
      pvec[i] = resNew[i] + beta * pvec[i];
      res[i] = resNew[i];
    }
  }
  
  if(k == maxIter) {
    printf("CG: Failed to converge iter = %d, rsq = %e\n", k, rsq); 
    //  Failed convergence 
  }
  
  Mphi(Mpvec, phi, NodeList, P);  
  for(int i=0; i < N ; i++) truersq += (Mpvec[i] - b[i])*(Mpvec[i] - b[i]);
  
  printf("# CG: Converged iter = %d, rsq = %e, truesq = %e\n",k,rsq,truersq);

  return truersq; // Convergence 
}

//- For each node n, with a link to another node,
//  it checks that the neighbour table on the linked
//  node contains the original node n as a neighbour.
void BooleanCheck(Graph &NodeList, Graph &AuxNodeList, Param P){

  int q = P.q;
  int Levels = P.Levels;
  int T = P.t;
  
  for(int n=0; n< T*(endNode(Levels,P)+1); n++) {
    for(int m=0; m<q+2; m++) {
      //Check that the link is valid
      if(NodeList[n].nn[m] != -1) {
	for(int p=0; p<q+2; p++) {
	  //cout<<"Node: "<<n<<"Link: "<<m<<"     Node:"<<NodeList[n].nn[m]<<" Link:"<<p<<endl;
	  //Loop over all links on the linked node,
	  //check if original node exists. 
	  if( n == NodeList[ NodeList[n].nn[m] ].nn[p] ) {
	    AuxNodeList[n].nn[m] = 1;
	  }
	}
      }
    }
  }
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
      for(int n = 0; n < endNode(0,P)+1; n++) {
	cout << endl<< " Node number = " << n + offset << " FL="<<NodeList[n].fwdLinks<<" : ";
	for(int i = 0; i < q+2; i++) cout << NodeList[offset + n].nn[i] << "  ";
      } 
    }
    for(int lev = 1; lev < Levels+1; lev++)  {
      cout << endl << "lev = " << lev << "  T = " << t << endl;
      for(int n = endNode(lev-1,P)+1; n < endNode(lev,P)+1; n++) {
	cout << endl<< " Node number = " << n + offset << " FL="<<NodeList[n].fwdLinks<<" : ";
	for(int i = 0; i < q+2; i++) cout << NodeList[offset + n].nn[i] << "  ";
      }
    }      
  }  
  cout<<endl;
}

void PrintComplexPositions(const vector<Vertex> NodeList, Param P) {

  int Levels = P.Levels;
  
  if(P.verbose) cout<<endl<<"#Printing for Level 0"<<endl;
  for(int n=0; n<endNode(0,P)+1; n++) {
    if(P.verbose) {
      cout<<"n="<<n<<" z="<<NodeList[n].z.real()<<","<<NodeList[n].z.imag();
      cout<<" |z|="<<abs(NodeList[n].z)<<" phi="<<arg(NodeList[n].z);
    }
  }
  for(int l=1; l<Levels+1; l++) {
    if(P.verbose) cout<<endl<<"Printing for Level "<<l<<endl;
    for(int n=endNode(l-1,P)+1; n<endNode(l,P)+1; n++) {
      if(P.verbose) {
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

  if(P.verbose) cout<<endl<<"Checking boundary areas"<<endl;
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
      if(P.verbose) cout<<"n="<<n<<" area "<<k+1<<" = "<<areaGeneral(P, length_01, length_02, length_12)<<endl;
      sig1 += pow(equi_area - areaGeneral(P, length_01, length_02, length_12),2);
      sig2 += pow(ave       - areaGeneral(P, length_01, length_02, length_12),2);
    }
  }

  sig1 /= count - 1;
  sig2 /= count - 1;
  sig1 = sqrt(sig1);
  sig2 = sqrt(sig2);

  cout<<"AREA EQUI = "<<equi_area<<endl;
  cout<<"AREA STD DEV EQUI = "<<sig1<<endl;

  cout<<"AREA AVE = "<<ave<<endl;  
  cout<<"AREA STD DEV AVE = "<<sig2<<endl;  

}

void CheckEdgeLength(const vector<Vertex> NodeList, Param P) {
  
  int q = P.q;
  int Levels = P.Levels;
  double length = 0.0;
  double sig = 0.0;
  int  nn_node;
  bool Vcentre = P.Vcentre;
  double length_0 = d12(NodeList[0].z, NodeList[1].z);
  double tol = 1e-9;

  //Level 0 is specific to how the graph is centred.
  if(Vcentre) {
    if(P.verbose) cout<<" lev =  " << 0 << endl;
    if(P.verbose) cout<<endl<<" Node number = "<<0<<" : "<<endl;
    for(int i = 0; i < q; i++){
      nn_node = NodeList[0].nn[i];
      length = d12(NodeList[0].z, NodeList[nn_node].z);
      if(P.verbose) cout<<" "<<NodeList[0].nn[i]<<" > "<<length<<"  ";
    }
  }
  else {
    if(P.verbose) cout<<" lev = "<<0<<endl;
    if(P.verbose) cout<<endl<<" Node number = "<<0<<" : "<<endl;
    for(int i = 0; i < 2; i++){
      nn_node = NodeList[0].nn[i];
      length = d12(NodeList[0].z, NodeList[nn_node].z);
      if(P.verbose) cout << NodeList[0].nn[i] << " >  " << length<< "  ";
    }
  }
  
  for(int lev = 1; lev < Levels+1; lev++)  {
    if(P.verbose) cout<<endl<<endl<<" lev = "<<lev<<endl;      
    for(int n = endNode(lev-1,P) + 1;n < endNode(lev,P) + 1 ;n++) {
      if(P.verbose) cout<<endl<<" Node number = "<<n<<":"<<endl;
      sig += pow( length_0 - d12(NodeList[n].z, NodeList[NodeList[n].nn[q-1]].z), 2);
      
      for(int i = 0; i <q; i++){
	nn_node = NodeList[n].nn[i];
	if(NodeList[n].nn[i] != -1 ) {
	  length = d12(NodeList[n].z, NodeList[nn_node].z);
	  if(P.verbose) {
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
  cout<<endl;
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
