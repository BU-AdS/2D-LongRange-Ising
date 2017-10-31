#ifndef ALG_H
#define ALG_H
#include <complex>

using namespace std;

#define I complex<double>(0.0,1.0)

#define Nlinks 7
//#define M_PI 3.1415926535897932385

struct Param{
  bool bc = true;
  int Levels = 3;
  int MaxIter = 1000;
  double tol = pow(10,-6);
  int q = Nlinks;
  int t = 20;
  double msqr = 0.1;
};

struct Vertex{
  int Number;
  int nn[Nlinks+2];
  int fwdLinks;
  complex<double> z;
};

typedef vector<Vertex> Graph;
#define I complex<double>(0.0,1.0)

complex<double> newVertex(complex<double> z,complex<double> z0,int k, int q);

#if 0
void GetComplexPositions(Graph &NodeList, int Levels, int q);
void PrintComplexPositions(const vector<Vertex> NodeList, int Levels, int q);
  
void BuildGraph(Graph &NodeList, int Levels, int q);
void BooleanCheck(Graph &NodeList, Graph &AuxNodeList, int Levels, int q);
void PrintNodeTables(vector<Vertex> NodeList, int  Levels, int q);

//CG routine
int  Mphi(vector<double> &phi, const  vector<double> phi0,
	  vector<Vertex> NodeList, int  Levels, int q , bool bc );
double Minv_phi(vector<double> &phi, const vector<double> phi0,
		const  vector<double>  b,  const vector<Vertex> NodeList,
		int  Levels, int q, bool bc);

#endif 
//Hyperbolic Geometry Routines


void CheckEdgeLength(const vector<Vertex> NodeList, Param P);
complex<double> T(complex<double> z,  complex<double> w);
complex<double> R(complex<double> z, complex<double> omega);
complex<double> flip(complex<double> z, complex<double> z1, complex<double> z2);
double s(complex<double> z);
double r(double s );
double d12(complex<double> z1, complex<double> z2);
double s3p(int q);
double area3q(int q);
double centralRad(double s);
complex<double> DisktoUHP(complex<double> z);
complex<double> UHPtoDisk(complex<double> u);
complex<double> inversion(complex<double> z0, double r);
complex<double> squareInversion(complex<double>z0,double r1,double r2 );
double greens2D(complex<double> z, complex<double> w);



//Using the formula c(n) = (q-4)*c(n-1) - c(n-2) where c is the
//number of nodes on circumference at level n, we can construct
//the address of the end node on a given level for triangulation q:
//EN(lev,q) = SUM c(n) n=0..level
inline long unsigned int endNode(int lev, int q) { 

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

//- Edge length from center z = 0
inline double edgeLength(int q) {
  return sqrt( 1 - 4*sin(M_PI/q)*sin(M_PI/q) );
}


//- Rotate z about z0 by 2*k*pi/q 
complex<double> newVertex(complex<double> z,complex<double> z0,int k, int q) {
  complex<double> w( 0.0, 2.0 * sin(k * M_PI/q) );
  complex<double> a( cos(k*M_PI/q)*(1 - norm(z0)), sin(k*M_PI/q)*(1 + norm(z0)) ); 
  w = w*z0;  
  return - (a*z - w)/(conj(w)*z - conj(a)); 
}

//- Get the z coordinates of every node on the Poincare disk 
void GetComplexPositions(Graph &NodeList, Param P){

  int q = P.q;
  int Levels = P.Levels;
  
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
    for(long unsigned int n=endNode(l-1,q)+1; n<endNode(l,q)+1; n++) {
      for(int k=0; k<q; k++) {
	if(NodeList[n].nn[k] != -1) {
	  NodeList[NodeList[n].nn[k]].z = newVertex(NodeList[NodeList[n].nn[0]].z, NodeList[n].z, k, q);	  
	}	
      }      
    }
  }  
}

//- Construct the nearest neighbour table, and deduce the z
//- coordinates of each node.
void BuildGraph(Graph &NodeList, Param P){

  int q = P.q;
  int Levels = P.Levels;
  int T = P.t;

  int offset = endNode(Levels,q) + 1;
    
  //Level 0 spatial: trivial
  for(int mu=1; mu<q+1; mu++) 
    NodeList[0].nn[mu-1] = mu;
  
  // CONVENTION: The link on the same level, going from node n to (n-1)
  //             is the 0th entry in the neighbour table. The next entries
  //             are the links going to the next higher level. The next
  //             entry is same level link going from n to (n+1). The last
  //             entries go to the lower level. This means all links in the 
  //             neighbour table are listed in anti-clockwise order.
  
  //Level 1
  for(long unsigned int n=endNode(0,q)+1; n<endNode(1,q)+1; n++){
    
    //This is the first node, treat it separately.
    n-1 == 0 ? NodeList[n].nn[0] = endNode(1,q) : NodeList[n].nn[0] = n-1;
    
    //Address of first new node on level l+1 from node a
    int x = endNode(1,q)+1 +(n-1)*(q-4); 
    
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
  NodeList[q].nn[q-3] = endNode(1,q)+1;
  
  //Level >=2
  for(int l=2; l<Levels+1; l++){
    
    //Loop over all nodes on this level
    //Get first new node on level l+1
    int x = endNode(l,q)+1;
    
    //Loop over all nodes on this level
    for(long unsigned int n=endNode(l-1,q)+1; n<endNode(l,q)+1; n++){      
      
      //Assign links on the same level 
      //Check if first node
      if(n == endNode(l-1,q)+1) {
	NodeList[n].nn[0] = endNode(l,q);
      } else {
	NodeList[n].nn[0] = n-1;
      }
      //Check if last node
      if(n == endNode(l,q)) {
	NodeList[n].nn[NodeList[n].fwdLinks]   = n+1;
	NodeList[n].nn[NodeList[n].fwdLinks+1] = endNode(l-1,q)+1;
      }
      else NodeList[n].nn[NodeList[n].fwdLinks+1] = n+1;
      
      //Loop over all links on node n to level l+1. If the node
      //has two back links, there are only (q-4) new links
      //to l+1, else there are (q-3).
      //By definition, the first link in this loop has two back links.
      if(l<Levels) {
	for(int i=1; i<NodeList[n].fwdLinks+1; i++) {
	  NodeList[n].nn[i] = x+i-1;
	  NodeList[x+i-1].nn[q-1] = n;
	  NodeList[x+i-1].fwdLinks = q-3;
	  if(i==1) {
	    NodeList[x+i-1].nn[q-2] = n;
	    NodeList[x+i-1].fwdLinks = q-4;
	    n == endNode(l-1,q)+1 ? NodeList[x+i-1].nn[q-1] = endNode(l-1,q) : NodeList[x+i-1].nn[q-1] = n-1;
	  }
	}
      }
      x += NodeList[n].fwdLinks-1;
      
      //fix link q-1 on start node
      NodeList[endNode(l-1,q)+1].nn[q-1]=endNode(l-1,q);
      //fix link q-2 on start node
      NodeList[endNode(l-1,q)+1].nn[q-2]=endNode(l-2,q)+1;
      //fix link q-3 on end node
      if(n == endNode(Levels,q)) NodeList[endNode(l,q)].nn[q-3] = -1;
      else NodeList[endNode(l,q)].nn[q-3] = endNode(l,q) + 1;
    }
  }

  //Populate temporal links on t=0 disk
  for(int n=0; n<endNode(Levels,q)+1; n++) {
    //Fwd link
    NodeList[n].nn[q  ] = n + offset;
    //Bkd link
    NodeList[n].nn[q+1] = (T-1)*offset + n;
  }
  
  //Construct disks and t links for 0 < t < T
  for(int t=1; t<T; t++)
    for(int n=0; n<endNode(Levels,q)+1; n++) {
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
  for(int n=0; n<endNode(Levels,q)+1; n++) {
    NodeList[t*offset + n].nn[q] = n;
  }
}



int Mphi_t(vector<double> &phi, const vector<double> phi0,
	   vector<Vertex> NodeList, Param p) {

  int Levels = p.Levels;
  int q = p.q;
  int InternalNodes = endNode(Levels-1,q);
  int TotNumber =  endNode(Levels,q)+1;
  int T = p.t;
  int offset = TotNumber;
  
  double msqr = p.msqr;
  bool bc = p.bc;

  for(int t = 0; t<T; t++) {

    //loop over Interior nodes on all disks    
    for(int i = t*offset; i < t*offset + InternalNodes+1; i++) {
      //cout<<"interior i="<<i<<" t="<<t<<endl;
      //mass term
      phi[i] = msqr * phi0[i];    
      
      for(int mu = 0; mu < q+2; mu++) {
	phi[i] += (phi0[i] - phi0[NodeList[i].nn[mu]]);
      }
    }
  
    // Dirichlet or Neuman at Exterior Nodes.  
    for(int i = t*offset + InternalNodes+1; i < t*offset + TotNumber; i++){

      //cout<<"Exterior i="<<i<<" t="<<t<<endl;
      //mass term
      phi[i] = msqr * phi0[i];
      if(bc == true) phi[i] += q * phi0[i];
      
      //the zeroth link is always on the same level
      if(bc == true) phi[i] +=  - phi0[NodeList[i].nn[0]];
      else           phi[i] += phi0[i]- phi0[NodeList[i].nn[0]];
      
      //The q-1 (and q-2) link(s) go back one level. 
      //The q-3 (or q-2) link is on the same level.
      //We use the member data fwdLinks to sum accordingly.
      for(int mu = q-1; mu > (NodeList[i].fwdLinks); mu--) {
	if(bc == true) phi[i] +=  - phi0[NodeList[i].nn[mu]];
	else           phi[i] += phi0[i] - phi0[NodeList[i].nn[mu]];
      }

      //Temporal links at exterior
      for(int mu = q; mu < q+2; mu++) {
	phi[i] += (phi0[i] - phi0[NodeList[i].nn[mu]]);
      }
    }
    
    
  }
  return 0;  
}


double Minv_phi_t(vector<double> &phi, const vector<double> phi0,
		  const vector<double> b, const vector<Vertex> NodeList,
		  Param p)
{
  // CG solutions to Mphi = b 
  //  see http://en.wikipedia.org/wiki/Conjugate_gradient_method
  int q = p.q;
  int Levels = p.Levels;
  int diskN = endNode(Levels,q) + 1;
  int N = p.t*diskN;
  
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
  
  int maxIter = p.MaxIter;
  double resEpsilon = 0.000001;
  // iterate till convergence
  rsqNew = 100.0;
  int k = 0;
  
  while( (k<maxIter)&&(sqrt(rsqNew) > resEpsilon*bsqrt) ){
    
    k++;
    rsq = 0;
    for (int i = 0; i < N; i++) rsq += res[i]*res[i];
    
    cout << endl << "In CG at iteration = "<<k <<" Residue Squared   = " << rsq << endl;

    //Mat-Vec operation
    Mphi_t(Mpvec, pvec, NodeList, p);  

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
    printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq); 
    //  Failed convergence 
  }
  
  Mphi_t(Mpvec, phi, NodeList, p );  
  for(int i=0; i < N ; i++) truersq += (Mpvec[i] - b[i])*(Mpvec[i] - b[i]);
  
  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  return truersq; // Convergence 
}


int Mphi(vector<double> &phi, const vector<double> phi0, vector<Vertex> NodeList, Param p) {
  int Levels = p.Levels;
  int q = p.q;
  int InternalNodes = endNode(Levels-1,q);
  int TotNumber =  endNode(Levels,q)+1;
  
  double msqr = p.msqr;
  bool bc = p.bc;
  
  // Interior nodes.
  for(int i = 0; i < InternalNodes+1; i++) {
    //mass term
    phi[i] = msqr * phi0[i];    
    
    for(int mu = 0; mu < q; mu++) {
      phi[i] += (phi0[i] - phi0[NodeList[i].nn[mu]]);
    }
    }
  
  // Dirichlet or Neuman  at Exterior Nodes.  
  for(int i = InternalNodes+1; i < TotNumber+1; i++){
    
    //mass term
    phi[i] = msqr * phi0[i];
    if(bc == true) phi[i] += q * phi0[i];
    
    //the zeroth link is always on the same level
    if(bc == true) phi[i] +=  - phi0[NodeList[i].nn[0]];
    else           phi[i] += phi0[i]- phi0[NodeList[i].nn[0]];
    
    //The last 2 or 1 link(s) go back one level. 
    //The q-3 (or q-2) link is on the same level.
    //We use the member data fwdLinks to sum accordingly.
    for(int mu = q-1; mu > (NodeList[i].fwdLinks); mu--) {
      if(bc == true) phi[i] +=  - phi0[NodeList[i].nn[mu]];
      else           phi[i] += phi0[i] - phi0[NodeList[i].nn[mu]];
    }
  }
  
  return 0;  
}


double Minv_phi(vector<double> &phi, const vector<double> phi0,
		const vector<double> b, const vector<Vertex> NodeList, Param p)
{
  // CG solutions to Mphi = b 
  //  see http://en.wikipedia.org/wiki/Conjugate_gradient_method
  int q = p.q;
  int Levels = p.Levels;
  int N = endNode(Levels,q) + 1;
  
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
  
  int maxIter = p.MaxIter;
  double resEpsilon = 0.000001;
  // iterate till convergence
  rsqNew = 100.0;
  int k = 0;
  
  while( (k<maxIter)&&(sqrt(rsqNew) > resEpsilon*bsqrt) ){
    
    k++;
    rsq = 0;
    for (int i = 0; i < N; i++) rsq += res[i]*res[i];
    
    cout << endl << "In CG at iteration = "<<k <<" Residue Squared   = " << rsq << endl;

    //Mat-Vec operation
    Mphi(Mpvec, pvec, NodeList, p);  

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
    printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq); 
    //  Failed convergence 
  }
  
  Mphi(Mpvec, phi, NodeList, p );  
  for(int i=0; i < N ; i++) truersq += (Mpvec[i] - b[i])*(Mpvec[i] - b[i]);
  
  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  return truersq; // Convergence 
}



//- For each node n, with a link to another node,
//  it checks that the neighbour table on the linked
//  node contains the original node n as a neighbour.
void BooleanCheck(Graph &NodeList, Graph &AuxNodeList, Param P){

  int q = P.q;
  int Levels = P.Levels;
  int T = P.t;
  
  for(int n=0; n< T*(endNode(Levels,q)+1); n++) {
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
  int offset = 0;
  
  for(int t=0; t<T; t++) {
    
    offset = t*( endNode(Levels,q) + 1 );
    
    cout << "lev = " << 0 << "  T = " << t << endl;
    cout << endl<< " Node number = " << 0 + offset << " : ";
    for(int i = 0; i < q+2; i++){
      cout << NodeList[offset + 0].nn[i] << "  ";
    }
    
    for(int lev = 1; lev < Levels+1; lev++)  {
      int  end = endNode(lev-1,q);
      int  endnew =   endNode(lev,q);
      cout << endl << endl  <<  "lev = " << lev << "  T = " << t << endl;
      for(int n = offset + end + 1; n < offset + endnew +1; n++) {
	cout<< endl << " Node number = " << n << " : ";    
	for(int i = 0; i <q+2; i++){
	  cout << NodeList[n].nn[i] << "  ";
	}
      }
    }
    cout<<endl;
  }
}

void PrintComplexPositions(const vector<Vertex> NodeList, Param P) {

  int q = P.q;
  int Levels = P.Levels;
  
  cout<<"#Printing for origin "<<endl;
  cout<<"n="<<0<<" z="<<NodeList[0].z.real()<<","<<NodeList[0].z.imag()<<" "<<"|z|="<<sqrt(norm(NodeList[0].z))<<endl;
  for(int l=1; l<Levels+1; l++) {
    cout<<"#Printing for level "<<l<<endl;
    for(int n=endNode(l-1,q)+1; n<endNode(l,q)+1; n++) {
      //      cout<<"n="<<n<<" z="<<NodeList[n].z.real()<<","<<NodeList[n].z.imag()<<" "<<"|z|="<<abs(NodeList[n].z)<<" phase  = " << arg(NodeList[n].z)  << endl;
      cout<<  n  << "        "<<NodeList[n].z.real()<<"         "<<NodeList[n].z.imag()<<"          " <<abs(NodeList[n].z)<<"           " << arg(NodeList[n].z)  << endl;
    }
  }
}

void CheckEdgeLength(const vector<Vertex> NodeList, Param P) {

  int q = P.q;
  int Levels = P.Levels;
  double length(10);
  int  nn_node;
  
  cout << " lev =  " << 0 << endl;
  cout << endl<< " Node  number = " << 0 << " : ";
  for(int i = 0; i < q; i++){
    nn_node = NodeList[0].nn[i];
    length = d12(NodeList[0].z, NodeList[nn_node].z);
    cout << NodeList[0].nn[i] << " >  " << length   << "  ";
  }
  
  for(int lev = 1; lev < Levels+1; lev++)  {
    int  end = endNode(lev-1,q);
    int  endnew =   endNode(lev,q);
    cout << endl << endl  <<  " lev =  " << lev << endl;      
    for(int n = end + 1;n < endnew +1;n++) {
      cout<< endl << " Node number = " << n << " : ";    
      for(int i = 0; i <q; i++){
	nn_node = NodeList[n].nn[i];
	if(NodeList[n].nn[i] != -1 ) {
	  length = d12(NodeList[n].z, NodeList[nn_node].z);
	  cout << NodeList[n].nn[i] << " >  " << length << "  ";
	}
      }
    }
  }
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
   return log((1+abs(z))/(1-abs(z)));
}

//Real ot Geodesic distance  s from orgina
double  r(double s )
{
  return tanh(s/2);
}

//Geodesic distandce from z1 to z2
double d12(complex<double> z1, complex<double> z2) {
  double length =  2.0*asinh(abs(z1 - z2)/sqrt((1.0 - norm(z1))*(1.0 - norm(z2))));
  return length;
}

// length of arc q fold triangle to origin.
double s3p(int q)
{  //vertex centered Arc lengeth
  return 2.0*acosh(1.0/sin(M_PI/(double)q));
}

// Area equilateral triangel with angles 2 pi/q
double area3q(int q)
{
  return M_PI - 3.0*(2.0*M_PI/(double)q);
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
  return (u -I)/(1.0 - I * u); 
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
