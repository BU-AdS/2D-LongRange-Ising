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

//Data file for lattice/analytical propagator data,
void dataDump(vector<Vertex> NodeList, double *phi, Param p, int level, int t_range, int shift) {

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
