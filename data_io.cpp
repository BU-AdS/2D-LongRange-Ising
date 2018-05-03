#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>

#include "util.h"
#include "hyp_util.h"
#include "data_proc.h"

using namespace std;

void writeObservables(double **ind_corr, double *run_corr, int *norm_corr,
		      double **ind_ft_corr, double *run_ft_corr,
		      int idx, observables obs, Param p){

  double norm = 1.0/idx;
  int S1 = p.S1;
  int Lt = p.Lt;
  int vol = S1*Lt;
  int x_len = S1/2 + 1;
  char fname[256];

  //The observables (spec heat, suscep, etc) are used only as a guide
  //to discern criticality, The real quantities of interest are
  //The critical exponents, especally \eta, and its dependence
  //on \sigma, hence MUST have proper error estimates of the
  //correlation function values.

  for(int dth = 0; dth <S1/2+1; dth++) {
    double *jk_err = (double*)malloc((Lt/2+1)*sizeof(double));
    jackknife(ind_corr, run_corr, jk_err, p.n_jkblock, idx, Lt/2+1, dth, p);
    
    sprintf(fname, "correlators_dth%d.dat", dth);
    ofstream file(fname);    
    for(int i=0; i<Lt/2+1; i++) {
      file<<i<<" "<<run_corr[dth + i*x_len]*norm/(norm_corr[dth + i*x_len]);
      file<<" "<<jk_err[i]*norm/(norm_corr[dth + i*x_len])<<endl;      
    }
    file.close();
    free(jk_err);
  }

  for(int dt = 0; dt <Lt/2+1; dt++) {
    double *jk_err = (double*)malloc((S1/2+1)*sizeof(double));
    jackknife(ind_corr, run_corr, jk_err, p.n_jkblock, idx, S1/2+1, dt, p);
    
    sprintf(fname, "correlators_dt%d.dat", dt);
    ofstream file(fname);
    for(int i=0; i<S1/2+1; i++) {
      file<<i<<" "<<run_corr[i + dt*x_len]*norm/(norm_corr[i + dt*x_len]);
      file<<" "<<jk_err[i]*norm/(norm_corr[i + dt*x_len])<<endl;
    }
    file.close();
    free(jk_err);
  }

  //Go to l=2
  for(int l=0; l<3; l++) {
    
    sprintf(fname, "correlators_FTl%d.dat", l);
    ofstream file(fname);    
    double *jk_err = (double*)malloc((Lt/2+1)*sizeof(double));
    jackknifeFT(ind_ft_corr, run_ft_corr, jk_err, p.n_jkblock, idx, l, p);
    
    for(int dt=0; dt<Lt/2+1; dt++) {
      file<<dt<<" "<<run_ft_corr[3*dt + l]*norm<<" "<<jk_err[dt]<<endl;
    }
    file.close();
    free(jk_err);
  }
  
  
  double J = p.J;
  if(p.theory_type != ISING) J = 1.0;
  
  //Jacknife the second moments of means.
  double jkErrSuscep;
  jkErrSuscep = (vol*J)*jackknifeVar(obs.PhiAb_arr, obs.Phi2_arr,
				     obs.Suscep[idx-1]/(vol*J), p.n_jkblock, idx, J);
  
  double jkErrSpecHeat;
  jkErrSpecHeat = J*J*jackknifeVar(obs.E_arr, obs.E2_arr,
				   obs.SpecHeat[idx-1]/(J*J), p.n_jkblock, idx, J*J);
  
  //Jacknife the Binder cumumlant
  double jkErrBinder;
  jkErrBinder = jackknifeBin(obs.Phi2_arr, obs.Phi4_arr,
			     obs.Binder[idx-1], p.n_jkblock, idx);

  //Dump jackknifed observables
  ofstream file_obsJK;
  file_obsJK.open("JKobservables.dat", ios_base::app);
  file_obsJK<<idx<<" "<<p.sigma<<" ";
  if(p.theory_type == ISING) file_obsJK<<p.J<<" ";
  if(p.theory_type == PHI4) file_obsJK<<p.musqr<<" "<<p.lambda<<" ";
  file_obsJK<<obs.SpecHeat[idx-1]<<" "<<jkErrSpecHeat<<" ";
  file_obsJK<<obs.Suscep[idx-1]<<" "<<jkErrSuscep<<" ";
  file_obsJK<<obs.Binder[idx-1]<<" "<<jkErrBinder<<endl;  
  file_obsJK.close();
  
  //Dump raw observable data for any further offline analysis
  ofstream file_obs;
  file_obs.open("observables.dat", ios_base::app);
  int meas = 0;
  for(int i=0; i<10; i++) {
    meas = idx - 10 + i;
    file_obs<<meas<<" "<<p.sigma<<" "<<p.J<<" "<<obs.E_arr[meas]<<" ";
    file_obs<<obs.E2_arr[meas]<<" "<<obs.PhiAb_arr[meas]<<" "<<obs.Phi_arr[meas]<<" ";
    file_obs<<obs.Phi2_arr[meas]<<" "<<obs.Phi4_arr[meas]<<" "<<obs.Suscep[meas]<<" ";
    file_obs<<obs.SpecHeat[meas]<<" "<<obs.Binder[meas]<<endl;
  }
  file_obs.close();  
}

void writePhi3(double **ind_corr_phi_phi3, double *run_corr_phi_phi3,
	       double **ind_corr_phi3_phi3, double *run_corr_phi3_phi3,
	       double **ind_ft_corr_phi_phi3, double *run_ft_corr_phi_phi3,
	       double **ind_ft_corr_phi3_phi3, double *run_ft_corr_phi3_phi3,
	       int *norm_corr, int idx, observables obs, Param p){
  
  double norm = 1.0/idx;
  int S1 = p.S1;
  int Lt = p.Lt;
  int x_len = S1/2 + 1;
  char fname[256];

  for(int dth = 0; dth <S1/2+1; dth++) {

    double *jk_err = (double*)malloc((Lt/2+1)*sizeof(double));
    ofstream file;
    
    //Phi Phi3
    jackknife(ind_corr_phi_phi3, run_corr_phi_phi3, jk_err, p.n_jkblock, idx, Lt/2+1, dth, p);
    sprintf(fname, "correlatorsPhiPhi3_dth%d.dat", dth);    
    file.open(fname);
    for(int i=0; i<Lt/2+1; i++) {
      file<<i<<" "<<run_corr_phi_phi3[dth + i*x_len]*norm/(norm_corr[dth + i*x_len]);
      file<<" "<<jk_err[i]*norm/(norm_corr[dth + i*x_len])<<endl;      
    }
    file.close();
    
    jackknife(ind_corr_phi3_phi3, run_corr_phi3_phi3, jk_err, p.n_jkblock, idx, Lt/2+1, dth, p);    
    sprintf(fname, "correlatorsPhi3Phi3_dth%d.dat", dth);
    file.open(fname);
    for(int i=0; i<Lt/2+1; i++) {
      file<<i<<" "<<run_corr_phi3_phi3[dth + i*x_len]*norm/(norm_corr[dth + i*x_len]);
      file<<" "<<jk_err[i]*norm/(norm_corr[dth + i*x_len])<<endl;      
    }
    file.close();
    free(jk_err);
  }

  for(int dt = 0; dt <Lt/2+1; dt++) {

    double *jk_err = (double*)malloc((S1/2+1)*sizeof(double));
    ofstream file;
    
    jackknife(ind_corr_phi_phi3, run_corr_phi_phi3, jk_err, p.n_jkblock, idx, S1/2+1, dt, p);    
    sprintf(fname, "correlatorsPhiPhi3_dt%d.dat", dt);
    file.open(fname);
    for(int i=0; i<S1/2+1; i++) {
      file<<i<<" "<<run_corr_phi_phi3[i + dt*x_len]*norm/(norm_corr[i + dt*x_len]);
      file<<" "<<jk_err[i]*norm/(norm_corr[i + dt*x_len])<<endl;
    }
    file.close();

    jackknife(ind_corr_phi3_phi3, run_corr_phi3_phi3, jk_err, p.n_jkblock, idx, S1/2+1, dt, p);    
    sprintf(fname, "correlatorsPhi3Phi3_dt%d.dat", dt);
    file.open(fname);
    for(int i=0; i<S1/2+1; i++) {
      file<<i<<" "<<run_corr_phi3_phi3[i + dt*x_len]*norm/(norm_corr[i + dt*x_len]);
      file<<" "<<jk_err[i]*norm/(norm_corr[i + dt*x_len])<<endl;
    }
    file.close();
    free(jk_err);    
  }

  //Go to l=2
  for(int l=0; l<3; l++) {

    double *jk_err = (double*)malloc((Lt/2+1)*sizeof(double));
    ofstream file;
    
    jackknifeFT(ind_ft_corr_phi_phi3, run_ft_corr_phi_phi3, jk_err, p.n_jkblock, idx, l, p);
    sprintf(fname, "correlatorsPhiPhi3_FTl%d.dat", l);
    file.open(fname);
    for(int dt=0; dt<Lt/2+1; dt++) {
      file<<dt<<" "<<run_ft_corr_phi_phi3[3*dt + l]*norm<<" "<<jk_err[dt]<<endl;
    }
    file.close();
    
    jackknifeFT(ind_ft_corr_phi3_phi3, run_ft_corr_phi3_phi3, jk_err, p.n_jkblock, idx, l, p);
    sprintf(fname, "correlatorsPhi3Phi3_FTl%d.dat", l);
    file.open(fname);
    for(int dt=0; dt<Lt/2+1; dt++) {
      file<<dt<<" "<<run_ft_corr_phi3_phi3[3*dt + l]*norm<<" "<<jk_err[dt]<<endl;
    }
    file.close();
    free(jk_err);
  }  
}


//Data file for lattice/analytical propagator data,
void dataDump(vector<Vertex> NodeList, double *phi, Param p, int level,
	      int t_range, int shift) {

  long unsigned int TotNumber = (endNode(p.Levels,p) + 1) * p.Lt;
  double norm = 0.0;
  for(long unsigned int i = 0;i < TotNumber; i++) norm += phi[i]*phi[i];
  for(long unsigned int i = 0;i < TotNumber; i++) phi[i] /= sqrt(norm); 
  
  long unsigned int j = p.src_pos;

  int T = p.Lt;
  int T_offset = 0;
  double theta = 0.0;
  //double msqr = (double)p.msqr + (double)p.delta_msqr*shift;
  double delta = p.Lt < 2 ? 0.5 + sqrt(0.25 + p.msqr) : 1.0 + sqrt(1 + p.msqr);
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
		p.Lt,
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
	  int t1 = j / (TotNumber/p.Lt);
	  int t2 = i / (TotNumber/p.Lt);
	  //Assume PBC.
	  int delta_t = (t2-t1);// > p.Lt/2 ? (t2-t1) - p.Lt : (t2-t1);
	  
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
  int T = P.Lt;
  int t_offset  = 0;
  T == 1 ? t_offset = 0 : t_offset = 2;
  
  for(int t=0; t<T; t++) {
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
    for(int lev = 1; lev < Levels+1; lev++)  {
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

void visualiserIsing(int *s, Param p) {  
  
  for(int i=0; i<p.S1; i++) {
    for(int j=0; j<p.Lt; j++) {
      if(s[i + p.S1*j] < 0) cout<<"\033[1;41m \033[0m";
      else cout<<"\033[1;44m \033[0m";
    }
    cout<<endl;
  }
}

void visualiserIsingCluster(int *s, bool *cluster, Param p) {  
  
  for(int i=0; i<p.S1; i++) {
    for(int j=0; j<p.Lt; j++) {
      if(cluster[i +p.S1*j]) cout<<"\033[1;42m \033[0m";
      else if(s[i + p.S1*j] < 0) cout<<"\033[1;41m \033[0m";
      else cout<<"\033[1;44m \033[0m";
      
    }
    cout<<endl;
  }
}

void visualiserSqr(double *phi_cyl, double barr, Param p) {  
  
  for(int i=0; i<p.S1; i++) {
    for(int j=0; j<p.Lt; j++) {
      if(phi_cyl[i + p.S1*j] < -1.5*barr) cout<<"\033[1;41m \033[0m";
      if(-1.5*barr < phi_cyl[i + p.S1*j] && phi_cyl[i + p.S1*j] < -0.5*barr) cout<<"\033[1;43m \033[0m";
      if(-0.5*barr < phi_cyl[i + p.S1*j] && phi_cyl[i + p.S1*j] <  0.5*barr) cout<<"\033[1;42m \033[0m";
      if( 0.5*barr < phi_cyl[i + p.S1*j] && phi_cyl[i + p.S1*j] <  1.5*barr) cout<<"\033[1;46m \033[0m";
      if( 1.5*barr < phi_cyl[i + p.S1*j]) cout<<"\033[1;44m \033[0m";
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
