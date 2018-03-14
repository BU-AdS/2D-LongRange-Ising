#ifdef USE_GPU
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <vector>

#include <gpu_mc_update_lr.cuh>
#include <mc_util.h>

int Check;

// this GPU kernel function is used to initialize the random states
__global__ void init(unsigned int seed, curandState_t* states) {

  // initialize the state 
  curand_init(seed,       // the seed 
	      blockIdx.x, // use thread id? 
	      0,          // How much extra we advance in the sequence per call
	      &states[blockIdx.x]);
}

// this GPU kernel takes an array of states, ints, and puts a random int into each 
__global__ void randoms(curandState_t* states, double* numbers) {
  
  // curand works like rand - except that it takes a state as a parameter 
  numbers[blockIdx.x] = curand_uniform(&states[blockIdx.x]);
  
}

__global__ void cluster_add(double *gpu_rands, double *gpu_phi, int *gpu_s,
			    double *gpu_LR_couplings, double phi_lc, int S1, int Lt,
			    int arr_len, int t1, int x1, int cSpin, int *added) {
  
  int idx = blockIdx.x;
  if(gpu_s[idx] == cSpin) {
    int t2,x2,dt,dx;
    
    //Index divided by circumference, using the int floor feature/bug,
    //gives the timeslice index.
    t2 = idx / S1;
    dt = abs(t2-t1) > Lt/2 ? Lt - abs(t2-t1) : abs(t2-t1);
    
    //The index modulo the circumference gives the spatial index.
    x2 = idx % S1;            
    dx = abs(x2-x1) > S1/2 ? S1 - abs(x2-x1) : abs(x2-x1);      
    
    if(gpu_rands[blockIdx.x] < 1 - exp(2*phi_lc*gpu_phi[idx]*
				       gpu_LR_couplings[dt+dx*arr_len])) {
      added[blockIdx.x] = 1;
    }       
  }
}

void MonteCarlo2DIsing::GPU_wolffClusterAddLR(int i, Param p, int cSpin,
					      double *gpu_rands, int *added) {
  
  double phi_lc = phi[i];
  int arr_len = p.surfaceVol/2+1;
  int S1 = p.S1;
  int Lt = p.Lt;

  //This implementation parallelises over the entire surface of the lattice.
  //It performs a boolean check to see if the candidate site has the
  //correct spin, then performs the probablisitic test to add the site.

  int t1,x1;

  t1 = i / S1;
  x1 = i % S1;

  //Invoke the kernel to get some random numbers 
  randoms<<<p.surfaceVol, 1>>>(states, gpu_rands);
  //Test bonds
  cluster_add<<<p.surfaceVol, 1>>>(gpu_rands, gpu_phi, gpu_s, gpu_LR_couplings,
				   phi_lc, S1, Lt, arr_len, t1, x1, cSpin, added);
  
  
}

void MonteCarlo2DIsing::GPU_wolffUpdateLR(Param p, int iter) {

  // allocate space on the GPU for the random states 
  cudaMalloc((void**) &states, p.surfaceVol*sizeof(curandState_t));  
  // invoke the GPU to initialize all of the random states 
  init<<<p.surfaceVol, 1>>>(1234, states);
  cudaMalloc((void**) &gpu_rands, p.surfaceVol*sizeof(double));
  
  //Choose a random spin.
  int i = 1;
  int cSpin = s[i];
  int added[p.surfaceVol];
  for(int j=0; j<p.surfaceVol; j++) added[j] = 0;
  added[i] = 1;
  Check = 1;
  
  // The site belongs to the cluster, so flip it.
  s[i] *= -1;
  phi[i] *= -1;
  
  cudaMemcpy(gpu_phi, phi, p.surfaceVol*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_s, s, p.surfaceVol*sizeof(int), cudaMemcpyHostToDevice);
  
  //This function is recursive and will call itself
  //until all attempts to increase the cluster size
  //have failed.
  int site = 0;
  while(Check > 0) {
    int j=0;
    while(site != 0) {
      if(added[j] !=0 ) {
	site = j;
	added[j] = 0;
      }
      j++;
    }
    GPU_wolffClusterAddLR(site, p, cSpin, gpu_rands, added);
    --Check;
  }
}

#endif
