#ifdef USE_GPU
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <stdio.h>

#include <gpu_mc_update_lr.cuh>
#include <mc_util.h>

bool Check;

int gpu_wc_ave = 0;
int gpu_wc_size = 0;
int gpu_wc_calls = 0;
int gpu_wc_poss = 0;
int gpu_wc_t_size = 0;
int gpu_wc_s_size = 0;

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
			    int x_len, int t1, int x1, int cSpin, bool *gpu_added) {
  
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

    double coupling = 1.0/pow((dx*dx + dt*dt), 12);
    
    if(gpu_rands[idx] < 1 - exp(2*phi_lc*gpu_phi[idx]*coupling)) {
      /*
	printf("Hello from block %d, thread %d : %f %f %f %f %f\n",
	blockIdx.x, threadIdx.x, gpu_rands[idx], phi_lc, gpu_phi[idx],
	coupling, 1 - exp(2*phi_lc*gpu_phi[idx]*coupling));
      */
      
      gpu_added[idx] = true;
      gpu_s[idx] *= -1;
      gpu_phi[idx] *= -1;
    }
  }
  __syncthreads();
}

void MonteCarlo2DIsing::GPU_initRand(Param p, int seed, curandState_t* states) {

  // invoke the GPU to initialize all of the random states 
  init<<<p.surfaceVol, 1>>>(seed, states);
  
}


void MonteCarlo2DIsing::GPU_wolffClusterAddLR(int i, Param p, int cSpin) {
  
  double phi_lc = -1.0*phi[i];
  int S1 = p.S1;
  int Lt = p.Lt;
  int x_len = S1/2 + 1;
  int t_len = Lt/2 + 1;
  int arr_len = x_len*t_len;
  int vol = p.surfaceVol;

  //This implementation parallelises over the entire surface of the lattice.
  //It performs a boolean check to see if the candidate site has the
  //correct spin, then performs the probablisitic test to add the site.

  int t1,x1;

  t1 = i / S1;
  x1 = i % S1;

  //Invoke the kernel to get some random numbers 
  randoms<<<vol, 1>>>(states, gpu_rands);
  double cpu_rands[vol];
  cudaMemcpy(cpu_rands, gpu_rands, vol*sizeof(double), cudaMemcpyDeviceToHost);
  for (int i = 0; i < vol; i++) {
    //printf("%i %f\n", i, cpu_rands[i]);
  }

  cudaMemcpy(gpu_added, added, vol*sizeof(bool), cudaMemcpyHostToDevice);  

  //Test bonds
  cluster_add<<<vol, 1>>>(gpu_rands, gpu_phi, gpu_s, gpu_LR_couplings,
			  phi_lc, S1, Lt, x_len, t1, x1, cSpin, gpu_added);

  cudaMemcpy(added, gpu_added, p.surfaceVol*sizeof(bool), cudaMemcpyDeviceToHost);
}    


void MonteCarlo2DIsing::GPU_wolffUpdateLR(Param p, int i) {
  
  gpu_wc_calls++;
  gpu_wc_s_size = 0;
  
  //Pre-chosen random spin.
  int cSpin = s[i];
  Check = true;
  //The site belongs to the cluster, so flip it.
  s[i] *= -1;
  //phi[i] *= -1;
  
  for(int j=0; j<p.surfaceVol; j++) added[j] = false;
  added[i] = true;

  cudaMemcpy(gpu_phi, phi, p.surfaceVol*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_s, s, p.surfaceVol*sizeof(int), cudaMemcpyHostToDevice);
  
  bool internalCheck = 0;
  while(Check) {
    internalCheck = false;
    for(int j=0; j<p.surfaceVol; j++) {      
      if(added[j] == true ) {
	internalCheck = true;
	added[j] = false;
	gpu_wc_size++;
	GPU_wolffClusterAddLR(j, p, cSpin);
      }
    }
    if(internalCheck == false) Check = false;
    //printf("Size = %d Check = %d ", gpu_wc_s_size, Check);    
  }
  
  cudaMemcpy(phi, gpu_phi, p.surfaceVol*sizeof(double), cudaMemcpyDeviceToHost);  
  cudaMemcpy(s, gpu_s, p.surfaceVol*sizeof(int), cudaMemcpyDeviceToHost);

  phi[i] *= -1;
  
  //printf("\nAve. cluster size = %f (%d/%d)\n", 1.0*gpu_wc_size/gpu_wc_calls,gpu_wc_size,gpu_wc_calls );
}

#endif
