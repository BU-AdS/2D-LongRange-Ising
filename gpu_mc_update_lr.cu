#ifdef USE_GPU
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <stdio.h>

#include <gpu_mc_update_lr.cuh>
#include <mc_util.h>

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

void MonteCarlo2DIsing::GPU_wolffClusterAddLR(Param p, int cSpin) {
  
  
  
  
}

void MonteCarlo2DIsing::GPU_wolffUpdateLR(Param p, int iter) {

  int N = p.surfaceVol;

  // CUDA's random number library uses curandState_t to keep track of the seed value
  // we will store a random state for every thread  
  //curandState_t* states;
  
  // allocate space on the GPU for the random states 
  cudaMalloc((void**) &states, N * sizeof(curandState_t));
  
  // invoke the GPU to initialize all of the random states 
  init<<<N, 1>>>(1234, states);

  // allocate an array of unsigned ints on the CPU and GPU 
  double cpu_nums[N];

  cudaMalloc((void**) &gpu_rands, N * sizeof(double));
  
  // invoke the kernel to get some random numbers 
  randoms<<<N, 1>>>(states, gpu_rands);
  
  /*
  // copy the random numbers back 
  cudaMemcpy(cpu_nums, gpu_nums, N * sizeof(double), cudaMemcpyDeviceToHost);

  printf("\n\n\n");
  
  // print them out 
  for (int i = 0; i < N; i++) {
    printf("%f ", cpu_nums[i]);
  }

  printf("\n\n\n");
  */

  // free the memory we allocated for the states and numbers 
  //cudaFree(states);
  //cudaFree(gpu_rands);
  
  //Choose a random spin.
  int i = 1;
  int cSpin = s[i];
  
  // The site belongs to the cluster, so flip it.
  s[i] *= -1;
  phi[i] *= -1;
  
  //This function is recursive and will call itself
  //until all attempts to increase the cluster size
  //have failed.
  GPU_wolffClusterAddLR(p, cSpin);
  
}


#endif
