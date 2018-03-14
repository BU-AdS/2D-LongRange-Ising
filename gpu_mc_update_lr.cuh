#ifdef USE_GPU
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <stdio.h>

#include "util.h"

// this GPU kernel function is used to initialize the random states
__global__ void init(unsigned int seed, curandState_t* states);


// this GPU kernel takes an array of states, ints, and puts a random int into each 
__global__ void randoms(curandState_t* states, double* numbers);


//void GPU_wolffClusterAddLR(int i, int *s, int cSpin, double *LR_couplings,
//			   double *gpu_rands, double *phi, Param p);

//void GPU_wolffUpdateLR(double *phi, int *s, Param p,
//		       double *LR_couplings, int iter);


#endif
