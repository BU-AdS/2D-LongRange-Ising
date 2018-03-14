#ifdef USE_GPU
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <stdio.h>

#include "util.h"

//Used to initialize the random states
__global__ void init(unsigned int seed, curandState_t* states);


//Takes an array of states, ints, and puts a random int into each 
__global__ void randoms(curandState_t* states, double* numbers);


__global__ void cluster_add(double *gpu_rands, double *gpu_phi, int *gpu_s,
			    double *gpu_LR_couplings, int S1, int Lt, int arr_len,
			    int t1, int x1, int cSpin, int *added);

#endif
