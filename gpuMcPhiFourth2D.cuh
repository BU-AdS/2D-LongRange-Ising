#ifdef USE_GPU
#ifndef GPUMCPHIFOURTH2D_H
#define GPUMCPHIFOURTH2D_H
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

#endif
#endif
