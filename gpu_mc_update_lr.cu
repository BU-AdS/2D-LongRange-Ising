#ifdef USE_GPU
//#define DEBUG
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <cooperative_groups.h>
#include <stdio.h>

#include <gpu_mc_update_lr.cuh>
#include <mc_util.h>

int gpu_wc_ave = 0;
int gpu_wc_size = 0;
int gpu_wc_calls = 0;
int gpu_wc_poss = 0;
int gpu_wc_t_size = 0;
int gpu_wc_s_size = 0;
int gpu_accept = 0;
int gpu_tries = 0;

const int sze = 512;
extern int CACHE_LINE_SIZE;
int chunk = CACHE_LINE_SIZE/sizeof(double);

using namespace cooperative_groups;


//--------------------------------------------------------------------------
//Simple CUDA error parser
//--------------------------------------------------------------------------
void Check_CUDA_Error(const char *message, int site, int iter) {
  cudaError_t error = cudaGetLastError();
  if(error!=cudaSuccess) {
    fprintf(stderr, "ERROR: %s at site %d iter %d: %s\n", message, site, iter, cudaGetErrorString(error) );
    exit(-1);
  }
}


//--------------------------------------------------------------------------
//Device function to reconstruct double prec value from
//texture cache (two single prec caches).
//--------------------------------------------------------------------------
__device__ __inline__ double fetch_double(int2 i){
  return __hiloint2double(i.y, i.x);
}


//--------------------------------------------------------------------------
//Kernel to convert booleans to a summed integer
//--------------------------------------------------------------------------
__global__ void count_added(bool *added, int *sum){

  __shared__ int sdata[sze];
  int bid = blockIdx.x;
  int tid = threadIdx.x;
  int bdim = blockDim.x;

  //Transfer bool to shared int array
  for(int a=0; a<bdim; a++) {
    if(added[bid*bdim + tid]) sdata[tid]++;
  }

  //Do block reduction in shared memory  
  for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
    if (tid < s) sdata[tid] += sdata[tid + s];
    __syncthreads();
  }
  
  auto g = this_thread_block();
  if (g.thread_rank() == 0) atomicAdd(sum, sdata[0]);

}


//--------------------------------------------------------------------------
//Construct and bind a texture array in double prec
//--------------------------------------------------------------------------
void texAndBind(cudaTextureObject_t &texObject, double *gpu_array, int vol) {
  
  cudaTextureDesc td;
  memset(&td, 0, sizeof(td));
  td.normalizedCoords = 0;
  td.addressMode[0] = cudaAddressModeClamp;
  td.readMode = cudaReadModeElementType;

  struct cudaResourceDesc resDesc;
  memset(&resDesc, 0, sizeof(resDesc));
  resDesc.resType = cudaResourceTypeLinear;
  resDesc.res.linear.devPtr = (void*)gpu_array;
  resDesc.res.linear.sizeInBytes = vol*sizeof(double);
  resDesc.res.linear.desc.f = cudaChannelFormatKindSigned;
  resDesc.res.linear.desc.x = 8*sizeof(int);
  resDesc.res.linear.desc.y = 8*sizeof(int);
  resDesc.res.linear.desc.z = 0;
  resDesc.res.linear.desc.w = 0;

  cudaCreateTextureObject(&texObject, &resDesc, &td, NULL);

  /*
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 1);
    int texel_size = 4 * sizeof(int);
    printf("size of int = %d\n", sizeof(int));
    unsigned long texels = resDesc.res.linear.sizeInBytes / texel_size;
    if (texels > (unsigned)prop.maxTexture1DLinear) {
    printf("Attempting to bind too large a texture %lu > %d\n", texels, prop.maxTexture1DLinear);
    } else {
    printf("Bound a texture of size %lu < %d\n", texels, prop.maxTexture1DLinear);
    }
  */
}



//--------------------------------------------------------------------------
//Copy arrays from device to host
//--------------------------------------------------------------------------
void MonteCarlo2DIsing::GPU_copyArraysToHost(Param p) {
  
  int vol = p.S1*p.Lt;
  cudaMemcpy(phi, gpu_phi, vol*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(s, gpu_s, vol*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(cpu_added, gpu_added, vol*sizeof(bool), cudaMemcpyDeviceToHost);
  
}


//--------------------------------------------------------------------------
//Copy arrays from host to device
//--------------------------------------------------------------------------
void MonteCarlo2DIsing::GPU_copyArraysToDevice(Param p) {
  
  int vol = p.S1*p.Lt;
  cudaMemcpy(gpu_phi, phi, vol*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_s, s, vol*sizeof(int), cudaMemcpyHostToDevice);  
  cudaMemcpy(gpu_added, cpu_added, vol*sizeof(bool), cudaMemcpyHostToDevice);
  
}


//--------------------------------------------------------------------------
//Given the array of changed sites from Metropolis, update the phi array
//--------------------------------------------------------------------------
__global__ void updatePhiAndSpin(double *gpu_rands, double *gpu_phi, int *gpu_s,				 
				 bool *gpu_added) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  
  if(gpu_added[idx]) {
    gpu_phi[idx] = gpu_rands[idx];
    gpu_rands[idx] < 0 ? gpu_s[idx] = -1 : gpu_s[idx] = 1;
  }
}


//--------------------------------------------------------------------------
//Initialise the GPU random number array from the CPU
//--------------------------------------------------------------------------
void MonteCarlo2DIsing::GPU_initRand(Param p, int seed,
				     curandState_t* states) {
  
  // invoke the GPU (from the CPU) to initialize the random states 
  init<<<p.surfaceVol/sze, sze>>>(seed, states);
  
}


//--------------------------------------------------------------------------
//Initialize the random states
//--------------------------------------------------------------------------
__global__ void init(unsigned int seed, curandState_t *states) {

  int tid = blockDim.x*blockIdx.x + threadIdx.x;
  // initialize the state 
  curand_init(seed,       // the seed 
	      tid,
	      0,          // How much extra in the sequence per call
	      &states[tid]);
}


//--------------------------------------------------------------------------
//Takes an array of states and puts a random number into the array
//--------------------------------------------------------------------------
__global__ void randoms(curandState_t *states, double *rands) {
  int tid = blockDim.x*blockIdx.x + threadIdx.x;
  rands[tid] = curand_uniform(&states[tid]);  
}


//--------------------------------------------------------------------------
//Takes an array of rands and computes candidate phi for metro
//--------------------------------------------------------------------------
__global__ void candidatePhi(curandState_t *states, double* rands,
			     double *phi, const double delta_phi) {
  int tid = blockDim.x*blockIdx.x + threadIdx.x;  
  rands[tid] = phi[tid] + delta_phi*(2.0*curand_uniform(&states[tid]) - 1.0);
}


//--------------------------------------------------------------------------
//Device function to compute a (custom) coupling function
//--------------------------------------------------------------------------
__device__ double couplingFunction(double dth, double dt, double sigma, int S1){
  
  dth *= M_PI/S1;  
  dt  *= M_PI/S1;
  
  return pow((cosh(dt) - cos(dth)) , -(1+sigma/2));
  
}


//--------------------------------------------------------------------------
//Add new sites to the Wolff cluster, given a newly added site i
//--------------------------------------------------------------------------
__global__ void clusterAdd(double *gpu_rands, double *gpu_phi, int *gpu_s,
			   double phi_i, int S1, int Lt, int t1, int x1,
			   int cSpin, bool *gpu_added,
			   double sigma, bool usePow, double couplingNorm) {

  int bid = blockIdx.x;
  int tid = threadIdx.x;
  int bdim = blockDim.x;  
  int idx = bid * bdim + tid;
  
  double coupling = 0;

  if(gpu_s[idx] == cSpin) {

    int t2,x2;
    double dt,dx;
    
    //Index divided by circumference, using the int floor feature/bug,
    //gives the timeslice index.
    t2 = idx / S1;
    dt = abs(t2-t1) > Lt/2 ? Lt - abs(t2-t1) : abs(t2-t1);
    
    //The index modulo the circumference gives the spatial index.
    x2 = idx % S1;            
    dx = abs(x2-x1) > S1/2 ? S1 - abs(x2-x1) : abs(x2-x1);      
    
    if(usePow) coupling = pow(sqrt(dx*dx + dt*dt), -(2+sigma));
    else coupling = couplingNorm*couplingFunction(dx, dt, sigma, S1);
    
#ifdef DEBUG
    printf("(TESTING) Hello from block %d, thread %d : %f %f %f %f %f\n",
	   blockIdx.x, threadIdx.x, gpu_rands[idx], phi_i, gpu_phi[idx],
	   coupling, 1 - exp(-2*fabs(phi_i*gpu_phi[idx])*coupling));
#endif
    
    if(gpu_rands[idx] < 1 - exp(-2*fabs(phi_i*gpu_phi[idx])*coupling)) {
      
#ifdef DEBUG
      printf("(ADDING) Hello from block %d, thread %d : %f %f %f %f %f\n",
	     blockIdx.x, threadIdx.x, gpu_rands[idx], phi_i, gpu_phi[idx],
	     coupling, 1 - exp(-2*fabs(phi_i*gpu_phi[idx])*coupling));
#endif      
      gpu_added[idx] = true;
      gpu_s[idx] *= -1;
      gpu_phi[idx] *= -1;
    }
  }  
}


//--------------------------------------------------------------------------
//Perform a single Wolff update
//--------------------------------------------------------------------------
void MonteCarlo2DIsing::GPU_wolffUpdateLR(Param p, int i, int iter, int n) {

  int S1 = p.S1;
  int Lt = p.Lt;
  int vol = p.surfaceVol;
  int t1,x1;
  
  gpu_wc_calls++;
  gpu_wc_s_size = 0;
  
  //Host-chosen random spin.
  int cSpin = s[i];

  //The site belongs to the cluster, so flip it.
  s[i] *= -1;
  phi[i] *= -1.0;
  
  //(re)initialise the the array of sites to check.  
  for(int j=0; j<vol; j++) cpu_added[j] = false;
  //cudaMemset(gpu_added, false, vol*sizeof(bool));
  
  //We start with the seed only.
  cpu_added[i] = true;

  //copy the phi and spin arrays to the GPU.
  //if(n==0) {
  //GPU_copyArraysToDevice(p);
  cudaMemcpy(gpu_phi+i, phi+i, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_s+i, s+i, sizeof(int), cudaMemcpyHostToDevice);
  // }

  int factor = 16;
  
  //This boolean keeps track of whether or not any sites are left to be
  //checked. If a loop over all the lattice shows that no sites are
  //left to be checked, the routine halts and the Wolff cluster is defined.
  bool internalCheck;

  //This boolean is depenedent on internal check. We set it to true so that
  //the while loop actually starts. Once `internalCheck` is false,
  //Check will evaluate to false.
  bool Check = true;
  while(Check) {
    internalCheck = false;
    for(int j=0; j<vol; j++) {      
      if(cpu_added[j] == true ) {
	internalCheck = true;
	cpu_added[j] = false;
	gpu_wc_size++;
#ifdef DEBUG
	printf("Checking site %d\n\n\n",j);
#endif
	
	//This implementation parallelises over the entire surface of the lattice.
	//It performs a boolean check to see if the candidate site has the
	//correct spin, then performs the probablisitic test to add the site.

	t1 = j / S1;
	x1 = j % S1;
	
	//Get some random numbers 
	randoms<<<factor*vol/sze, sze/factor>>>(states, gpu_rands);
	
	cudaMemcpy(gpu_added, cpu_added, vol*sizeof(bool), cudaMemcpyHostToDevice);
	clusterAdd<<<vol/sze, sze>>>(gpu_rands, gpu_phi, gpu_s, phi[j],
				     S1, Lt, t1, x1, cSpin, gpu_added,
				     p.sigma, p.usePowLaw, LR_couplings[0]);  
	cudaMemcpy(cpu_added, gpu_added, vol*sizeof(bool), cudaMemcpyDeviceToHost);
      }
    }
    if(internalCheck == false) Check = false;
    //printf("Size = %d Check = %d ", gpu_wc_s_size, Check);    
  }

  //copy the phi and spin arrays to the host.
  GPU_copyArraysToHost(p);
  
  if( iter%p.n_skip == 0 && n == 0) {
    printf("Average (GPU) cluster size at iter %d,  %d/%d = %f = %f%\n",
	   iter, gpu_wc_size, gpu_wc_calls, 1.0*gpu_wc_size/(gpu_wc_calls), 100.0*gpu_wc_size/(gpu_wc_calls*p.surfaceVol));
  }
}


//--------------------------------------------------------------------------
//Compute the change in KE and perform global sum
//--------------------------------------------------------------------------
__global__ void KEreduce(int i, const int S1, const int Lt,
			 const double sigma, const bool usePow,
			 const double couplingNorm,
			 const double lambda_p, const double musqr_p, 
			 double *gpu_phi,
			 double *gpu_rands_aux,
			 bool *gpu_added,
			 double *result,
			 cudaTextureObject_t tex_phi,
			 cudaTextureObject_t tex_cands){ 
  
  int t1,x1,t2,x2;
  double dx,dt;
  double phi_lc;
  double coupling;
  __shared__ double sdata[sze];
  
  int bid = blockIdx.x;
  int tid = threadIdx.x;
  int bdim = blockDim.x;
  int idx = bid*bdim + tid;  

  sdata[tid] = 0.0;
  
  //printf("%d %d %d %d\n", bid, tid, bdim, idx);

  __shared__ double phi_i;
  __shared__ double phi_new;
  __shared__ double pmpn;
  __shared__ double pnsmps;

  //int2 rval;
  //rval = tex1Dfetch<int2>(tex_phi, i);
  //phi_i = fetch_double(rval);
  //rval = tex1Dfetch<int2>(tex_cands, i);
  //phi_new = fetch_double(rval);

  if(gpu_added[idx]) {
    //rval = tex1Dfetch<int2>(tex_cands, idx);
    //phi_lc = fetch_double(rval);
    phi_lc = gpu_rands_aux[idx];
    //if(idx < 10 && i < 10) printf("GPU idx = %d phi[idx]=%f %f\n", idx, phi_lc, gpu_rands_aux[idx]);
  }
  else {
    //rval = tex1Dfetch<int2>(tex_phi, idx);
    //phi_lc = fetch_double(rval);
    phi_lc = gpu_phi[idx];
    //if(idx < 10 && i < 10) printf("GPU idx = %d phi[idx]=%f %f\n", idx, phi_lc, gpu_phi[idx]);
  }

  if(tid == 0) {
    phi_i = gpu_phi[i];
    phi_new = gpu_rands_aux[i];
    pmpn = phi_i - phi_new;
    pnsmps = 0.5*(phi_new*phi_new - phi_i*phi_i);
    __syncthreads();
  }
  
  t1 = i/S1;
  x1 = i%S1;    
    
  if(idx != i) {
    //Index divided by circumference, using the int floor feature/bug,
    //gives the timeslice index.
    t2 = idx / S1;
    dt = abs(t2-t1) > Lt/2 ? Lt - abs(t2-t1) : abs(t2-t1);      
    //The index modulo the circumference gives the spatial index.
    x2 = idx % S1;            
    dx = abs(x2-x1) > S1/2 ? S1 - abs(x2-x1) : abs(x2-x1);      
    
    if(usePow) coupling = pow(sqrt(dx*dx + dt*dt), -(2+sigma));
    else coupling = couplingNorm*couplingFunction(dx, dt, sigma, S1);
    
    sdata[tid] = (pmpn*phi_lc + pnsmps)*coupling*coupling;
    
    //if(idx < 10 && i < 10) printf("KE GPU = %.10f, i %d, idx %d\n", sdata[tid], i, idx);
    
#ifdef DEBUG    
    if(idx < 10) printf("KE GPU = %.10f, %d %d, %f %f %f %f %f\n", sdata[tid], i, idx, dx, dt, pmpn, pnsmps, coupling);
#endif
  }
  __syncthreads();
  
  //Do block reduction in shared memory  
  for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
    if (tid < s) sdata[tid] += sdata[tid + s];
    __syncthreads();
  }
  //alt: if(bid != 2**k)
  /*
    int dim = blockDim.x;
    while (dim % 2 == 0) {
    dim /= 2;
    if (tid < dim) sdata[tid] += sdata[tid + dim];
    __syncthreads();
    }
  */
  
  auto g = this_thread_block();
  if (g.thread_rank() == 0) atomicAdd(&result[i], sdata[0]);    
}


//--------------------------------------------------------------------------
//Metro step on the total change in energy
//--------------------------------------------------------------------------
__global__ void KEmetroStep(const double lambda_p, const double musqr_p, const int i, 
			    const double *gpu_rands,
			    const double *gpu_rands_aux,
			    const double *gpu_phi,
			    const double *result,
			    bool *gpu_added,
			    cudaTextureObject_t tex_phi,
			    cudaTextureObject_t tex_cands) {
  
  double DeltaE = 0.0;
  
  //int2 rval = tex1Dfetch<int2>(tex_cands, i);
  //double phi_new = fetch_double(rval);
  //rval = tex1Dfetch<int2>(tex_phi, i);
  //double phi_i = fetch_double(rval);
  
  double phi_i_sq = gpu_phi[i];
  double phi_new_sq = gpu_rands_aux[i];    

  phi_i_sq *= phi_i_sq;
  phi_new_sq *= phi_new_sq;
  
  //KE
  DeltaE += result[i];

  //PE	
  DeltaE += lambda_p*(phi_new_sq*phi_new_sq - phi_i_sq*phi_i_sq);
  DeltaE +=  musqr_p*(phi_new_sq            - phi_i_sq);
  
  //if(i < 10) printf("DeltaE GPU = %.8f, result = %.8f %d\n", DeltaE, *result, i);	
  
  if(DeltaE < 0.0) {
    //printf("Acepted %d\n", i);
    gpu_added[i] = true;
  }  
  else if (gpu_rands[i] < exp(-DeltaE) ) {
    //printf("Acepted %d\n", i);
    gpu_added[i] = true;
  }
}


//--------------------------------------------------------------------------
//Perfrom a full Metropolis sweep
//--------------------------------------------------------------------------
double MonteCarlo2DIsing::GPU_metropolisUpdateLR(Param p, int iter) {

  int Lt = p.Lt;
  int S1 = p.S1;
  int vol = p.surfaceVol;
  
  double lambda_p = 0.25*p.lambda;
  double musqr_p  = 0.50*p.musqr;
  
  double sigma = p.sigma;
  double delta_phi = p.delta_phi;
  bool usePowLaw = p.usePowLaw;
  bool doMetroCheck = p.doMetroCheck;
  double couplingNormTheta = LR_couplings[0];

  //copy the phi array to the GPU if this is the very first metro step.
  if(iter == 0) GPU_copyArraysToDevice(p);
  
  //Get some random numbers
  randoms<<<vol/sze, sze>>>(states, gpu_rands);
  //DEBUG copy to host
  if(doMetroCheck) cudaMemcpy(debug_arr1, gpu_rands, vol*sizeof(double), cudaMemcpyDeviceToHost);
  
  //Store all candidate phi on device in gpu_rands_aux
  candidatePhi<<<vol/sze, sze>>>(states, gpu_rands_aux, gpu_phi, delta_phi);
  //DEBUG copy to host
  if(doMetroCheck) cudaMemcpy(debug_arr2, gpu_rands_aux, vol*sizeof(double), cudaMemcpyDeviceToHost);
    
  //create texture objects
  cudaTextureObject_t tex_phi;
  texAndBind(tex_phi, gpu_phi, vol);
  cudaTextureObject_t tex_cands;
  texAndBind(tex_cands, gpu_rands_aux, vol);
  
  cudaMemset(gpu_added, false, vol*sizeof(bool));
  cudaMemset(gpu_result, 0.0, vol*sizeof(double));
    
  for(int i=0; i<vol; i++) {
    
    //Reduce KE array    
    KEreduce<<<vol/sze, sze>>>(i, S1, Lt,
			       sigma, usePowLaw,
			       couplingNormTheta,
			       lambda_p, musqr_p,
			       gpu_phi,
			       gpu_rands_aux,
			       gpu_added,
			       gpu_result,
			       tex_phi, tex_cands);
    
    Check_CUDA_Error("KEreduce Execution Failed!", i, iter);
    
    //Launch a single thread to perform the Metropolis step.
    KEmetroStep<<<1, 1>>>(lambda_p, musqr_p, i,
			  gpu_rands,
			  gpu_rands_aux,
			  gpu_phi,
			  gpu_result,
			  gpu_added,
			  tex_phi, tex_cands);
    
    Check_CUDA_Error("KEmetroStep Execution Failed!", i, iter);
    
  }
  
  cudaDestroyTextureObject(tex_phi);
  cudaDestroyTextureObject(tex_cands);
  
  //Update the phi array and the spin arrays the device using the added array.
  updatePhiAndSpin<<<vol/sze, sze>>>(gpu_rands_aux, gpu_phi, gpu_s, gpu_added);
  
  cudaMemcpy(cpu_added, gpu_added, vol*sizeof(bool), cudaMemcpyDeviceToHost);
  
  for(int a=0; a<vol; a++) if(cpu_added[a]) gpu_accept++;
  //gpu_accept += count_added<<<vol/sze, sze>>>(gpu_added, sum); 
  gpu_tries += vol;
  
  if( (iter < p.n_therm + p.n_metro_cool) && (iter%p.n_skip == 0 && iter > 0)) {
    printf("At iter %d the (GPU) Metro acceptance rate is (%d/%d) = %.10f\n",
	   iter, gpu_accept, gpu_tries,
	   (double)gpu_accept/(double)gpu_tries);
  }
  //exit(0);
  return 0.0;
  
}
#endif

