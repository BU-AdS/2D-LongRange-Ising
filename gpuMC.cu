#ifdef USE_GPU
//#define DEBUG
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <cooperative_groups.h>
#include <stdio.h>

#include <gpuMC.cuh>
#include <phiFourth2D.h>

int gpu_wc_ave = 0;
int gpu_wc_size = 0;
int gpu_wc_calls = 0;
int gpu_wc_poss = 0;
int gpu_accept = 0;
int gpu_tries = 0;

const int sze = 32;
extern int CACHE_LINE_SIZE;
int chunk = CACHE_LINE_SIZE/sizeof(double);

using namespace cooperative_groups;

void gpu_corr_wolffClusterAddLR(int i, int *s, int cSpin, double *LR_couplings,
				double *phi, bool *cluster, Param p);

//--------------------------------------------------------------------------
//Simple CUDA error parser
//--------------------------------------------------------------------------
void Check_CUDA_Error(const char *message, int site, int iter) {
  cudaError_t error = cudaGetLastError();
  if(error!=cudaSuccess) {
    fprintf(stderr, "ERROR: %s at site %d iter %d: %s\n", message, site,
	    iter, cudaGetErrorString(error) );
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
//Device function to emulate double prec atomic adds
//--------------------------------------------------------------------------
__device__ double atomicAddCustom(double* address, double val) {
  unsigned long long int *address_as_ull = (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
		    __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
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
  if(tid == 0) {
    for(int a=0; a<bdim; a++) {
      if(added[bid*bdim + a]) sdata[a] = 1;
      else sdata[a] = 0;
    }
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
void PhiFourth2D::GPU_copyArraysToHost(Param p) {
  
  int vol = p.S1*p.Lt;
  cudaMemcpy(s, gpu_s, vol*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(cpu_added, gpu_added, vol*sizeof(bool), cudaMemcpyDeviceToHost);
  cudaMemcpy(phi, gpu_phi, vol*sizeof(double), cudaMemcpyDeviceToHost);
  
}


//--------------------------------------------------------------------------
//Copy arrays from host to device
//--------------------------------------------------------------------------
void PhiFourth2D::GPU_copyArraysToDevice(Param p) {
  
  int vol = p.S1*p.Lt;
  cudaMemcpy(gpu_s, s, vol*sizeof(int), cudaMemcpyHostToDevice);  
  cudaMemcpy(gpu_added, cpu_added, vol*sizeof(bool), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_phi, phi, vol*sizeof(double), cudaMemcpyHostToDevice);
}

//--------------------------------------------------------------------------
//Copy arrays from device to host
//--------------------------------------------------------------------------
void Ising2D::GPU_copyArraysToHostI(Param p) {
  
  int vol = p.S1*p.Lt;
  cudaMemcpy(s, gpu_s, vol*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(cpu_added, gpu_added, vol*sizeof(bool), cudaMemcpyDeviceToHost);
}


//--------------------------------------------------------------------------
//Copy arrays from host to device
//--------------------------------------------------------------------------
void Ising2D::GPU_copyArraysToDeviceI(Param p) {
  
  int vol = p.S1*p.Lt;
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
void Ising2D::GPU_initRand(Param p, int seed, curandState_t* states) {
  
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
	      20,          // How much extra in the sequence per call
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
__global__ void clusterAddILR(const double *gpu_rands, const double *gpu_isingProb,
			      int *gpu_s,
			      const int S1, const int Lt, const int t1, const int x1,
			      const double J, const double t_scale, const int cSpin,
			      bool *gpu_added,
			      const double sigma, const bool usePow,
			      const double couplingNorm) {
  
  int bid = blockIdx.x;
  int tid = threadIdx.x;
  int bdim = blockDim.x;  
  int idx = bid * bdim + tid;
  
  double coupling = 0;

  if(gpu_s[idx]*cSpin > 0) {

    int t2,x2,dt,dx;
    
    //Index divided by circumference, using the int floor feature/bug,
    //gives the timeslice index.
    t2 = idx / S1;
    dt = abs(t2-t1) > Lt/2 ? Lt - abs(t2-t1) : abs(t2-t1);
    
    //The index modulo the circumference gives the spatial index.
    x2 = idx % S1;            
    dx = abs(x2-x1) > S1/2 ? S1 - abs(x2-x1) : abs(x2-x1);      
    
    //if(usePow) coupling = pow(sqrt(dx*dx + dt*dt), -(2+sigma));
    //else coupling = couplingNorm*couplingFunction(dx, t_scale*dt, sigma, S1);
    
#ifdef DEBUG
    printf("(TESTING) block %d, thread %d, idx %d : %d %d %f %f %f\n",
	   blockIdx.x, threadIdx.x, idx, cSpin, gpu_s[idx], gpu_rands[idx], coupling, 1 - exp(-2.0*J*coupling));
#endif
    
    //if(gpu_rands[idx] < (1 - exp(-2.0*J*coupling))) {
    if(gpu_rands[idx] < gpu_isingProb[dx + (S1/2+1)*dt]) {
      
#ifdef DEBUG
      printf("(ADDING) block %d, thread %d, idx %d : %d %d %f %f %f\n",
	     blockIdx.x, threadIdx.x, idx, cSpin, gpu_s[idx], gpu_rands[idx], coupling, 1 - exp(-2.0*J*coupling));
#endif      
      gpu_added[idx] = true;
      gpu_s[idx] *= -1;
    }
  }  
}


//--------------------------------------------------------------------------
//Perform a single Wolff update
//--------------------------------------------------------------------------
void Ising2D::GPU_wolffUpdateILR(Param p, int i, int iter) {

  int S1 = p.S1;
  int Lt = p.Lt;
  int vol = p.surfaceVol;
  int t1,x1;
  double t_scale = p.t_scale;
  double J = p.J;

  //copy the spin array to the GPU if first iteration.
  if(iter == 0) {
    printf("\nCopying Host spin array to device\n\n"); 
    cudaMemcpy(gpu_s, s, vol*sizeof(int), cudaMemcpyHostToDevice);
  }
  
  //Reset averages to exclude thermalisation data.
  if(iter == p.n_therm) {
    printf("\n\nResetting Cluster Stats.\n\n");
    gpu_wc_calls = 0;
    gpu_wc_size = 0;
  }    
  
  gpu_wc_calls++;
  
  //Host-chosen random spin.
  int cSpin = s[i];

  //The site belongs to the cluster, so flip it.
  s[i] *= -1;
  //Mirror the flip on the device
  cudaMemcpy(gpu_s+i, s+i, sizeof(int), cudaMemcpyHostToDevice);
  
  //(re)initialise the the array of sites to check.  
  for(int j=0; j<vol; j++) cpu_added[j] = false;
  
  //We start with the seed only.
  cpu_added[i] = true;
  
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
	//printf("\nChecking site %d\n\n",j);
	
	//This implementation parallelises over the entire surface of the lattice.
	//It performs a boolean check to see if the candidate site has the
	//correct spin, then performs the probablisitic test to add the site.

	t1 = j / S1;
	x1 = j % S1;
	
	//Get some random numbers 
	randoms<<<vol/sze, sze>>>(states, gpu_rands);
	
	cudaMemcpy(gpu_added, cpu_added, vol*sizeof(bool), cudaMemcpyHostToDevice);
	clusterAddILR<<<vol/sze, sze>>>(gpu_rands, gpu_isingProb, gpu_s, S1, Lt, t1, x1,
					J, t_scale, cSpin, gpu_added,
					p.sigma, p.usePowLaw, LR_couplings[0]); 
	cudaMemcpy(cpu_added, gpu_added, vol*sizeof(bool), cudaMemcpyDeviceToHost);
      }
    }
    if(internalCheck == false) Check = false;
  }
  cudaMemcpy(s, gpu_s, vol*sizeof(int), cudaMemcpyDeviceToHost);
  
  if( iter%p.n_skip == 0 ) {
    printf("Average (GPU) cluster size at iter %d,  %d/%d = %f = %f% \n",
	   iter, gpu_wc_size, gpu_wc_calls, (1.0*gpu_wc_size)/gpu_wc_calls, (100.0*gpu_wc_size)/(gpu_wc_calls*p.surfaceVol) );
  }
}




//--------------------------------------------------------------------------
//Add new sites to the Wolff cluster, given a newly added site i
//--------------------------------------------------------------------------
__global__ void clusterAddLR(const double *gpu_rands,
			     double *gpu_phi, int *gpu_s,
			     const double phi_i, const int S1, const int Lt,
			     const int t1, const int x1, const double t_scale,
			     const int cSpin, bool *gpu_added,
			     const double sigma, const bool usePow,
			     const double couplingNorm) {

  int bid = blockIdx.x;
  int tid = threadIdx.x;
  int bdim = blockDim.x;  
  int idx = bid * bdim + tid;
  
  double coupling = 0;

  if(gpu_phi[idx]*cSpin > 0) {

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
    else coupling = couplingNorm*couplingFunction(dx, t_scale*dt, sigma, S1);
    
#ifdef DEBUG
    printf("(TESTING) block %d, thread %d, idx %d cSpin %d : %f %f %f %f %f\n",
	   blockIdx.x, threadIdx.x, idx, cSpin, gpu_rands[idx], phi_i, gpu_phi[idx],
	   coupling, 1 - exp(-2*fabs(phi_i*gpu_phi[idx])*coupling));
#endif
    
    if(gpu_rands[idx] < (1 - exp(-2*fabs(phi_i*gpu_phi[idx])*coupling)) ) {
      
#ifdef DEBUG
      printf("(ADDING) block %d, thread %d, idx %d cSpin %d : %f %f %f %f %f\n",
	     blockIdx.x, threadIdx.x, idx, cSpin, gpu_rands[idx], phi_i, gpu_phi[idx],
	     coupling, 1 - exp(-2*fabs(phi_i*gpu_phi[idx])*coupling));
#endif      
      gpu_added[idx] = true;
      gpu_phi[idx] *= -1;
      gpu_s[idx] *= -1;
    }
  }  
}


//--------------------------------------------------------------------------
//Perform a single Wolff update
//--------------------------------------------------------------------------
void PhiFourth2D::GPU_wolffUpdateLR(Param p, int i, int iter, int n) {

  int S1 = p.S1;
  int Lt = p.Lt;
  int vol = p.surfaceVol;
  int t1,x1;
  double t_scale = p.t_scale;

  //Reset averages to exclude thermalisation.
  if(iter == p.n_therm && n == 0) {
    printf("\n\nResetting Cluster Stats.\n\n");
    gpu_wc_calls = 0;
    gpu_wc_size = 0;
  }    
  
  gpu_wc_calls++;
  
  //Host-chosen random spin.
  int cSpin = s[i];

  //The site belongs to the cluster, so flip it.
  s[i] *= -1;
  phi[i] *= -1.0;
  //Mirror the flip on the device
  cudaMemcpy(gpu_phi+i, phi+i, sizeof(double), cudaMemcpyHostToDevice);
  cudaDeviceSynchronize();
  cudaMemcpy(gpu_s+i, s+i, sizeof(int), cudaMemcpyHostToDevice);
  cudaDeviceSynchronize();
  
  //(re)initialise the the array of sites to check.  
  for(int j=0; j<vol; j++) cpu_added[j] = false;
  
  //We start with the seed only.
  cpu_added[i] = true;

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
      if(cpu_added[j] == true ) {  //We found a newly added site
	internalCheck = true;      //We must therfore continue the algorithm
	cpu_added[j] = false;      //We no longer veiw site j as newly added
	gpu_wc_size++;             //Increase the cluster size
	//printf("Checking site %d\n\n\n",j);
	
	//This implementation parallelises over the entire surface of the lattice.
	//It performs a boolean check to see if the candidate site has the
	//correct spin, then performs the probablisitic test to add the site.

	t1 = j / S1;
	x1 = j % S1;
	
	//Get some random numbers 
	randoms<<<vol/sze, sze>>>(states, gpu_rands);
	
	cudaMemcpy(gpu_added, cpu_added, vol*sizeof(bool), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	
	clusterAddLR<<<vol/sze, sze>>>(gpu_rands, gpu_phi, gpu_s, phi[j],
				       S1, Lt, t1, x1, t_scale, cSpin, gpu_added,
				       p.sigma, p.usePowLaw, LR_couplings[0]);
	cudaDeviceSynchronize();
	
	cudaMemcpy(cpu_added, gpu_added, vol*sizeof(bool), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
      }
    }
    if(internalCheck == false) Check = false;
  }

  //copy the phi and spin arrays to the host.
  GPU_copyArraysToHost(p);
  
  if( iter%p.n_skip == 0 && n == 0) {
    printf("Average (GPU) cluster size at iter %d,  %d/%d = %f = %f% \n",
	   iter, gpu_wc_size, gpu_wc_calls, (1.0*gpu_wc_size)/gpu_wc_calls, (100.0*gpu_wc_size)/(gpu_wc_calls*p.surfaceVol));
  }
}


//--------------------------------------------------------------------------
//Compute the change in KE and perform global sum
//--------------------------------------------------------------------------
__global__ void KEreduce(int i, const int S1, const int Lt,
			 const double sigma, const bool usePow,
			 const double couplingNorm,
			 const double t_scale,
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
    else coupling = couplingNorm*couplingFunction(dx, t_scale*dt, sigma, S1);
    
    sdata[tid] = (pmpn*phi_lc + pnsmps)*coupling*coupling;
    
    //if(idx < 10 && i < 10) printf("KE GPU = %.10f, i %d, idx %d\n", sdata[tid], i, idx);    
    //if(idx < 10) printf("KE GPU = %.10f, %d %d, %f %f %f %f %f\n", sdata[tid], i, idx, dx, dt, pmpn, pnsmps, coupling);
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
  if (g.thread_rank() == 0) atomicAddCustom(&result[i], sdata[0]);
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
double PhiFourth2D::GPU_metropolisUpdateLR(Param p, int iter) {

  int Lt = p.Lt;
  int S1 = p.S1;
  int vol = p.surfaceVol;

  double t_scale = p.t_scale;
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
			       t_scale,
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
  //int *cpu_sum;
  //*cpu_sum = 0;
  //count_added<<<vol/sze, sze>>>(gpu_added, gpu_sum);  
  //cudaMemcpy(cpu_sum, gpu_sum, sizeof(int), cudaMemcpyDeviceToHost);  
  //gpu_accept += *cpu_sum;  
  gpu_tries += vol;
  
  if( (iter < p.n_therm + p.n_metro_cool) && (iter%p.n_skip == 0 && iter > 0)) {
    printf("At iter %d the (GPU) Metro acceptance rate is (%d/%d) = %.10f\n",
	   iter, gpu_accept, gpu_tries,
	   (double)gpu_accept/(double)gpu_tries);
  }
  //exit(0);
  return 0.0;  
}

#if 0

//Correlation Functions.
int gpu_corr_sqnl_wc_size = 0;
int gpu_corr_sqnl_wc_ave = 0;
int gpu_corr_sqnl_wc_calls = 0;

void PhiFourth2D::GPU_correlatorsImpWolff(int i, int meas, double avePhi, Param p){

  int S1 = p.S1;
  int Lt = p.Lt;
  int vol = p.surfaceVol;
  int t1,x1;
  double t_scale = p.t_scale;
  
  gpu_corr_sqnl_wc_calls++;
  gpu_corr_sqnl_wc_size = 0;
  
  //Create copies of the phi and spin arrays for this calculation, and a boolean
  //array to identify the cluster.
  for(int a=0; a<p.surfaceVol; a++) {
    s_cpy[a] = s[a];
    phi_cpy[a] = phi[a];
    cpu_added[a] = false;
  }   
 
  //Host-chosen random spin.
  int cSpin = s_cpy[i];

  // The site belongs to the cluster, so flip it.
  s_cpy[i] *= -1;
  phi_cpy[i] *= -1;

  //Initialise the the array of sites to check.  
  for(int j=0; j<vol; j++) {
    cpu_added[j] = false;
  }
  cudaMemcpy(gpu_cluster, cpu_added, vol*sizeof(bool), cudaMemcpyHostToDevice);
  
  //We start with the seed only.
  cpu_added[i] = true;
  
  //copy the phi and spin arrays to the GPU.
  cudaMemcpy(gpu_phi_cpy, phi_cpy, vol*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_s_cpy, s_cpy, vol*sizeof(int), cudaMemcpyHostToDevice);

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
	gpu_corr_sqnl_wc_size++;
#ifdef DEBUG
	printf("Checking site %d\n\n\n",j);
#endif
	
	//This implementation parallelises over the entire surface of the lattice.
	//It performs a boolean check to see if the candidate site has the
	//correct spin, then performs the probablisitic test to add the site.

	t1 = j / S1;
	x1 = j % S1;
	
	//Get some random numbers 
	randoms<<<vol/sze, sze>>>(states, gpu_rands);
	
	cudaMemcpy(gpu_added, cpu_added, vol*sizeof(bool), cudaMemcpyHostToDevice);
	corrClusterAdd<<<vol/sze, sze>>>(gpu_rands, gpu_phi_cpy, gpu_s_cpy,
					 phi_cpy[j], S1, Lt, t1, x1, t_scale, cSpin,
					 gpu_added, gpu_cluster, p.sigma,
					 p.usePowLaw, LR_couplings[0]);  
	cudaMemcpy(cpu_added, gpu_added, vol*sizeof(bool), cudaMemcpyDeviceToHost);
      }
    }
    if(internalCheck == false) Check = false;
  }

  gpu_corr_sqnl_wc_ave += gpu_corr_sqnl_wc_size;
  
  printf("Average (GPU) Corr. cluster size at measurement %d = %d/%d = %.4f = %.4f%\n",
	 meas+1, gpu_corr_sqnl_wc_ave, gpu_corr_sqnl_wc_calls,
	 1.0*gpu_corr_sqnl_wc_ave/(gpu_corr_sqnl_wc_calls*p.surfaceVol),
	 100.0*gpu_corr_sqnl_wc_ave/(gpu_corr_sqnl_wc_calls*p.surfaceVol));


  //We now have an array 'gpu_cluster' on the GPU that defines a Wolff cluster.
  //We use it to perform the weighted correlation function calculation.

  cudaMemset(gpu_ind_corr_s, 0.0, (S1/2+1)*sizeof(double));
  corrClusterCalcSpatial<<<Lt, S1>>>(gpu_phi_cpy, S1, Lt, gpu_cluster, gpu_ind_corr_s);
  Check_CUDA_Error("corrClusterCalcSpatial Execution Failed!", i, meas);
  
  cudaMemset(gpu_ind_corr_t, 0.0, (Lt/2+1)*sizeof(double));
  corrClusterCalcTemporal<<<S1, Lt>>>(gpu_phi_cpy, S1, Lt, gpu_cluster, gpu_ind_corr_t);  
  Check_CUDA_Error("corrClusterCalcTemporal Execution Failed!", i, meas);

  //copy individual correlation functions to Host, to running sum on host.
  cudaMemcpy(ind_corr_t[meas], gpu_ind_corr_t, (Lt/2+1)*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(ind_corr_s[meas], gpu_ind_corr_s, (S1/2+1)*sizeof(double), cudaMemcpyDeviceToHost);
  
  for(int a=0; a<(S1/2+1); a++) {
    ind_corr_s[meas][a] /= gpu_corr_sqnl_wc_size;
    run_corr_s[a] += ind_corr_s[meas][a];
  }
  
  for(int a=0; a<(Lt/2+1); a++) {
    ind_corr_t[meas][a] /= gpu_corr_sqnl_wc_size;
    run_corr_t[a] += ind_corr_t[meas][a];
  }  
}

//--------------------------------------------------------------------------
//Add new sites to the Wolff cluster, given a newly added site i
//--------------------------------------------------------------------------
__global__ void corrClusterAdd(double *gpu_rands, double *gpu_phi, int *gpu_s,
			       double phi_i, int S1, int Lt, int t1, int x1, double t_scale,
			       int cSpin, bool *gpu_added, bool *gpu_cluster,
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
    else coupling = couplingNorm*couplingFunction(dx, t_scale*dt, sigma, S1);
    
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
      gpu_cluster[idx] = true;
      gpu_s[idx] *= -1;
      gpu_phi[idx] *= -1;
    }
  }  
}

__global__ void corrClusterCalcSpatial(double *gpu_phi_cpy, int S1, int Lt,
				       bool *gpu_cluster, double *ind_corr) {

  int bid = blockIdx.x;
  int tid = threadIdx.x;
  int bdim = blockDim.x;  
  int idx = bid * bdim + tid;
  __shared__ double sink_phi[32];
  double corr_lc[32];
  for(int a=0; a<32; a++) corr_lc[a] = 0.0;
  
  //Check if this site is in the Wolff cluster.
  if(gpu_cluster[idx]) {
        
    int x1,x2,dx;
    int t_slice = S1*bid;
    
    //Use the thread index for the spatial coordinate.
    x1 = tid;
    //Each unique thread now has a unique phi
    double phi_lc = gpu_phi_cpy[idx];
    
    //Each block acquires data in the spatial direction

    sink_phi[tid] = gpu_phi_cpy[idx];
    __syncthreads();
    
    //Each thread loops over all sinks.
    for(x2=0; x2<bdim; x2++) {
      
      //Check if the sink is in the cluster
      if(gpu_cluster[t_slice + x2]) {
	
	//Compute dx
	dx = abs(x2-x1) > S1/2 ? S1 - abs(x2-x1) : abs(x2-x1);            
	
	//Compute \phi * \phi correlation value.
	corr_lc[dx] += phi_lc*sink_phi[x2];
      }
    }
    __syncthreads();
    
    auto g = this_thread_block();
    for(int a=0; a<32; a++) atomicAddCustom(&ind_corr[a], corr_lc[a]);    
  }      
}


__global__ void corrClusterCalcTemporal(double *gpu_phi_cpy, int S1, int Lt,
					bool *gpu_cluster, double *ind_corr) {
  
  int bid = blockIdx.x;
  int tid = threadIdx.x;
  int bdim = blockDim.x;
  int gdim = gridDim.x;  
  int idx = tid * gdim + bid;
  __shared__ double sink_phi[32];
  double corr_lc[32];
  for(int a=0; a<32; a++) corr_lc[a] = 0.0;
  
  //Check if this site is in the Wolff cluster.
  if(gpu_cluster[idx]) {
        
    int t1,t2,dt;
    
    //Use the thread index for the temporal coordinate.
    t1 = tid;
    //Each unique thread now has a unique phi
    double phi_lc = gpu_phi_cpy[idx];
    
    //Each block acquires data in the temporal direction
    sink_phi[tid] = gpu_phi_cpy[idx];
    __syncthreads();
    
    //Each thread loops over all sinks.
    for(t2=0; t2<bdim; t2++) {
      
      //Check if the sink is in the cluster
      if(gpu_cluster[t2*gdim + bid]) {
	
	//Compute dt
	dt = abs(t2-t1) > Lt/2 ? Lt - abs(t2-t1) : abs(t2-t1);            
	
	//Compute \phi * \phi correlation value.
	corr_lc[dt] += phi_lc*sink_phi[t2];
      }
    }
    __syncthreads();

    auto g = this_thread_block();
    if (g.thread_rank() == tid && tid < Lt/2 + 1) atomicAddCustom(&ind_corr[tid], corr_lc[tid]);    
  }
}


#endif

#endif


