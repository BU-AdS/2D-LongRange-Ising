#!/bin/bash

#One must pass the CL variables to the executables in the specific order
#given. Options are listed as comments above the variable.

export OMP_NUM_THREADS=1

#d=Dirichlet, n=Neumann
BC='d'

#v=vertex centred, c=circumcentred.
CENTRE='v'

#v=verbose,q=quiet
VERBOSITY='q'

#2D  = rectangular
#ADS = AdS_{2+1}
LATTICE='2D'

#wolff = Wolff algorithm
#sw    = Swednsen Wang algorithm
CLUSTER='WOLFF'

#GPU/CPU
CLUSTER_ARCH='GPU'
METRO_ARCH='GPU'
#Perform the Metropolis step on the CPU with the same random numbers
#used on the GPU.
METRO_CHECK='false'

#wolff = Wolff algorithm
#sw    = Swednsen Wang algorithm
CLUSTER='WOLFF'

#SR = Short range
#POW = 1/|x-y|^{2+sigma} type
#RAD = 1/(cosh(dt) - cos(dtheta))^{1+sigma/2}
COUPLING_TYPE='rad'

TIMESLICES=$1
CIRCUMFERENCE=$2
T_SCALE=$3

N_METRO_COOL=100
N_THERM=500
N_MEAS=10000
N_SKIP=20
N_CLUSTER=4
DELTA_PHI=1.5

MUSQR=$5
LAMBDA=1.0
SIGMA=$4

make -j 12



#COMMAND="nvprof --metrics all -o metrics_512.nvvp ./adsrun --BC ${BC} --centre ${CENTRE} --verbosity ${VERBOSITY}
COMMAND="./adsrun --BC ${BC} --centre ${CENTRE} --verbosity ${VERBOSITY} 
                  --latType ${LATTICE} --couplingType ${COUPLING_TYPE} 
		  --Lt ${TIMESLICES} --S1 ${CIRCUMFERENCE} --tScale ${T_SCALE}
                  --nTherm ${N_THERM} --nMeas ${N_MEAS} --nSkip ${N_SKIP} 
                  --nCluster ${N_CLUSTER} --nMetroCool ${N_METRO_COOL} --deltaPhi ${DELTA_PHI} 
                  --muSqr ${MUSQR} --lambda ${LAMBDA} --sigma ${SIGMA} --clusterAlg ${CLUSTER} 
                  --clusterArch ${CLUSTER_ARCH} --metroArch ${METRO_ARCH} 
                  --metroCheck ${METRO_CHECK} "
#--nBlocks ${NBLOCKS} --nThreads ${NTHREADS} "


echo ""
echo "Command given:"
echo ${COMMAND}
echo ""

${COMMAND}
