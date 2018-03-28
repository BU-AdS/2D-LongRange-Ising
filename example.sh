#!/bin/bash

#One must pass the CL variables to the executables in the specific order
#given. Options are listed as comments above the variable.

export OMP_NUM_THREADS=8

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
METRO_CHECK='false'

#wolff = Wolff algorithm
#sw    = Swednsen Wang algorithm
CLUSTER='WOLFF'

#SR = Short range
#POW = 1/|x-y|^{2+sigma} type
#RAD = 1/(cosh(dt) - cos(dtheta))^{1+sigma/2}
COUPLING_TYPE='POW'

NBLOCKS=$3
NTHREADS=$4
TIMESLICES=$1
CIRCUMFERENCE=$2

N_METRO_COOL=0
N_THERM=1000
N_MEAS=100
N_SKIP=100
N_CLUSTER=4

MUSQR=-1.285
LAMBDA=1.0
SIGMA=1000

make -j 12



#COMMAND="nvprof --metrics all -o metrics_512.nvvp ./adsrun --BC ${BC} --centre ${CENTRE} --verbosity ${VERBOSITY}
COMMAND="./adsrun --BC ${BC} --centre ${CENTRE} --verbosity ${VERBOSITY} 
                  --latType ${LATTICE} --couplingType ${COUPLING_TYPE} \
		  --Lt ${TIMESLICES} --S1 ${CIRCUMFERENCE} 
                  --nTherm ${N_THERM} --nMeas ${N_MEAS} --nSkip ${N_SKIP} 
                  --nCluster ${N_CLUSTER} --nMetroCool ${N_METRO_COOL} 
                  --muSqr ${MUSQR} --lambda ${LAMBDA} --sigma ${SIGMA} --clusterAlg ${CLUSTER} \
                  --clusterArch ${CLUSTER_ARCH} --metroArch ${METRO_ARCH} \
                  --metroCheck ${METRO_CHECK} "
#--nBlocks ${NBLOCKS} --nThreads ${NTHREADS} "


echo ""
echo "Command given:"
echo ${COMMAND}
echo ""

${COMMAND}
