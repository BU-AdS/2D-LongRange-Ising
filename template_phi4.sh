#!/bin/bash

export OMP_NUM_THREADS=8

#Options are listed as comments above the variable.

#v=verbose,q=quiet
VERBOSITY='q'

#2D  = rectangular
#ADS = AdS_{2+1}
LATTICE='2D'

#phi4 or Ising
THEORY='phi4'

#wolff = Wolff algorithm
#sw    = Swednsen Wang algorithm
CLUSTER='WOLFF'

#GPU/CPU
CLUSTER_ARCH='CPU'
METRO_ARCH='CPU'

#Perform the Metropolis step on the CPU with the same random numbers used on the GPU.
METRO_CHECK='true'

#SR = Short range
#POW = 1/|x-y|^{2+sigma} type
#RAD = 1/(cosh(dt) - cos(dtheta))^{1+sigma/2}
COUPLING_TYPE='RAD'

CIRCUMFERENCE=16
TIMESLICES=64

N_METRO_COOL=100
N_THERM=1000
N_MEAS=1000
N_SKIP=50
N_CLUSTER=4
DELTA_PHI=1.5

MUSQR=-1.275
LAMBDA=1.0
SIGMA=10.0

COMMAND="./adsrun --verbosity ${VERBOSITY} --theory ${THEORY} --latType ${LATTICE} 
		  --couplingType ${COUPLING_TYPE} --Lt ${TIMESLICES} 
		  --S1 ${CIRCUMFERENCE} --nTherm ${N_THERM} --nMeas ${N_MEAS} 
		  --nSkip ${N_SKIP} --nCluster ${N_CLUSTER} 
		  --nMetroCool ${N_METRO_COOL} --deltaPhi ${DELTA_PHI} --muSqr ${MUSQR}
		  --lambda ${LAMBDA} --sigma ${SIGMA} --clusterAlg ${CLUSTER} 
                  --clusterArch ${CLUSTER_ARCH} --metroArch ${METRO_ARCH} 
                  --metroCheck ${METRO_CHECK} "

echo ""
echo "Command given:"
echo ${COMMAND}
echo ""

${COMMAND}