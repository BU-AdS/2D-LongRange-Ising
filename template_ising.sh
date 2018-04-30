#!/bin/bash

export OMP_NUM_THREADS=8

#Options are listed as comments above the variable.

#v=verbose,q=quiet
VERBOSITY='q'

#2D  = rectangular
#ADS = AdS_{2+1}
LATTICE='2D'

#phi4 or Ising
THEORY='ISING'

#wolff = Wolff algorithm
#sw    = Swednsen Wang algorithm
CLUSTER='WOLFF'

#GPU/CPU
CLUSTER_ARCH='CPU'

#SR = Short range
#POW = 1/|x-y|^{2+sigma} type
#RAD = 1/(cosh(dt) - cos(dtheta))^{1+sigma/2}
COUPLING_TYPE='POW'

#Ising J
#J=0.1844
#J=0.4406867935
J=0.43
#Ising h
H=0.0

SIGMA=200

CIRCUMFERENCE=32
TIMESLICES=96
N_THERM=1000
N_MEAS=10000
N_SKIP=50

COMMAND="./adsrun --verbosity ${VERBOSITY} --theory ${THEORY} --latType ${LATTICE} 
		  --couplingType ${COUPLING_TYPE} --J ${J} --h ${H}
		  --Lt ${TIMESLICES} --S1 ${CIRCUMFERENCE}
                  --nTherm ${N_THERM} --nMeas ${N_MEAS} --nSkip ${N_SKIP} 
   		  --sigma ${SIGMA} --clusterAlg ${CLUSTER} "

echo ""
echo "Command given:"
echo ${COMMAND}
echo ""

${COMMAND}
