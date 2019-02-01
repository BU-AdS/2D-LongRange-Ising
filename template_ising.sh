#!/bin/bash

export OMP_NUM_THREADS=8

#Options are listed as comments above the variable.

#v=verbose,q=quiet
VERBOSITY='q'

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
COUPLING_TYPE='RAD'

#Ising J
J=0.4406867935

#Ising h
H=0.0

SIGMA=1.5

CIRCUMFERENCE=16
TIMESLICES=64
N_THERM=500
N_MEAS=10000
N_SKIP=50
N_JKBLK=100
N_WRITE=100

COMMAND="./long_range --verbosity ${VERBOSITY} --theory ${THEORY} 
		      --couplingType ${COUPLING_TYPE} --J ${J} --h ${H}
		      --Lt ${TIMESLICES} --S1 ${CIRCUMFERENCE}
                      --nTherm ${N_THERM} --nMeas ${N_MEAS} 
      		      --nWrite ${N_WRITE} --nJkBlock ${N_JKBLK} --nSkip ${N_SKIP} 
       		      --sigma ${SIGMA} --clusterAlg ${CLUSTER} --clusterArch ${CLUSTER_ARCH} "

echo ""
echo "Command given:"
echo ${COMMAND}
echo ""

${COMMAND}
