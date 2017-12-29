#!/bin/bash

#One must pass the CL variables to the executables in the specific order
#given. Options are listed as comments above the variable.

#d=Dirichlet, n=Neumann
BC='d'

#v=vertex centred, c=circumcentred.
CENTRE='v'

#v=verbose,q=quiet
VERBOSITY='q'

Q=7
MAX_ITER=100000
TOL=1e-20
TIMESLICES=50
MSQR=1.0
LEVELS=3
SRC_POS=-1
g_MSQR=1.0
g_LATT=1.0
delta_MSQR=0.01
N_SHIFT=1

N_THERM=500000
N_MEAS=100
N_SKIP=10000
MUSQR=0.7
LAMBDA=0.5

make

rm -rf data_dump
mkdir data_dump

COMMAND="./adsrun ${BC} ${CENTRE} ${VERBOSITY} \
	 	  ${MAX_ITER} ${TOL} ${TIMESLICES} ${MSQR} ${delta_MSQR} \
	 	  ${LEVELS} ${SRC_POS} ${g_MSQR} ${g_LATT} ${Q} ${N_SHIFT} \
		  ${N_THERM} ${N_MEAS} ${N_SKIP} ${MUSQR} ${LAMBDA} "

echo ${COMMAND}

${COMMAND}
