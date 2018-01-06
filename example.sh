#!/bin/bash

#One must pass the CL variables to the executables in the specific order
#given. Options are listed as comments above the variable.

#d=Dirichlet, n=Neumann
BC='d'

#v=vertex centred, c=circumcentred.
CENTRE='v'

#v=verbose,q=quiet
VERBOSITY='q'

Q=8
MAX_ITER=100000
TOL=1e-10
TIMESLICES=160
MSQR=1.0
LEVELS=2
SRC_POS=-1
g_MSQR=1.0
g_LATT=1.0
delta_MSQR=0.00
N_SHIFT=1

#Ensure these values are sensible!
#Currently set for testing only.
N_THERM=50000
N_MEAS=1000
N_SKIP=5000
N_WOLFF=20
MUSQR=-0.1
LAMBDA=0.3

make

rm -rf data_dump
mkdir data_dump

COMMAND="./adsrun ${BC} ${CENTRE} ${VERBOSITY} \
	 	  ${MAX_ITER} ${TOL} ${TIMESLICES} ${MSQR} ${delta_MSQR} \
	 	  ${LEVELS} ${SRC_POS} ${g_MSQR} ${g_LATT} ${Q} ${N_SHIFT} \
		  ${N_THERM} ${N_MEAS} ${N_SKIP} ${N_WOLFF} ${MUSQR} ${LAMBDA} "

echo ${COMMAND}

${COMMAND}

# mu=-0.010, lambda=0.010, mass=0.16742524
# mu=-0.015, lambda=0.015, mass=1.545
# mu=-0.020, lambda=0.020, mass=200.5
