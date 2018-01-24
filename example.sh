#!/bin/bash

#One must pass the CL variables to the executables in the specific order
#given. Options are listed as comments above the variable.

#d=Dirichlet, n=Neumann
BC='d'

#v=vertex centred, c=circumcentred.
CENTRE='v'

#v=verbose,q=quiet
VERBOSITY='q'

#sq=Square Lattice, ads=AdS Lattice
LATTICE='ads'

Q=7
MAX_ITER=100000
TOL=1e-30
TIMESLICES=56
MSQR=-0.1
LEVELS=3
SRC_POS=-1
g_MSQR=1.0
g_LATT=1.0
delta_MSQR=0.00
N_SHIFT=1

#Ensure these values are sensible!
#Currently set for testing only.
N_THERM=500000
N_MEAS=5000
N_SKIP=5000
N_WOLFF=10
MUSQR=-9.2
LAMBDA=1.0

make

rm -rf data_dump
mkdir data_dump

COMMAND="./adsrun ${BC} ${CENTRE} ${VERBOSITY} ${LATTICE} \
	 	  ${MAX_ITER} ${TOL} ${TIMESLICES} ${MSQR} ${delta_MSQR} \
	 	  ${LEVELS} ${SRC_POS} ${g_MSQR} ${g_LATT} ${Q} ${N_SHIFT} \
		  ${N_THERM} ${N_MEAS} ${N_SKIP} ${N_WOLFF} ${MUSQR} ${LAMBDA} "

echo ${COMMAND}

${COMMAND}

