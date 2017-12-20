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
TOL=1e-10
TIMESLICES=100
MSQR=-0.2
LEVELS=8
SRC_POS=-1
g_MSQR=1.0
g_LATT=1.0
LAMBDA=0.0
SCALE=1.0
N_SHIFT=400

make

COMMAND="./adsrun ${BC} ${CENTRE} ${VERBOSITY} \
	 	  ${MAX_ITER} ${TOL} ${TIMESLICES} ${MSQR} ${LAMBDA} \
	 	  ${LEVELS} ${SRC_POS} ${g_MSQR} ${g_LATT} ${Q} ${N_SHIFT}"

echo ${COMMAND}

${COMMAND}
