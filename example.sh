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
MAX_ITER=1000
TOL=1e-10
TIMESLICES=1
MASS_SQUARED=1.0
LEVELS=5
SRC_POS=0
C_MSQR=2.0
N_LATT=0.5
LAMBDA=0.0

./adsrun ${BC} ${CENTRE} ${VERBOSITY} \
	 ${MAX_ITER} ${TOL} ${TIMESLICES} ${MASS_SQUARED} ${LAMBDA} ${LEVELS} \
	 ${SRC_POS} ${C_MSQR} ${N_LATT} ${Q}
