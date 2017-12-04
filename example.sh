#!/bin/bash

#One must pass the CL variables to the executables in the specific order
#given. Options are listed as comments above the variable.

#d=Dirichlet, n=Neumann
BC='d'

#v=vertex centred, c=circumcentred.
CENTRE='c'

#v=verbose,q=quiet
VERBOSITY='q'

MAX_ITER=1000
TOL=1e-10
TIMESLICES=1
MASS_SQUARED=$1
LEVELS=4
SRC_POS=$2

./adsrun ${BC} ${CENTRE} ${VERBOSITY} \
	 ${MAX_ITER} ${TOL} ${TIMESLICES} ${MASS_SQUARED} ${LEVELS} \
	 ${SRC_POS} 
