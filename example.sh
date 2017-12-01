#!/bin/bash

#One must pass the CL variables to the executables in the specific order
#given. Options are listed as comments above the variable.

#d=Dirichlet, n=Neumann
BC='d'

#v=vertex centred, c=circumcentred.
CENTRE='v'

#verbose,quiet
VERBOSITY='verbose'

MAX_ITER=100000
TOL=1e-10
TIMESLICES=1
MASS_SQUARED=0.1
LEVELS=9
SRC_POS=20

./adsrun ${BC} ${CENTRE} ${VERBOSITY} \
	 ${MAX_ITER} ${TOL} ${TIMESLICES} ${MASS_SQUARED} ${LEVELS} \
	 ${SRC_POS} 
