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
TOL=1e-10
TIMESLICES=56
MSQR=0.1
LEVELS=3
SRC_POS=-1
g_MSQR=1.0
g_LATT=1.0
delta_MSQR=0.00
N_SHIFT=1

#Ensure these values are sensible!
#Currently set for testing only.
N_THERM=100000
N_MEAS=5000
N_SKIP=5000
N_WOLFF=3
MUSQR=-0.50371
LAMBDA=0.25

make

rm -rf data_dump
mkdir data_dump

COMMAND="./adsrun ${BC} ${CENTRE} ${VERBOSITY} ${LATTICE} \
	 	  ${MAX_ITER} ${TOL} ${TIMESLICES} ${MSQR} ${delta_MSQR} \
	 	  ${LEVELS} ${SRC_POS} ${g_MSQR} ${g_LATT} ${Q} ${N_SHIFT} \
		  ${N_THERM} ${N_MEAS} ${N_SKIP} ${N_WOLFF} ${MUSQR} ${LAMBDA} "

echo ${COMMAND}

${COMMAND}

#<phi(x) - <phi(x)>> * <phi(y) - <phi(y)>> = B / (|x-y|)**delta + D

#Surface = T = 32 

#M_{AdS}=1.0 \mu^2=-0.00 \lambda=1.00 unstable

#M_{AdS}=1.0 \mu^2=-0.30 \lambda=1.00 permafrost!



#M_{AdS}=1.0 \mu^2=-0.16 \lambda=1.00 critical binder = 0.25
#B               = 0.0973079        +/- 0.001401     (1.44%)
#delta           = 0.557269         +/- 0.01088      (1.953%)

#M_{AdS}=1.0 \mu^2=-0.15 \lambda=1.00 critical binder = 0.21
#B               = 0.0924327        +/- 0.001697     (1.836%)
#delta           = 0.630542         +/- 0.01496      (2.372%)
