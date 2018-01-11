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
TIMESLICES=32
MSQR=$1
LEVELS=$4
SRC_POS=-1
g_MSQR=1.0
g_LATT=1.0
delta_MSQR=0.00
N_SHIFT=1

#Ensure these values are sensible!
#Currently set for testing only.
N_THERM=100000
N_MEAS=250
N_SKIP=5000
N_WOLFF=20
MUSQR=$2
LAMBDA=$3

make

mkdir run_${MSQR}_${MUSQR}_${LAMBDA}_${LEVELS}
cp adsrun run_${MSQR}_${MUSQR}_${LAMBDA}_${LEVELS}/.

COMMAND="./adsrun ${BC} ${CENTRE} ${VERBOSITY} \
	 	  ${MAX_ITER} ${TOL} ${TIMESLICES} ${MSQR} ${delta_MSQR} \
	 	  ${LEVELS} ${SRC_POS} ${g_MSQR} ${g_LATT} ${Q} ${N_SHIFT} \
		  ${N_THERM} ${N_MEAS} ${N_SKIP} ${N_WOLFF} ${MUSQR} ${LAMBDA} "

echo ${COMMAND}

(cd run_${MSQR}_${MUSQR}_${LAMBDA}_${LEVELS};
 mkdir data_dump;
 ./${COMMAND} >& log.log &)
