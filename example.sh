#!/bin/bash

#One must pass the CL variables to the executables in the specific order
#given. Options are listed as comments above the variable.

#d=Dirichlet, n=Neumann
BC='d'

#v=vertex centred, c=circumcentred.
CENTRE='v'

#v=verbose,q=quiet
VERBOSITY='q'

#sq_local    = Square Local
#sq_nonlocal = Square Non-Local
#sq_ads      = Square AdS
#ads_local   = AdS Lattice
LATTICE='sq_nonlocal'

Q=8
MAX_ITER=100000
TOL=1e-80
TIMESLICES=32
MSQR=1.0
LEVELS=2
SRC_POS=-1
g_MSQR=1.0
g_LATT=1.0
delta_MSQR=0.0
N_SHIFT=1

#Ensure these values are sensible!
#Currently set for testing only.
N_THERM=2000
N_MEAS=500
N_SKIP=100
N_WOLFF=10
MUSQR=-1.2725
LAMBDA=1.0
SIGMA=$2

TWS=$1

make

rm ads_wisdom

rm -rf data_dump
mkdir data_dump

COMMAND="./adsrun ${BC} ${CENTRE} ${VERBOSITY} ${LATTICE} \
	 	  ${MAX_ITER} ${TOL} ${TIMESLICES} ${MSQR} ${delta_MSQR} \
	 	  ${LEVELS} ${SRC_POS} ${g_MSQR} ${g_LATT} ${Q} ${N_SHIFT} \
		  ${N_THERM} ${N_MEAS} ${N_SKIP} ${N_WOLFF} ${MUSQR} ${LAMBDA} \
		  ${SIGMA} ${TWS}"

echo ${COMMAND}

${COMMAND}

#Q=7, L=3, M=1
#N_latt  = 0.09435 C_msqr = 4.383623689

#M=2
#N_latt  = 0.08435 C_msqr = 3.959558488

#M=0.1
#N_latt  = 1.318378533 C_msqr = 85.51228172
