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

#wolff = Wolff algorithm
#sw    = Swednsen Wang algorithm
CLUSTER='wolff'

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
N_MEAS=5000
N_SKIP=20
N_CLUSTER=4
MUSQR=-1.2725
LAMBDA=1.0
SIGMA=$2
TWS=$1

make -j 12

rm ads_wisdom

rm -rf data_dump
mkdir data_dump

COMMAND="./adsrun ${BC} ${CENTRE} ${VERBOSITY} ${LATTICE} \
	 	  ${MAX_ITER} ${TOL} ${TIMESLICES} ${MSQR} ${delta_MSQR} \
	 	  ${LEVELS} ${SRC_POS} ${g_MSQR} ${g_LATT} ${Q} ${N_SHIFT} \
		  ${N_THERM} ${N_MEAS} ${N_SKIP} ${N_CLUSTER} ${MUSQR} 
                  ${LAMBDA} ${SIGMA} ${TWS} ${CLUSTER} "

echo ${COMMAND}

${COMMAND}

#Q=7, L=3, M=1
#N_latt  = 0.09435 C_msqr = 4.383623689

#M=2
#N_latt  = 0.08435 C_msqr = 3.959558488

#M=0.1
#N_latt  = 1.318378533 C_msqr = 85.51228172
