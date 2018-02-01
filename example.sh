#!/bin/bash

#One must pass the CL variables to the executables in the specific order
#given. Options are listed as comments above the variable.

#d=Dirichlet, n=Neumann
BC='d'

#v=vertex centred, c=circumcentred.
CENTRE='v'

#v=verbose,q=quiet
VERBOSITY='q'

#sql=Square Local, sqnl=Square Non-Local, ads=AdS Lattice
LATTICE='sq_nonlocal'

Q=7
MAX_ITER=100000
TOL=1e-80
TIMESLICES=21
MSQR=-0.7
LEVELS=2
SRC_POS=-1
g_MSQR=1.0
g_LATT=0.1
delta_MSQR=0.0
N_SHIFT=1

#Ensure these values are sensible!
#Currently set for testing only.
N_THERM=10000
N_MEAS=5000
N_SKIP=1000
N_WOLFF=10
MUSQR=1.0
LAMBDA=0.30
SIGMA=1.33

TWS=1.75

make

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
