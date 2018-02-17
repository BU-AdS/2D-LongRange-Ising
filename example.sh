#!/bin/bash

#One must pass the CL variables to the executables in the specific order
#given. Options are listed as comments above the variable.

export OMP_NUM_THREADS=16

#d=Dirichlet, n=Neumann
BC='d'

#v=vertex centred, c=circumcentred.
CENTRE='v'

#v=verbose,q=quiet
VERBOSITY='q'

#2D  = rectangular
#ADS = AdS_{2+1}
LATTICE='2D'

#wolff = Wolff algorithm
#sw    = Swednsen Wang algorithm
CLUSTER='wolff'

#SR = Short range
#LR = Long range
COUPLING='LR'

#POW    = 1/r^a type
#RADIAL = 1/(cosh(t) = cos(theta))
COUPLING_TYPE='POW'

Q=8
LEVELS=2
TIMESLICES=32
CIRCUMFERENCE=32

SRC_POS=-1
MAX_ITER=100000
TOL=1e-80

MSQR=1.0
delta_MSQR=0.0
g_MSQR=1.0
g_LATT=1.0
N_SHIFT=1

#Ensure these values are sensible!
#Currently set for testing only.
N_THERM=500
N_MEAS=5000
N_SKIP=100
N_CLUSTER=4

MUSQR=$2
LAMBDA=1.0
SIGMA=$1

TWS=1.0

make -j 12

rm ads_wisdom

rm -rf data_dump
mkdir data_dump

COMMAND="./adsrun ${BC} ${CENTRE} ${VERBOSITY} ${LATTICE} ${COUPLING} \
	 	  ${MAX_ITER} ${TOL} ${TIMESLICES} ${CIRCUMFERENCE} \
                  ${MSQR} ${delta_MSQR} ${LEVELS} ${SRC_POS} \
                  ${g_MSQR} ${g_LATT} ${Q} ${N_SHIFT} ${N_THERM} \
                  ${N_MEAS} ${N_SKIP} ${N_CLUSTER} ${MUSQR} \
                  ${LAMBDA} ${SIGMA} ${TWS} ${CLUSTER} ${COUPLING_TYPE} "

echo ${COMMAND}

${COMMAND}

#Q=7, L=3, M=1
#N_latt  = 0.09435 C_msqr = 4.383623689

#M=2
#N_latt  = 0.08435 C_msqr = 3.959558488

#M=0.1
#N_latt  = 1.318378533 C_msqr = 85.51228172
