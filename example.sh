#!/bin/bash

#One must pass the CL variables to the executables in the specific order
#given. Options are listed as comments above the variable.

export OMP_NUM_THREADS=4

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
CLUSTER='WOLFF'

#SR = Short range
#POW = 1/|x-y}^{2+sigma} type
#RAD = 1/(cosh(dt) - cos(dtheta))^{1+sigma/2}
COUPLING_TYPE='RAD'

TIMESLICES=16
CIRCUMFERENCE=16
#Ensure these values are sensible!
#Currently set for testing only.
N_THERM=1000
N_MEAS=5000
N_SKIP=100
N_CLUSTER=4

MUSQR=$2
LAMBDA=1.0
SIGMA=$1
TWS=1.0

Q=7
LEVELS=3
SRC_POS=-1
MAX_ITER=100000
TOL=1e-16

MSQR=$3
delta_MSQR=0.0
g_MSQR=1.0
g_LATT=1.0
N_SHIFT=1

make -j 12

COMMAND="./adsrun --BC ${BC} --centre ${CENTRE} --verbosity ${VERBOSITY} --latType ${LATTICE} --couplingType ${COUPLING_TYPE} \
		  --maxIter ${MAX_ITER} --tol ${TOL} --Lt ${TIMESLICES} --S1 ${CIRCUMFERENCE} \
                  --mSqr ${MSQR} --deltaMsqr ${delta_MSQR} --levels ${LEVELS} --srcPos ${SRC_POS} \
                  --cMsqr ${cMSQR} --cLatt ${cLATT} --q ${Q} --nShift ${N_SHIFT} --nTherm ${N_THERM} \
                  --nMeas ${N_MEAS} --nSkip ${N_SKIP} --nCluster ${N_CLUSTER} --muSqr ${MUSQR} \
                  --lambda ${LAMBDA} --sigma ${SIGMA} --tScale ${TWS} --clusterAlg ${CLUSTER} "

echo ""
echo "Command  given:"
echo ${COMMAND}
echo ""

${COMMAND}
