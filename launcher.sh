#!/bin/bash

make

#One must pass the CL variables to the executables in the specific order
#given. Options are listed as comments above the variable.

#d=Dirichlet, n=Neumann
BC='d'
BC_STRING="Dirichlet"
if [ "$BC" = "n" ]; then
    BC_STRING="Neumann"
fi

#v=vertex centred, c=circumcentred.
CENTRE='v'
CENTRE_STRING="Vertex"
if [ "$CENTRE" = "c" ]; then
    CENTRE_STRING="Centre"
fi

#v=verbose,q=quiet
VERBOSITY='q'

MAX_ITER=1000
TOL=1e-10
TIMESLICES=100
LEVELS=$1
MASS=$2
SRC_POS=$3
C_MASS=$4
N_LATT=$5
LAMBDA=$6
Q=$7
SCALE=1.0

COMMAND="./adsrun ${BC} ${CENTRE} ${VERBOSITY} \
	 ${MAX_ITER} ${TOL} ${TIMESLICES} ${MASS} ${LAMBDA} ${LEVELS} \
	 ${SRC_POS} ${C_MASS} ${N_LATT} ${Q} ${SCALE}"

echo ${COMMAND}
${COMMAND}

cat q${Q}_Lev${LEVELS}_T${TIMESLICES}_msqr${MASS}_src*_sinkLev${LEVELS}_${BC_STRING}_${CENTRE_STRING}.dat > lattice.dat
rm *${CENTRE_STRING}.dat

cp plotter.p plotter_t.p

sed -i '' s/_C_MASS_/${C_MASS}/g plotter_t.p
sed -i '' s/_N_LATT_/${N_LATT}/g plotter_t.p
sed -i '' s/_MASS_/${MASS}/g plotter_t.p
sed -i '' s/_Q_/${Q}/g plotter_t.p
sed -i '' s/_LEVELS_/${LEVELS}/g plotter_t.p
sed -i '' s/_CENTRE_STRING_/${CENTRE_STRING}/g plotter_t.p
sed -i '' s/_BC_STRING_/${BC_STRING}/g plotter_t.p

gnuplot plotter_t.p
rm plotter_t.p
rm lattice.dat
rm fit.log

WRITE=$8

if [ $WRITE -eq 1 ]; then
    epstopdf q=${Q}_m2=${MASS}_Cm2=${C_MASS}_CN=${N_LATT}.eps --autorotate=All
    rm q=${Q}_m2=${MASS}_Cm2=${C_MASS}_CN=${N_LATT}.eps
    open q=${Q}_m2=${MASS}_Cm2=${C_MASS}_CN=${N_LATT}.pdf
    mkdir PDFs
    mv q=${Q}_m2=${MASS}_Cm2=${C_MASS}_CN=${N_LATT}.pdf PDFs    
fi

rm q=${Q}_m2=${MASS}_Cm2=${C_MASS}_CN=${N_LATT}.eps
