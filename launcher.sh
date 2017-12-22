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
LEVELS=$1
MASS=$2
SRC_POS=$3
C_MASS=$4
N_LATT=$5
Q=$6
TIMESLICES=$7
SINK_T=$8

LAMBDA=0.0
SCALE=1.0

COMMAND="./adsrun ${BC} ${CENTRE} ${VERBOSITY} \
	 ${MAX_ITER} ${TOL} ${TIMESLICES} ${MASS} ${LAMBDA} ${LEVELS} \
	 ${SRC_POS} ${C_MASS} ${N_LATT} ${Q} ${SCALE}"

echo ${COMMAND}
${COMMAND}

cat ./data_dump/q${Q}_Lev${LEVELS}_T${TIMESLICES}_msqr${MASS}_srct0_*_sinkt${SINK_T}Lev${LEVELS}_${BC_STRING}_${CENTRE_STRING}.dat > lattice.dat

cp plotter.p plotter_t.p

sed -i '' s/_C_MASS_/${C_MASS}/g plotter_t.p
sed -i '' s/_N_LATT_/${N_LATT}/g plotter_t.p
sed -i '' s/_MASS_/${MASS}/g plotter_t.p
sed -i '' s/_Q_/${Q}/g plotter_t.p
sed -i '' s/_LEVELS_/${LEVELS}/g plotter_t.p
sed -i '' s/_CENTRE_STRING_/${CENTRE_STRING}/g plotter_t.p
sed -i '' s/_BC_STRING_/${BC_STRING}/g plotter_t.p
sed -i '' s/_SINK_T_/${SINK_T}/g plotter_t.p
sed -i '' s/_TIMESLICES_/${TIMESLICES}/g plotter_t.p

if [ ${TIMESLICES} -ge 1 ]; then
    sed -i '' s/_d_/2/g plotter_t.p
else
    sed -i '' s/_d_/1/g plotter_t.p
fi
    
gnuplot plotter_t.p

WRITE=$9

echo $WRITE

fname=q${Q}_msqr${MASS}_gmsqr${C_MASS}_gN${N_LATT}_T${TIMESLICES}_dT${SINK_T}

if [ $WRITE -eq 1 ]; then
    epstopdf q${Q}_msqr${MASS}_gmsqr${C_MASS}_gN${N_LATT}_T${TIMESLICES}_dT${SINK_T}.eps --autorotate=All
    mkdir EPS
    mv q${Q}_msqr${MASS}_gmsqr${C_MASS}_gN${N_LATT}_T${TIMESLICES}_dT${SINK_T}.eps ./EPS/.
    mkdir PDFs
    mv q${Q}_msqr${MASS}_gmsqr${C_MASS}_gN${N_LATT}_T${TIMESLICES}_dT${SINK_T}.pdf ./PDFs/.
fi

rm q${Q}_msqr${MASS}_gmsqr${C_MASS}_gN${N_LATT}_T${TIMESLICES}_dT${SINK_T}.eps
rm fit.log
rm lattice.dat
rm plotter_t.p
