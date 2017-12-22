#!/bin/bash

#construct mass

LEV=$1
Q=$2
T=$3
T_S=$4

rm analysis/PDFs/*
rm analysis/EPS/*


EXP=2
MANT=-9    
while [ $MANT -le 9 ]; do
    ./launcher.sh ${LEV} ${MANT}.000e-0${EXP} -1 1.0 1.0 ${Q} ${T} ${T_S} 1
    let MANT=MANT+1
done

EXP=1
MANT=-1
while [ $MANT -le 9 ]; do
    DP1=0
    while [ $DP1 -le 9 ]; do	    
	./launcher.sh ${LEV} ${MANT}.${DP1}00e-0${EXP} -1 1.0 1.0 ${Q} ${T} ${T_S} 1
	let DP1=DP1+1
    done
    let MANT=MANT+1
done

EXP=0
while [ $EXP -le 0 ]; do
    MANT=1
    while [ $MANT -le 4 ]; do
	DP1=0
	while [ $DP1 -le 9 ]; do
	    DP2=0
	    while [ $DP2 -le 9 ]; do
		./launcher.sh ${LEV} ${MANT}.${DP1}${DP2}0e+0${EXP} -1 0.0 1.0 ${Q} ${T} ${T_S} 1
		let DP2=DP2+1
	    done
	    let DP1=DP1+1
	done
	let MANT=MANT+1	
    done
    let EXP=EXP+1
done

mkdir criticalValues
mv crits* criticalValues/.
