#!/bin/bash

#construct mass

LEV=$1

EXP=0
while [ $EXP -ge 1 ]; do
    MANT=1    
    while [ $MANT -le 9 ]; do
	DP1=0
	while [ $DP1 -le 9 ]; do	    
	    ./launcher.sh ${LEV} ${MANT}.${DP1}00e-0${EXP} -1 1.0 1.0 0 7 0
	    let DP1=DP1+1
	done
	let MANT=MANT+1
    done
    let EXP=EXP-1
done

EXP=0
while [ $EXP -le 0 ]; do
    MANT=1    
    while [ $MANT -le 9 ]; do
	DP1=0
	while [ $DP1 -le 9 ]; do	    
	    ./launcher.sh ${LEV} ${MANT}.${DP1}00e+0${EXP} -1 1.0 1.0 0 7 0
	    let DP1=DP1+1
	done
	let MANT=MANT+1	
    done
    let EXP=EXP+1
done

mkdir criticalValues
mv crits* criticalValues/.
