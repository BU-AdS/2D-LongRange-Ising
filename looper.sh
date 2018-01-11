#!/bin/bash

MSQR="0"
MUSQR=15
LAMBDA=1.0
LEV=3

while [ ${MUSQR} -le 25 ]; do
    
    ./searcher.sh 0.${MSQR} -0.${MUSQR} ${LAMBDA} ${LEV} >& log.${MSQR}_${MUSQR}_${LAMBDA}_${LEV} &    

    let MUSQR=MUSQR+2
    
done
