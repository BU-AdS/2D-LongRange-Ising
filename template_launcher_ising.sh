#!/bin/bash

#This script will create an LxT directory, subdirectories for sigma, then lambda and mu.
#Before using it, ensure that in `template_ising.sh` you have the following CLI structure:
#L=$1
#T=$2
#J=$3
#SIGMA=$4

L=24
T=96
MU=0
J=0
P_J=2269
SIGMA=1
P_SIGMA=75

while [ ${P_J} -le 2270 ]; do

    DIR=${L}x${T}/sig${SIGMA}p${P_SIGMA}/J${p}p${P_J}
    mkdir -p ${DIR}
    cp template_ising.sh ${DIR}/.
    cp ~/2D-LongRange-Ising/adsrun ${DIR}/.
    (cd ${DIR}; ./template_ising.sh ${L} ${T} ${J}.${P_J} ${SIGMA}.${P_SIGMA} >& log.log &)
    let P_J=P_J+2

done
