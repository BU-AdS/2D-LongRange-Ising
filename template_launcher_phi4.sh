#!/bin/bash

#This script will create an LxT directory, subdirectories for sigma, then lambda and mu.
#Before using it, ensure that in `template_phi4.sh` you have the following CLI structure:
#MU=$1
#LAMBDA=$2
#SIGMA=$3

L=24
T=96
MU=0
P_MU=2
LAMBDA=1
P_LAMBDA=0
SIGMA=1
P_SIGMA=75

while [ ${P_MU} -le 8 ]; do

    DIR=${L}x${T}/sig${SIGMA}p${P_SIGMA}/L${LAMBDA}p${P_LAMBDA}_Mum${MU}p0${P_MU}/
    mkdir -p ${DIR}
    cp template_phi4.sh ${DIR}/.
    cp ~/2D-LongRange-Ising/adsrun ${DIR}/.
    (cd ${DIR}; ./template_phi4.sh ${T} ${L} -${MU}.0${P_MU} ${LAMBDA}.${P_LAMBDA} ${SIGMA}.${P_SIGMA} >& log.log &)
    let P_MU=P_MU+2

done
