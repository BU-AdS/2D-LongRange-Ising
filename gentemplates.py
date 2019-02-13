print(
"#!/bin/bash
export OMP_NUM_THREADS=4
VERBOSITY='q'
LATTICE='2D'
THEORY='ISING'
CLUSTER='WOLFF'
CLUSTER_ARCH='CPU'
COUPLING_TYPE='SR'"+
"J=0.4406867935"+
"H=0.0
SIGMA=1.5
CIRCUMFERENCE=64
TIMESLICES=64
N_THERM=500
N_MEAS=10000
N_SKIP=50
N_JKBLK=100
N_WRITE=100
COMMAND=\"./adsrun --verbosity ${VERBOSITY} --theory ${THEORY} --latType ${LATTICE}
		  --couplingType ${COUPLING_TYPE} --J ${J} --h ${H}
		  --Lt ${TIMESLICES} --S1 ${CIRCUMFERENCE}
                  --nTherm ${N_THERM} --nMeas ${N_MEAS}
      		  --nWrite ${N_WRITE} --nJkBlock ${N_JKBLK} --nSkip ${N_SKIP}
       		  --sigma ${SIGMA} --clusterAlg ${CLUSTER} --clusterArch ${CLUSTER_ARCH} \"
echo \"\"
echo \"Command given:\"
echo ${COMMAND}
echo \"\"
${COMMAND}")
