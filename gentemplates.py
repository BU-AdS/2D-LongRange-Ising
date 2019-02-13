import sys

sys.stdout = open('output.sh','wt')

print(
"#!/bin/bash\n"\
"export OMP_NUM_THREADS=4\n"\
"VERBOSITY='q'\n"
"LATTICE='2D'\n"\
"THEORY='ISING'\n"\
"CLUSTER='WOLFF'\n"\
"CLUSTER_ARCH='CPU'\n"\
"COUPLING_TYPE='SR'\n"+\
"J=0.4406867935"+\
"H=0.0\n"\
"SIGMA=1.5\n"\
"CIRCUMFERENCE=64\n"\
"TIMESLICES=64\n"\
"N_THERM=500\n"\
"N_MEAS=10000\n"\
"N_SKIP=50\n"\
"N_JKBLK=100\n"\
"N_WRITE=100\n"\
"COMMAND=\"./adsrun --verbosity ${VERBOSITY} --theory ${THEORY} --latType ${LATTICE} --couplingType ${COUPLING_TYPE} --J ${J} --h ${H} --Lt ${TIMESLICES} --S1 ${CIRCUMFERENCE} --nTherm ${N_THERM} --nMeas ${N_MEAS} --nWrite ${N_WRITE} --nJkBlock ${N_JKBLK} --nSkip ${N_SKIP} --sigma ${SIGMA} --clusterAlg ${CLUSTER} --clusterArch ${CLUSTER_ARCH} \"\n"\
"echo \"\"\n"\
"echo \"Command given:\"\n"\
"echo ${COMMAND}\n"\
"echo \"\"\n"\
"${COMMAND}\n")
