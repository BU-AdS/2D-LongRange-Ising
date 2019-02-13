import sys
import os
import numpy as np
from shutil import copyfile

execute = False

# Critical temperature for infinite volume system
J0 = 0.4406867935

Jlist = np.linspace(0.9, 1.1, 11)*J0

print(Jlist)

for J in Jlist:
    dname = "J={}".format(J)

    if not os.path.exists(dname):
        os.makedirs(dname)

    execname = "{}/adsrun".format(dname)
    copyfile("adsrun", execname)
    os.chmod(execname, 0o755)

    fname = "{}/template_ising.sh".format(dname)
    sys.stdout = open(fname,'wt')

    print(
    "#!/bin/bash\n"\
    "export OMP_NUM_THREADS=4\n"\
    "VERBOSITY='q'\n"
    "LATTICE='2D'\n"\
    "THEORY='ISING'\n"\
    "CLUSTER='WOLFF'\n"\
    "CLUSTER_ARCH='CPU'\n"\
    "COUPLING_TYPE='SR'\n"+\
    "J={}\n".format(J)+\
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

    os.chmod(fname, 0o755)


    sys.stdout.close()
