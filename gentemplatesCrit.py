import sys
import os
import numpy as np
from shutil import copyfile
from util import *

execute = False

nthreads = 1

print(Jlist)

for L in sideList:
    for J in JcritDict[L]:
        dname = dirname(J, L)

        if not os.path.exists(dname):
            os.makedirs(dname)

        execname = "long_range"
        copyfile(execname, "{}/{}".format(dname,execname))
        os.chmod("{}/{}".format(dname,execname), 0o755)

        fname = "{}/template_ising.sh".format(dname)
        sys.stdout = open(fname,'wt')

        print(
        "#!/bin/bash\n"\
        "export OMP_NUM_THREADS="+str(nthreads)+"\n"\
        "VERBOSITY='q'\n"
        "THEORY='ISING'\n"\
        "CLUSTER='WOLFF'\n"\
        "CLUSTER_ARCH='CPU'\n"\
        "COUPLING_TYPE='SR'\n"+\
        "J={}\n".format(J)+\
        "H=0.0\n"\
        "SIGMA=1.5\n"\
        "CIRCUMFERENCE={}\n".format(L)+\
        "TIMESLICES={}\n".format(L)+\
        "N_THERM=500\n"\
        "N_MEAS=10000\n"\
        "N_SKIP=50\n"\
        "N_JKBLK=100\n"\
        "N_WRITE=100\n"\
        "COMMAND=\"./"+execname+\
        " --verbosity ${VERBOSITY} --theory ${THEORY} --couplingType ${COUPLING_TYPE} --J ${J} --h ${H} --Lt ${TIMESLICES} --S1 ${CIRCUMFERENCE} --nTherm ${N_THERM} --nMeas ${N_MEAS} --nWrite ${N_WRITE} --nJkBlock ${N_JKBLK} --nSkip ${N_SKIP} --sigma ${SIGMA} --clusterAlg ${CLUSTER} --clusterArch ${CLUSTER_ARCH} \"\n"\
        "echo \"\"\n"\
        "echo \"Command given:\"\n"\
        "echo ${COMMAND}\n"\
        "echo \"\"\n"\
        "${COMMAND}\n")

        os.chmod(fname, 0o755)


    sys.stdout.close()
