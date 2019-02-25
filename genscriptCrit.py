import sys
import numpy as np
import os
import subprocess
from util import *

execute = True

maxhours = 120
maxmem = 125

for L in sideList:
    for J in JcritDict[L]:
        dname = dirname(J, L)
        fname = "template_ising.sh"

        scriptname = "{}/script.sh".format(dname)
        f = open(scriptname,'wt')

        jobname = dname

        f.write("#!/bin/bash -l\n\n")
        f.write("#$ -l h_rt={}:00:00\n".format(maxhours))
        f.write("#$ -l mem_total={}G\n".format(maxmem))
        f.write("#$ -m a\n")
        f.write("#$ -j y\n")
        #f.write("#$ -pe omp 2\n")
        f.write("#$ -N {}\n\n".format(jobname))

        f.write("module load gsl\n")
        f.write("module load gcc/7.2.0\n")

        # f.write("{}/{}\n\n".format(dname,fname))
        f.write("./{} 1>/dev/null\n\n".format(fname))
        f.close()

        os.chmod(scriptname, 0o755)

        if execute:
                os.chdir(dname)
                subprocess.run(["qsub", "script.sh"])
                os.chdir('..')
