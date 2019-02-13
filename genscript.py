import sys
import numpy as np
import os
import subprocess

execute = True

maxhours = 120
maxmem = 125

# Critical temperature for infinite volume system
J0 = 0.4406867935

Jlist = np.linspace(0.9, 1.1, 11)*J0
Jlist = np.linspace(1, 1, 1)*J0

for J in Jlist:
    dname = "J={}".format(J)
    fname = "template_ising.sh"

    scriptname = "{}/script.sh".format(dname)
    f = open(scriptname,'wt')

    jobname = "J={}".format(J)

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
