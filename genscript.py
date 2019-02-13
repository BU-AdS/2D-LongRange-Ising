import sys

maxhours = 120
maxmem = 125

# Critical temperature for infinite volume system
J0 = 0.4406867935

Jlist = np.linspace(0.9, 1.1, 11)*J0

for J in Jlist:
    dname = "J={}".format(J)
    fname = "{}/template_ising.sh".format(dname)
    sys.stdout = open(fname,'wt')

 

fname = "script.sh"


jobname = "k={},L={},ET={}".format(k,L,ET)

f = open(fname, "w")

f.write("#!/bin/bash -l\n\n")
f.write("#$ -l h_rt={}:00:00\n".format(maxhours))
f.write("#$ -l mem_total={}G\n".format(maxmem))
f.write("#$ -m a\n")
f.write("#$ -j y\n")
#f.write("#$ -pe omp 2\n")
f.write("#$ -N {}\n\n".format(jobname))

f.write("module load anaconda3\n")
f.write("source activate py3\n\n")

f.write("python Phi4/eigsvsG.py {} {} {}\n\n".format(k, L, ET))

f.close()



