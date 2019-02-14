from matplotlib import pyplot as plt
import numpy as np
from numpy import array, log, pi, exp, sqrt, e
from scipy.optimize import curve_fit
from sys import exit
import os

# Critical temperature for infinite volume system
J0 = 0.4406867935

fname = "JKobservables.dat"
Jlist = np.linspace(0.9, 1.1, 11)*J0

for J in Jlist:
    print("{:.8f}".format(J))

binder = []

for J in Jlist:
    dname = "J={:.8f}".format(J)
    os.chdir(dname)

    # Read Binder cumulant
    with open(fname, 'r') as f:

        # Read last line from file, corresponding to the most accurate estimator
        for line in f:
            pass
        last = line

        # Binder cumulant is second to last object according to data_io.cpp
        binder.append(float(last.split()[-2]))

    os.chdir('..')

binder = np.array(binder)

plt.plot(Jlist, binder, linestyle='--', marker='o', color='b' )
plt.axvline(J0)

plt.xlabel("J")
plt.ylabel("Binder")
plt.savefig("Binder.pdf")
