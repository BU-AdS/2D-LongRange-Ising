from matplotlib import pyplot as plt
import numpy as np
from numpy import array, log, pi, exp, sqrt, e
from scipy.optimize import curve_fit
from sys import exit, argv
from util import *
import os
import sys
from paramplots import *

form = "png"

fname = "correlators_dth0.dat"
fname_mag = "JKobservables.dat"
miny = 1

if len(sys.argv)<2:
    print("{} <C>".format(argv[0]))
    exit(0)

C = int(argv[1])


for J in Jlist[1::2]:
    dname = "J={:.8f}_C={}_v2".format(J,C)
    os.chdir(dname)


    with open(fname, 'r') as f:
        corr = array([float(x.split()[1]) for x in f.readlines()])
        x = range(len(corr))

    with open(fname_mag, 'r') as f:
        # Read last line from file, corresponding to the most accurate estimator
        for line in f:
            pass
        last = line
        # Magnetization is third to last object according to data_io.cpp
        # FIXME
        mag = float(last.split()[-2])

    y = corr-mag**2

    print("J={}".format(J))
    print("corr={}".format(corr))
    print("y={}".format(y))

    c = next(color_cycle)
    # plt.loglog(x, y, linestyle='', marker='o', label=r'J={:.8f}'.format(J), c=c)
    # plt.loglog(x, corr, c=c)
    plt.plot(x, log(y), linestyle='', marker='o', label=r'J={:.8f}'.format(J), c=c)
    # miny = min(miny, min(y))
    miny = min(miny, min(log(y)))

    os.chdir('..')


plt.xlabel(r"r")
plt.ylabel(r"$\langle s(r) s(0)\rangle -\langle s\rangle^2$")
plt.xlim(min(x)-0.1,max(x)+0.1)
plt.ylim(miny-0.1, 1+0.1)
# plt.ylim(0.1, 1+0.1)
plt.legend(loc=3)
plt.title("C={}".format(C))
plt.savefig("corr_C={}.{}".format(C,form))
