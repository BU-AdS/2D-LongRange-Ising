from matplotlib import pyplot as plt
import numpy as np
from numpy import array, log, pi, exp, sqrt, e
from scipy.optimize import curve_fit
from sys import exit
from util import *
import os

def fitfun(x, a, b, c):
    return a*x**(-b)+c

fname = "correlators_dth0.dat"
miny = 1

for J in Jlist[::2]:
    dname = "J={:.8f}".format(J)
    os.chdir(dname)

    with open(fname, 'r') as f:
        corr = array([float(x.split()[1]) for x in f.readlines()])
        x = range(len(corr))

        plt.loglog(x, corr, linestyle='--', marker='o', label=dname)
        miny = min(miny, min(corr))

    os.chdir('..')

plt.xlabel("r")
plt.ylabel("<s(r) s(0)>")
plt.xlim(min(x)-0.1,max(x)+0.1)
plt.ylim(miny-0.1, 1+0.1)
plt.legend(loc=3)
plt.savefig("corr.pdf")
