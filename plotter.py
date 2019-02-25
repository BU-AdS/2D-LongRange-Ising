from matplotlib import pyplot as plt
import numpy as np
from numpy import array, log, pi, exp, sqrt, e
from scipy.optimize import curve_fit
from sys import exit, argv
from util import *
import os
import sys
import pandas as pd
from paramplots import *

form = "png"

fname = "correlators_dth0.dat"
fname_mag = "JKobservables.csv"
miny = 1

# Index of J at which we perform the fit
# ifit = 1
Jfit = {64:0.436279925565, 128:0.44068, 256:0.44068}
Jrange = {64:[0.43, 0.441], 128:[0.435,0.448], 256:[0.435,0.448]}

if len(sys.argv)<2:
    print("{} <C>".format(argv[0]))
    exit(0)

C = int(argv[1])


for J in [x for x in Jlist if x>Jrange[C][0] and x<Jrange[C][1]]:
    dname = dirname(J,C)
    os.chdir(dname)

    with open(fname, 'r') as f:
        corr = array([float(x.split()[1]) for x in f.readlines()])
        x = range(len(corr))

    data = pd.read_csv(fname_mag)
    # Time series of total magnetization
    magseries = data['avePhi'].tolist()
    # magseries = data['avePhiAb'].tolist()
    # Integrated time series of total magnetization
    intmagseries = np.cumsum(magseries) / np.array(range(1,len(magseries)+1))
    # Statistical average of the total magnetization over the time series.
    # XXX To get rid of initial condition we should probably discard the initial data
    mag = intmagseries[-1]

    # Subtract disconnected component
    y = corr-mag**2

    print("J={}".format(J))
    print("corr={}".format(corr))
    print("y={}".format(y))

    c = next(color_cycle)
    plt.loglog(x, y, linestyle='', marker='o', label=r'J={:.8f}'.format(J), c=c)
    plt.loglog(x, corr, c=c)
    miny = min(miny, min(y))

    if abs(J-Jfit[C])<10**-5:
        print("J={} selected for fit".format(J))
        xx = range(1,10)
        y = [corr[xi]-mag**2 for xi in xx]
        xx = array(xx)
        coef = np.polyfit(log(xx), log(y), 1)
        print("Power of the fit: {:.8f}".format(coef[0]))
        xfit = np.linspace(0,max(log(x)),100)
        plt.loglog(exp(xfit), exp(coef[1]+coef[0]*xfit), label="fit")
        #, label={r"$C*x^{{{:.4f}}}$".format(coef[0])})
        # plt.text(x=2, y=1, s=r"$x^{{{:.4f}}}$".format(coef[0]), fontsize=14)

    # plt.plot(x, log(y), linestyle='', marker='o', label=r'J={:.8f}'.format(J), c=c)
    # miny = min(miny, min(log(y)))

    os.chdir('..')


plt.xlabel(r"r")
plt.ylabel(r"$\langle s(r) s(0)\rangle -\langle s\rangle^2$")
plt.xlim(min(x)-0.1,max(x)+0.1)
# plt.ylim(miny-0.1, 1+0.1)
# plt.ylim(0.1, 1+0.1)
plt.legend(loc=3)
plt.title("C={}, pow={:.4f}".format(C, coef[0]))
plt.savefig("corr_C={}.{}".format(C,form))
