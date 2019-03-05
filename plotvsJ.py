from matplotlib import pyplot as plt
import numpy as np
from numpy import array, log, pi, exp, sqrt, e
from scipy.optimize import curve_fit
from sys import exit
import os
from paramplots import *
from util import *
import sys, os
import pandas as pd
import copy


form = "png"
fname = "JKobservables.csv"

color_cycle2 = copy.deepcopy(color_cycle)

for C in sideList:

    binder = []
# The susceptibility here is defined as <M^2> - <|M|>^2
    suscep = []

    for J in Jlist:
        dname = dirname(J, C)
        os.chdir(dname)

        data = pd.read_csv(fname)

        binder.append(data["Binder"].tolist()[-1])

        suscep.append(data["Suscep"].tolist()[-1])

#         # Time series of the total magnetization
        # magTS = data["avePhi"]
        # magSqTS= magTS**2
        # magAbsTS = abs(magTS)

        # T = len(magTS)
        # # print("T = {}".format(T))
        # suscep.append(1/T*sum(magSqTS)-(1/T*sum(magAbsTS))**2)


#         # Read Binder cumulant
        # with open(fname, 'r') as f:

            # # Read last line from file, corresponding to the most accurate estimator
            # for line in f:
                # pass
            # last = line

            # # Binder cumulant is second to last object according to data_io.cpp
            # binder.append(float(last.split()[-2]))

        os.chdir('..')

    binder = np.array(binder)
    suscep = np.array(suscep)

    plt.figure(1)
    c = next(color_cycle)
    plt.plot(Jlist, binder, linestyle='--', marker='o', c=c,
            label="C={}".format(C))


    plt.figure(2)
    c2 = next(color_cycle2)
    plt.plot(Jlist, suscep, linestyle='--', marker='o', c=c2,
            label="C={}".format(C))



plt.figure(1)
plt.ylim(0.55,0.7)
plt.xlim(0.42,0.46)
plt.axvline(J0)
plt.xlabel("J")
plt.ylabel("Binder")
plt.legend()
plt.savefig("Binder.{}".format(form))

plt.figure(2)
plt.axvline(J0)
plt.xlabel("J")
plt.ylabel("Suscep")
plt.legend()
plt.savefig("Suscep.{}".format(form))
