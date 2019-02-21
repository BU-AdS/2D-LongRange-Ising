from matplotlib import pyplot as plt
import numpy as np
from numpy import array, log, pi, exp, sqrt, e
from scipy.optimize import curve_fit
from sys import exit, argv
from paramplots import *
from util import *
import sys, os
import pandas as pd

form = "png"

fname = "JKobservables.csv"
miny = 1

if len(sys.argv)<2:
    print("{} <C>".format(argv[0]))
    exit(0)

C = int(argv[1])


for J in Jlist[1::2]:
    dname = "J={:.8f}_C={}_v3".format(J,C)
    os.chdir(dname)

    data = pd.read_csv(fname)

    # with open(fname, 'r') as f:
        # t = array([int(line.split()[0]) for line in f.readlines()])
    # with open(fname, 'r') as f:
        # mag = array([float(line.split()[-2]) for line in f.readlines()])

    t = data['idx'].tolist()
    mag = data['avePhi'].tolist()

    print(t)
    print(mag)

    os.chdir('..')

    avemag = np.cumsum(mag) / np.array(range(1,len(mag)+1))

    print("J={}, Average magnetization: {}".format(J, avemag[-1]))

    plt.plot(t, avemag, linestyle='--', marker='o')

    plt.xlabel(r"$t$")
    plt.ylabel(r"$\langle s\rangle$")
    plt.savefig("mag_J={}_C={}.png".format(J,C))
    plt.clf()
