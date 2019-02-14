from matplotlib import pyplot as plt
import numpy as np
from numpy import array, log, pi, exp, sqrt, e
from scipy.optimize import curve_fit
from sys import exit

fname = "observables.dat"

with open(fname, 'r') as f:
    mag = array([list(map(float, x.split())) for x in f.readlines()])

t = mag[:,0]
avemag = np.cumsum(mag[:,1]) / np.array(range(1,t.size+1))

print("Average magnetization: {}".format(avemag[-1]))

plt.plot(t, avemag, linestyle='--', marker='o')

plt.xlabel("t")
plt.ylabel("<s>")
plt.savefig("mag.pdf")
