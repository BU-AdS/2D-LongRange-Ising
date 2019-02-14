from matplotlib import pyplot as plt
import numpy as np
from numpy import array, log, pi, exp, sqrt, e
from scipy.optimize import curve_fit
from sys import exit

def fitfun(x, a, b, c):
    return a*x**(-b)+c

fname = "correlators_dth0.dat"

with open(fname, 'r') as f:
    corr = array([float(x.split()[1]) for x in f.readlines()])

print(corr)
n = 10

x = range(1,len(corr)+1)
plt.loglog(x, corr, linestyle='--', marker='o')


plt.xlabel("log r")
plt.ylabel("log <s(r) s(0)>")
plt.savefig("corr.pdf")

exit(0)

xlog = log(range(1,corr.size+1))
ylog = log(corr)
m,b = np.polyfit(xlog[:n], ylog[:n], 1)


popt, _ = curve_fit(fitfun, x, corr)
print(popt)

plt.loglog(x, corr-popt[2], linestyle='--', marker='o')

xs = np.linspace(min(x), max(x), 100)
plt.loglog(xs, fitfun(xs, *popt)-popt[2])

# plt.loglog(corr, basex=e, basey=e, linestyle='--', marker='o')

# plt.plot(x, y)
# plt.plot(xs, m*xs+b, label='fit')


