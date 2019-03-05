import numpy as np
from numpy import array, log, exp, sqrt

# Critical temperature for infinite volume system
J0 = 0.4406867935


# Jlist = np.linspace(0.8, 1.2, 11)*J0
Jlist = np.linspace(0.9, 1.1, 21)*J0

# Side lenght of square box
sideList = [16,32,64,128,256]

# Value of critical J for given box size, extracted from linear fit
def Jcrit(L):
    c = 0.402989
    return log(1+sqrt(2))/2-c/L

JcritDict = {L: array([0.95, 1.00, 1.05])*Jcrit(L) for L in sideList}

# Directory name based on value of J and Box size
def dirname(J, C):
    return "J={:.8f}_C={}_v4".format(J,C)
    # return "J={:.8f}_C={}_v5".format(J,C)

