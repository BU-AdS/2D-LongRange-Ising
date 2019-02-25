import numpy as np

# Critical temperature for infinite volume system
J0 = 0.4406867935


# Jlist = np.linspace(0.8, 1.2, 11)*J0
Jlist = np.linspace(0.9, 1.1, 21)*J0

# Side lenght of square box
sideList = [16,32,64,128,256]

# Directory name based on value of J and Box size
def dirname(J, C):
    return "J={:.8f}_C={}_v4".format(J,C)

