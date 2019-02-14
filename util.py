import numpy as np

# Critical temperature for infinite volume system
J0 = 0.4406867935

Jlist = np.linspace(0.8, 1.2, 21)*J0

# Side lenght of square box
sideList = [8,16,32,64,128]

# Directory name based on value of J and Box size
def dirname(J, side):
    return "J={:.8f}_C={}".format(J,circ)

