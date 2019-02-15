import matplotlib.pyplot as plt
from matplotlib import rc
from cycler import cycler
import matplotlib as mpl
from itertools import cycle

paramd = {'lines.linestyle': ['dashed', 'dashdot', 'solid', 'dotted'],
        'lines.marker': ['o','v','.','x'],
        'lines.markersize': [1,4,1,5]}

# print(mpl.rcParams.keys())
plt.style.use('ggplot')

# rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


color_cycle = cycle(['g', 'b', 'c', 'm', 'y', 'k'])


# plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
    # cycler('linestyle', ['-', '--', ':', '-.'])))


params = {'legend.fontsize': 8, 'figure.autolayout': True}
plt.rcParams.update(params)


def setparams(idx):
    for k,v in paramd.items():
        mpl.rcParams[k] = v[idx]
