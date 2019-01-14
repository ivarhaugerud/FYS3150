import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os

#PLOT THE 3D data for alpha-vs-beta graph

def get_data(filename, variables):
    df = pd.read_csv(filename,\
                    delim_whitespace=True, \
                    engine='python', \
                    names=variables)
    return df
    #using pandas to read the data files

amin = 0.78
amax = 1.415
delta = 0.005

alpha = np.arange(amin, amax, delta)

data = {} #loop over all the text file
for i in range(len(alpha)):
    data["%4.3f" % alpha[i]] = get_data("../data/parallel_%4.3f.txt" % alpha[i], ["b", "E", "EE"])

beta = np.array(data["%4.3f" % alpha[0]]["b"])

E = np.zeros((len(alpha), len(beta)))
EE = np.zeros((len(alpha), len(beta)))
var = np.zeros((len(alpha), len(beta)))

for i in range(len(alpha)):
    current = data["%4.3f" % alpha[i]]
    E[i, :] = np.array(current["E"])
    EE[i, :] = np.array(current["EE"])
    var[i, :] = np.array(current["EE"]) - np.array(current["E"])**2


current_lowest = 50
alpha_index = 0
beta_index = 0
for j in range(len(alpha)):
    for i in range(len(beta)):
        if E[j, i] < current_lowest:
            current_lowest = E[j, i]
            alpha_index = j
            beta_index = i

print(current_lowest, (EE[alpha_index, beta_index]-E[alpha_index, beta_index]**2), alpha[alpha_index], beta[beta_index])



import scipy.misc as scm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import matplotlib
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(beta, alpha)

plt.xlabel(r"$\beta$", fontsize=14)
plt.ylabel(r"$\alpha$", fontsize=14)
ax.set_zlabel(r"$E$ [a.u.]", fontsize=14)
wframe = None


#rotate the figure in 3D
for phi in np.linspace(-50, -25, 10):
    if wframe:
        ax.collections.remove(wframe)
    wframe = ax.plot_surface(X, Y, E, rstride=2, cstride=2, cmap=cm.coolwarm)
    ax.view_init(30, phi)
    plt.pause(.0001)

plt.savefig("../figures/3D.pdf")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/3D.pdf", "../figures/3D.pdf"))
plt.show()
