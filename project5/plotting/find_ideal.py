import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os

plt.style.use("bmh")
sns.color_palette("hls", 1)

import matplotlib
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

def get_data(filename, variables):
    df = pd.read_csv(filename,\
                    delim_whitespace=True, \
                    engine='python', \
                    names=variables)
    return df
    #using pandas to read the data files

alpha_min = 0.78
alpha_max = 1.415
alpha_step = 0.005
alpha_values = np.arange(alpha_min, alpha_max, alpha_step)

beta_min = 0
beta_max = 0.9
beta_step = 0.005
beta_values = np.arange(beta_min, beta_max, beta_step)

values = np.zeros( (len(beta_values), len(alpha_values)) )

for counter_alpha, alpha in enumerate(alpha_values):
    print(alpha)
    for counter_beta, beta in enumerate(beta_values):
        data = get_data("../big_data/parallel_expectationvalues_%4.3f_%4.3f.txt" % (alpha, beta), ["E", "EE", "r12"])

        values[counter_beta, counter_alpha] = np.mean(np.array(data["E"]))



import scipy.misc as scm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(alpha_values, beta_values)

plt.xlabel(r"$\beta$")
plt.ylabel(r"$\alpha$")
ax.set_zlabel(r"$E$")
wframe = None

current_lowest = 50
alpha_index = 0
beta_index = 0
for i in range(len(alpha_values)):
    for j in range(len(beta_values)):
        if values[j, i] < current_lowest:
            current_lowest = values[j, i]
            alpha_index = i
            beta_index = j

print(current_lowest, alpha_values[alpha_index], beta_values[beta_index])


#rotate the figure in 3D
for phi in np.linspace(-50, -25, 200):
    if wframe:
        ax.collections.remove(wframe)
    wframe = ax.plot_surface(X, Y, values, rstride=2, cstride=2, cmap=cm.coolwarm)
    ax.view_init(30, phi)
    plt.pause(.0001)

plt.savefig("../figures/3D.pdf")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/3D.pdf", "../figures/3D.pdf"))
plt.show()
