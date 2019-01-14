import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

#PLOT THE 3D data for rho-vs-n graph

def get_data(filename, variables):
    df = pd.read_csv(filename,\
                    delim_whitespace=True, \
                    engine='python', \
                    names=variables)
    return df
    #using pandas to read the data files

q_min = 2.5 #lowest  value of rho_max
q_max = 7.5 #largest value of rho_max
n_min = 100 #minimum n value
n_max = 200 #maximum n value
number_of_q = 20                #number of iterations
q = np.zeros(number_of_q + 6)   #add some extra elements
q[:-6] = np.linspace(q_min, q_max, number_of_q)

#add some extra values
q[-6] = 4.45
q[-5] = 4.50
q[-4] = 4.55
q[-3] = 4.60
q[-2] = 4.65
q[-1] = 4.75

index_change = np.argsort(q)

n = np.linspace(n_min, n_max, number_of_q+6) #the n values we will use


data = {} #loop over all the text file
for i in range(len(n)):
    data["%4.3f" % n[i]] = get_data("data/data_%4.3f.txt" % n[i], ["q", "eig_0", "eig_1", "eig_2", "eig_3", "eig_4"])

Z = np.zeros((number_of_q+6, number_of_q+6))

for i in range(len(n)):
    current = data["%4.3f" % n[i]] #calculate the average deviation of the three lowest energy states eigenvalues
    Z[i, :] = np.log10((current["eig_0"] + current["eig_1"] + current["eig_2"])/3)

print(q[np.argmin(Z[-1, :])], "Ideal value of rho_max")


import scipy.misc as scm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(q, n)

plt.ylabel(r"matrix size $ N$")
plt.xlabel(r"maximal position $\rho_{max}$")
ax.set_zlabel(r"$\log10{(\langle\lambda_{num}-\lambda_{ana}\rangle)}$")
wframe = None
#rotate the figure in 3D
for phi in np.linspace(0, 70, 200):
    if wframe:
        ax.collections.remove(wframe)
    wframe = ax.plot_surface(X[:, index_change], Y, Z[:, index_change], rstride=2, cstride=2, cmap=cm.coolwarm)
    ax.view_init(30, phi)
    plt.pause(.0001)

plt.savefig("../figures/rho_vs_n.pdf")
plt.show()
