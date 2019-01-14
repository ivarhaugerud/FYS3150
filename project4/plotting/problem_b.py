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

def linear_regresion(x, y):
    n = float(len(x))
    D = float(np.sum(np.square(x)) - (np.sum(x)**2)/n)
    E = float(np.sum(x*y) - np.sum(x)*np.sum(y)/n)
    F = float(np.sum(np.square(y)) - (np.sum(y)**2)/n)

    delta_m = np.sqrt((1/(n-2))*(D*F-E**2)/(D**2))
    delta_c = np.sqrt(1/(n-2)*(D/n+np.mean(x)**2)*(D*F-E**2)/(D**2))
    m = E/D
    c = np.mean(y)-m*np.mean(x)

    return m, c, delta_m, delta_c
    #using linear regression from Squires, with uncertainty to find slope and constant term

def expec_value_E(T):
    J = 1
    beta = 1/T
    Z = 4*(3+np.cosh(8*beta*J))
    return -32*J*np.sinh(8*beta*J)/Z

def expec_value_EE(T):
    J = 1
    beta = 1/T
    Z = 4*(3+np.cosh(8*beta*J))
    return 4*64*J*J*np.cosh(8*beta*J)/Z

def expec_value_M(T):
    J = 1
    beta = 1/T
    Z = 4*(3+np.cosh(8*beta*J))
    return 8*(2 + np.exp(8*beta*J))/Z

def expec_value_MM(T):
    J = 1
    beta = 1/T
    Z = 4*(3+np.cosh(8*beta*J))
    return 32*(1 + np.exp(beta*8*J))/Z # I actually find 32, oh well, this is per particle at least

def heat_cap(T):
    J = 1
    k = 1
    beta = 1/T
    return 1/(k*T*T)*(expec_value_EE(T)-expec_value_E(T)**2)


def sucept(T):
    J = 1
    k = 1
    beta = 1/T
    return 1/(k*T)*(expec_value_MM(T)-expec_value_M(T)**2)

temp = np.arange(0.5, 4.00, 0.1)

lent = len(temp)
average_E = np.zeros((lent, 2))
average_M = np.zeros((lent, 2))
average_EE = np.zeros((lent, 2))
average_MM = np.zeros((lent, 2))

counter = 0

for T in temp:
    results = get_data("../final_data/task_b_%4.3f_expect.txt" % (T), ["E", "EE", "M", "MM", "absM"])
    halfway_index = int(0.5*len(np.array(results["E"])))

    average_E[counter, 0] = 4*np.mean(np.array(results["E"])[halfway_index:])
    average_M[counter, 0] = 4*np.mean(np.array(results["absM"])[halfway_index:])
    average_EE[counter, 0] = 4*4*np.mean(np.array(results["EE"])[halfway_index:])
    average_MM[counter, 0] = 4*4*np.mean(np.array(results["MM"])[halfway_index:])

    average_E[counter, 1] = np.std(np.array(results["E"])[halfway_index:])
    average_M[counter, 1] = np.std(np.array(results["absM"])[halfway_index:])
    average_EE[counter, 1] = np.std(np.array(results["EE"])[halfway_index:])
    average_MM[counter, 1] = np.std(np.array(results["MM"])[halfway_index:])

    counter += 1

plt.figure(num=None, figsize=(10, 5), dpi=80, facecolor='w', edgecolor='k')

plt.subplot(2,3,1)
plt.errorbar(temp, average_E[:, 0]/4, fmt='o', yerr=average_E[:, 1], label="numerical E", zorder=1)
plt.plot(temp, expec_value_E(temp)/4, label="analytical E", zorder=2)
plt.ylabel(r"$ \langle E \rangle $ $[J]$", fontsize=14)

plt.subplot(2,3,2)
plt.errorbar(temp, average_EE[:, 0]/16, fmt='o', yerr=average_EE[:, 1], label="numerical M", zorder=1)
plt.plot(temp, expec_value_EE(temp)/16, label="analytical", zorder=2)
plt.ylabel(r"$\langle E^2 \rangle $ $[J^2]$", fontsize=14)

plt.subplot(2,3,4)
plt.errorbar(temp, average_M[:, 0]/4, fmt='o', yerr=average_M[:, 1], label="numerical E", zorder=1)
plt.plot(temp, expec_value_M(temp)/4, label="analytical E", zorder=2)
plt.ylabel(r"$\langle |M| \rangle $ ", fontsize=14)

plt.subplot(2,3,5)
plt.errorbar(temp, average_MM[:, 0]/16, fmt='o', yerr=average_MM[:, 1], label="numerical M", zorder=1)
plt.plot(temp, expec_value_MM(temp)/16, label="analytical", zorder=2)
plt.ylabel(r"$\langle M^2 \rangle $", fontsize=14)

plt.subplot(2,3,3)
plt.plot(temp, (average_EE[:, 0]-average_E[:, 0]**2)/(4*temp*temp), "o", zorder=1)
plt.plot(temp, heat_cap(temp)/4, zorder=2)
plt.ylabel(r"$C_V$ [$k_B$]", fontsize=14)

plt.subplot(2,3,6)
plt.plot(temp, (average_MM[:, 0]-average_M[:, 0]**2)/(4*temp), "o")
plt.plot(temp, sucept(temp)/4)
plt.ylabel(r"$\chi$ $[J^{-1}]$", fontsize=14)

plt.subplot(234)
plt.xlabel(r"$T$ $[k_B/J]$", fontsize=14)
plt.subplot(235)
plt.xlabel(r"$T$ $[k_B/J]$", fontsize=14)
plt.subplot(236)
plt.xlabel(r"$T$ $[k_B/J]$", fontsize=14)

plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
#plt.savefig("../figures/taskb.pdf", bbox_inches="tight")
plt.show()
#os.system('pdfcrop %s %s &> /dev/null &'%("../figures/taskb.pdf", "../figures/taskb.pdf"))
