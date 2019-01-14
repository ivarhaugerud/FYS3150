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


temp = np.arange(2.02, 2.49, 0.01)
temp2 = np.arange(2.2090, 2.3265, 0.0025)
combined = np.concatenate((temp, temp2))
combined = np.sort(combined)
L_vals = [40, 60, 80, 100, 120, 200]

average_E = np.zeros(len(combined))
average_EE = np.zeros(len(combined))
average_M = np.zeros(len(combined))
average_MM = np.zeros(len(combined))
average_Mabs = np.zeros(len(combined))

counter = 0

for L in L_vals:
    for T in combined:
        try:
            results = get_data("../fresh_data/%i_%5.4f_expect.txt" % (L, T), ["MCs", "E", "EE", "M", "MM", "absM"])
        except:
            results = get_data("../fresh_data/%i_%4.3f_expect.txt" % (L, T), ["MCs", "E", "EE", "M", "MM", "absM"])

        half_way_index = int(len(results["E"])*0)
        number_MCs = np.array(results["MCs"])[-1]

        average_E[counter]    = np.mean(np.array(results["E"])[half_way_index:])
        average_EE[counter]   = np.mean(np.array(results["EE"])[half_way_index:])
        average_M[counter]    = np.mean(np.array(results["M"])[half_way_index:])
        average_MM[counter]   = np.mean(np.array(results["MM"])[half_way_index:])
        average_Mabs[counter] = np.mean(np.array(results["absM"])[half_way_index:])

        counter += 1

    plt.figure(1)
    plt.xlabel(r"$T$ $[J/k_B]$", fontsize=14)
    plt.ylabel(r"$\chi$ $[J^{-1}]$", fontsize=14)
    plt.plot(combined, (average_MM-average_Mabs**2)/(combined)*L*L, "o", label = r"$L = %d$" % L)

    plt.figure(2)
    plt.xlabel(r"$T$ $[J/k_B]$", fontsize=14)
    plt.ylabel(r"$C_V$ $[k_B]$", fontsize=14)
    plt.plot(combined, (average_EE-average_E**2)/(combined*combined)*L*L, "o", label = r"$L = %d$" % L)

    plt.figure(3)
    plt.xlabel(r"$T$ $[J/k_B]$", fontsize=14)
    plt.ylabel(r"$\langle | M | \rangle$", fontsize=14)
    plt.plot(combined, average_Mabs, "o", label = r"$L = %d$" % L)

    plt.figure(4)
    plt.xlabel(r"$T$ $[J/k_B]$", fontsize=14)
    plt.ylabel(r"$E$ $[J]$", fontsize=14)
    plt.plot(combined, average_E, "o", label = r"$L = %d$" % L)

    counter = 0
name = ["suscept", "heatcap", "magne", "energy"]
for i in range(1, 5):
    plt.figure(i)
    plt.legend(loc="best", fontsize=13)
    plt.savefig("../figures/" + name[i-1] + ".pdf", bbox_inches = "tight")
    os.system('pdfcrop %s %s &> /dev/null &'%("../figures/" + name[i-1] + ".pdf", "../figures/" + name[i-1] + ".pdf"))
plt.show()
