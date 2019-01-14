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


results = get_data("../data/task_d_1.000_energy.txt", ["E"])
results2 = get_data("../data/task_d_2.400_energy.txt", ["E"])

energy1 = np.array(results["E"])/(20*20)
energy2 = np.array(results2["E"])/(20*20)

first_15 = int(0.0*len(energy1))

plt.yscale("log")

print(np.mean(energy2[first_15:]), "mean energy T=2.4")
print(np.std(energy2[first_15:]),  "STD energy  T=2.4")

print(np.mean(energy1[first_15:]), "mean energy T=1")
print(np.std(energy1[first_15:]), "STD energy, T=1")

plt.hist(energy1[first_15:], bins=6, label=r"T=1.0", histtype='bar', ec='black', normed=True)
plt.hist(energy2[first_15:], bins=50, alpha=0.60, label=r"T=2.4", histtype='bar', ec='black', normed=True)

plt.legend(loc="best", fontsize=13)
plt.xlabel(r"$E$ $[J]$", fontsize=14)
plt.ylabel(r"Probability $[\%]$", fontsize=14)

#plt.savefig("../figures/problem_d.pdf", bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%("../figures/problem_d.pdf", "../figures/problem_d.pdf"))

plt.show()


results  = get_data("../final_data/c_unorder_1.000_expect.txt", ["MC", "E", "EE", "M", "MM", "absM"])
results2 = get_data("../final_data/c_unorder_2.400_expect.txt", ["MC", "E", "EE", "M", "MM", "absM"])

plt.plot(results["EE"]-results["E"]**2)
plt.plot(results2["EE"]-results2["E"]**2)
plt.show()
