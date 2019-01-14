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


results = get_data("../final_data/timer_timer.txt", ["L", "time_para", "time_non_para"])
L = np.array(results["L"])
m, c, delta_m, delta_c = linear_regresion( np.log10(L), np.log10(results["time_para"]))

print("SLOPE PARALELL: ", m)
print(delta_m)

plt.plot(L, results["time_para"], "o", label="Parallel")
plt.plot(L, 10**( c + m*np.log10(L)))

m, c, delta_m, delta_c = linear_regresion( np.log10(L), np.log10(results["time_non_para"]))
print("SLOPE NON-PARALLEL: ", m)
print(delta_m)

plt.plot(L, results["time_non_para"], "o", label="Non-parallel")
plt.plot(L, 10**( c + m*np.log10(L)))

plt.xlabel("$L$", fontsize=14)
plt.ylabel("Time used [s]", fontsize=14)
plt.legend(loc="best", fontsize=13)
plt.savefig("../figures/timer.pdf", bbox_inches="tight")
plt.show()
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/timer.pdf", "../figures/timer.pdf"))

plt.show()

mean = np.mean(results["time_non_para"]/results["time_para"])
std  =  np.std(results["time_non_para"]/results["time_para"])

print("mean:", mean)
print("std:",  std)
