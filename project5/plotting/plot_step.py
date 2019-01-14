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

alpha_min = 0.3
alpha_max = 2.5
alpha_step = 0.1
alpha_values = np.arange(alpha_min, alpha_max, alpha_step)

best_step = np.zeros(len(alpha_values))
counter = 0

for counter, i in enumerate(alpha_values):
    data = get_data("../data/find_step_accepts_%3.2f_1.00.txt" % i, ["step", "accepts"])

    if counter % 3 == 0:
        plt.plot(data["step"], data["accepts"]*100, "-", label=r"$\alpha=$%3.2f" % i)

    best_step_argument = np.argmin( abs( data["accepts"]-0.5 ))
    best_step[counter] = np.array(data["step"])[best_step_argument]

plt.legend(loc="best", fontsize=13, ncol=2)
plt.plot(data["step"], np.ones(len(data["step"]))*50, "--", color="k")
plt.xlabel(r"Step length $\delta$ ", fontsize=14)
plt.ylabel(r"Percentage of accepts [%]", fontsize=14)

plt.savefig("../figures/acceptence_rate.pdf", bbox_inches = "tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/acceptence_rate.pdf", "../figures/acceptence_rate.pdf"))
plt.show()

m, c, delta_m, delta_c = linear_regresion(np.log10(alpha_values), np.log10(best_step))

plt.plot(alpha_values, 10**c*10**(m*np.log10(alpha_values)))
plt.plot(alpha_values, best_step, "o")

plt.xlabel(r"Variational parameter $\alpha$", fontsize=14)
plt.ylabel(r"Ideal step length $\delta$", fontsize=14)
plt.legend(loc="best", fontsize=13)
plt.xscale("log")
plt.yscale("log")

plt.savefig("../figures/ideal_step.pdf", bbox_inches = "tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/ideal_step.pdf", "../figures/ideal_step.pdf"))
plt.show()

print("m: ", m)
print("c: ", c)

print("m: ", delta_m)
print("c: ", delta_c)
