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

alpha = 0.995
beta  = 0.28

data  = get_data("../data/virial_E2_0.995.txt", ["beta", "omega", "E", "EE", "r12", "pot", "1/r12"])
data2 = get_data("../data/virial_E0_1.000.txt", ["beta", "omega", "E", "EE", "r12", "pot", "1/r12"])

beta = np.array(data["beta"])
omega = np.array(data["omega"])
E = np.array(data["E"])
EE = np.array(data["EE"])
r12 = np.array(data["r12"])
pot = np.array(data["pot"])
kin = E - pot
r12inv = np.array(data["1/r12"])

beta2 = np.array(data2["beta"])
omega2 = np.array(data2["omega"])
E2 = np.array(data2["E"])
EE2 = np.array(data2["EE"])
r122 = np.array(data2["r12"])

pot2 = np.array(data2["pot"])
kin2 = E2 - pot2


pHO = (pot-r12inv)/pot
pC = r12inv/pot

plt.scatter(omega, pot/kin,    s=22, label="interacting")
plt.scatter(omega2, pot2/kin2, s=22, label="noninteracting")
#plt.plot(omega, (1*pHO-0.5*pC)**(-1), "o")

plt.legend(loc="best", fontsize=14)

plt.xlabel(r"Characteristic frequency $\omega$ ", fontsize=14)
plt.ylabel(r"$\langle V \rangle/\langle T\rangle$", fontsize=14)

plt.savefig("../figures/virial.pdf", bbox_inches = "tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/virial.pdf", "../figures/virial.pdf"))

plt.show()

plt.scatter(omega, (pot-r12inv)/kin,    s=22, label=r"$\langle V_{H0} \rangle /\langle T\rangle $")
plt.scatter(omega, r12inv/kin,        s=22, label=r"$\langle V_{C}  \rangle/\langle T\rangle$")

plt.legend(loc="best", fontsize=15)

plt.xlabel(r"Characteristic frequency $\omega$ ", fontsize=14)
plt.ylabel(r"Relative potential energy contribution", fontsize=14)
plt.axis([-0.01, 1.01, 0.4, 2.0])
plt.savefig("../figures/H0_vs_C.pdf", bbox_inches = "tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/H0_vs_C.pdf", "../figures/H0_vs_C.pdf"))

plt.show()
