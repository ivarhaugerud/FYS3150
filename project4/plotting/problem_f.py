import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.signal as scs
import os

plt.style.use("bmh")
sns.color_palette("hls", 1)

import matplotlib
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

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


def find_TC(T_C, L):
    nu = 1.0

    delta_Tc = T_C[-1]-T_C[-2]
    delta_L = 1.0/(L[-1])-1/(L[-2])
    a = delta_Tc/delta_L
    return T_C - a/L

def get_data(filename, variables):
    df = pd.read_csv(filename,\
                    delim_whitespace=True, \
                    engine='python', \
                    names=variables)
    return df

Ls = np.array([40, 60, 80, 100, 120, 200])
Linvs = np.array([1.0/40, 1.0/60, 1.0/80, 1.0/100, 1.0/120, 1.0/200])
Ts = np.arange(2.2090, 2.3265, 0.0025)
T_c = 2.269185

TC_s = np.zeros(len(Ls))
average_M = np.zeros(len(Ts))
average_MM = np.zeros(len(Ts))


for L_i, L in enumerate(Ls):
    for T_i, T in enumerate(Ts):
        results = get_data("../fresh_data/%i_%5.4f_expect.txt" % (L, T), ["MCs", "E", "EE", "M", "MM", "absM"])

        first_5_percent = int(0.5*len(results["E"]))

        average_M[T_i]  = np.mean(np.array(results["E"])[first_5_percent:])
        average_MM[T_i] = np.mean(np.array(results["EE"])[first_5_percent:])

    suscept = (average_MM-average_M**2)/(Ts*Ts)*L*L
    filtered = scs.savgol_filter(suscept, 47, 7)

    plt.plot(Ts, suscept, "o", c=sns.color_palette()[L_i], alpha = 0.5)
    plt.plot(Ts, filtered, "-", c=sns.color_palette()[L_i], label = r"$L = %d$" % L)
    TC_s[L_i] = Ts[np.argmax(filtered)]

plt.xlabel(r"$T$ $[J/k_B]$", fontsize=14)
plt.ylabel(r"$C_V$ $[k_B]$", fontsize=14)
plt.legend(loc = 'best', fontsize=13)
plt.savefig("../figures/fitted_suscept.pdf")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/fitted_suscept.pdf", "../figures/fitted_suscept.pdf"))
plt.show()



print(TC_s)
print(find_TC(TC_s, Ls))
print("Relative error: ", (find_TC(TC_s, Ls) - T_c)/T_c)

[m, c, dm, dc] = linear_regresion(Linvs, TC_s)

#plt.plot(Linvs, np.ones(len(Ls))*T_c, "--")
plt.plot(Linvs, TC_s, "o")
plt.plot(np.linspace(0, Linvs[0], 100), m*np.linspace(0, Linvs[0], 100) + c)
plt.fill_between(np.linspace(0, Linvs[0], 100), (m+dm)*np.linspace(0, Linvs[0], 100) + c + dc, (m-dm)*np.linspace(0, Linvs[0], 100) + c -dc, alpha=0.3, color=sns.color_palette()[1])
plt.plot(0, T_c, '*')
plt.xlabel(r"$1/L$", fontsize=14)
plt.ylabel(r'$T$ $[J/k_B]$', fontsize=14)
plt.savefig("../figures/TC.pdf")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/TC.pdf", "../figures/TC.pdf"))
plt.show()

print(c, " +- ", dc)
