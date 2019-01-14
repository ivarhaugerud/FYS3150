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

T1_order = get_data("../final_data/c_order_1.000_expect.txt", ["MCs", "E", "EE", "M", "MM", "absM"])
T1_unorder = get_data("../final_data/c_unorder_1.000_expect.txt", ["MCs", "E", "EE", "M", "MM", "absM"])
T24_order = get_data("../final_data/c_order_2.400_expect.txt", ["MCs", "E", "EE", "M", "MM", "absM"])
T24_unorder = get_data("../final_data/c_unorder_2.400_expect.txt", ["MCs", "E", "EE", "M", "MM", "absM"])


plt.plot(np.array(T1_order["MCs"]), np.array(T1_order["E"]), label="$T=1.0J/k_B$ ordered")
plt.plot(np.array(T1_unorder["MCs"]), np.array(T1_unorder["E"]), label="$T=1.0J/k_B$ unordered")

plt.plot(np.array(T24_order["MCs"]), np.array(T24_order["E"]), label="$T=2.4J/k_B$ ordered")
plt.plot(np.array(T24_unorder["MCs"]), np.array(T24_unorder["E"]), label="$T=2.4J/k_B$ unordered")

plt.xscale("log")
plt.legend(loc="best", fontsize=13)
plt.xlabel(r"MC cycles", fontsize=14)
plt.ylabel(r"$\langle E \rangle $ $[J]$", fontsize=14)
plt.savefig("../figures/problem_c_E.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/problem_c_E.pdf", "../figures/problem_c_E.pdf"))
plt.show()

###

plt.plot(np.array(T1_order["MCs"]), np.array(T1_order["absM"])    , label="$T=1J/k_B$ ordered")
plt.plot(np.array(T1_unorder["MCs"]), np.array(T1_unorder["absM"]), label="$T=1J/k_B$ unordered")

plt.plot(np.array(T24_order["MCs"]), np.array(T24_order["absM"]),   label="$T=2.4J/k_B$ ordered")
plt.plot(np.array(T24_unorder["MCs"]), np.array(T24_unorder["absM"]), label="$T=2.4J/k_B$ unordered")

plt.xscale("log")
plt.legend(loc="best", fontsize=13)
plt.xlabel(r"MC cycles", fontsize=14)
plt.ylabel(r"$\langle |M| \rangle $", fontsize=14)
plt.savefig("../figures/problem_c_M.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/problem_c_M.pdf", "../figures/problem_c_M.pdf"))
plt.show()


temps = np.arange(0.2, 4.8001, 0.1)
final_vales = np.zeros(len(temps))

for counter, T in enumerate(temps):
    data = get_data("../final_data/task_c_2_%4.3f_accepts.txt" % T, ["accepts", "total"])
    nr_accepted = np.array(data["accepts"])[-1]

    #nr_possible = abs(np.array(data["total"])[-1])
    nr_possible = 20*20*1e6
    final_vales[counter] = nr_accepted/nr_possible*100

plt.xlabel(r"Possible accepted states")
plt.ylabel(r"Accepted states $[\%]$")

plt.plot(temps, final_vales, "o")
plt.xlabel(r"$T$ [$J/k_B$]", fontsize=14)
plt.ylabel(r"Accepted flips [%]", fontsize=14)
plt.savefig("../figures/problem_c_2.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/problem_c_2.pdf", "../figures/problem_c_2.pdf"))

plt.show()
