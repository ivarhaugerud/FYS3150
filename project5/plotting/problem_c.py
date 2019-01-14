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

def insort(a, b, kind='mergesort'):
    # took mergesort as it seemed a tiny bit faster for my sorted large array try.
    c = np.concatenate((a, b)) # we still need to do this unfortunatly.
    c.sort(kind=kind)
    flag = np.ones(len(c), dtype=bool)
    np.not_equal(c[1:], c[:-1], out=flag[1:])
    return c[flag]

for i in [0.6, 0.8, 1.0, 1.2, 1.4]:
    #data = get_data("../data/stability_expectationvalues_%3.2f_1.00.txt" % i, ["MCC", "E", "EE", "r12"])
    data = get_data("../data/final_data_2_expectationvalues_%4.3f_1.000.txt" % i , ["MCC", "E", "EE", "r12"])
    MCC = np.array(data["MCC"])
    E = np.array(data["E"])
    EE = np.array(data["EE"])
    r12 = np.array(data["r12"])

    plt.figure(3)
    plt.scatter(MCC, E, s=20.0, label=r"$\alpha = $%3.2f" % i)
    plt.xscale("log")

    plt.figure(4)
    plt.scatter(MCC, EE - E*E, s=20.0, label=r"$\alpha = $%3.2f" % i)
    plt.xscale("log")


plt.figure(3)
plt.legend(loc="best", fontsize=13)
plt.xlabel(r"Monte Carlo cycles", fontsize=14)
plt.ylabel(r"Expectation value energy $\langle E \rangle $ [a.u]", fontsize=14)
plt.savefig("../figures/stability_energy_E1.pdf", bbox_inches = "tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/stability_energy_E1.pdf", "../figures/stability_energy_E1.pdf"))

plt.figure(4)
plt.legend(loc="best", fontsize=13)
plt.xlabel(r"Monte Carlo cycles", fontsize=14)
plt.ylabel(r"Variance energy $σ_E^2$ [a.u$^2$]", fontsize=14)
plt.savefig("../figures/stability_variance_E1.pdf", bbox_inches = "tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/stability_variance_E1.pdf", "../figures/stability_variance_E1.pdf"))
plt.show()

alpha1 = np.arange(0.3, 2, 0.05)
alpha2 = np.arange(0.9, 1.1, 0.01)

alpha = insort(alpha1, alpha2)

E  = np.zeros(len(alpha))
EE = np.zeros(len(alpha))
r12 = np.zeros(len(alpha))

for counter, i in enumerate(alpha):
    data = get_data("../data/analytical_E0_%3.2f.txt" % i, ["MCC", "E", "EE", "r12"])

    E[counter]  = np.mean(np.array(data["E"]))
    EE[counter] = np.mean(np.array(data["EE"]))
    r12[counter] = np.mean(np.array(data["r12"]))


alpha1 = np.arange(0.3, 2, 0.05)
alpha2 = np.arange(0.7, 1.0, 0.01)
alpha3 = insort(alpha1, alpha2)

E_2  = np.zeros(len(alpha3))
EE_2 = np.zeros(len(alpha3))
r12_2 =  np.zeros(len(alpha3))

for counter, i in enumerate(alpha3):
    data = get_data("../data/analytical_E1_%3.2f.txt" % i, ["MCC", "E", "EE", "r12"])

    E_2[counter]   = np.mean(np.array(data["E"]))
    EE_2[counter]  = np.mean(np.array(data["EE"]))
    r12_2[counter] = np.mean(np.array(data["r12"]))

plt.figure(1)
plt.plot(alpha, E, "o", label="noninteracting")
plt.plot(alpha3, E_2, "o", label="interacting")

plt.legend(loc="best", fontsize=13)
plt.xlabel(r"Variational parameter $\alpha$", fontsize=14)
plt.ylabel(r"Average energy $\langle E \rangle $ [a.u] ", fontsize=14)
plt.savefig("../figures/energy_c.pdf", bbox_inches = "tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/energy_c.pdf", "../figures/energy_c.pdf"))

plt.figure(2)
plt.plot(alpha, EE-E*E, "o", label="noninteracting")
plt.plot(alpha3, EE_2-E_2*E_2, "o", label="interacting")

plt.legend(loc="best", fontsize=13)
plt.xlabel(r"Variational parameter $\alpha$", fontsize=14)
plt.ylabel(r"Variance energy $σ_E^2$ [a.u$^2$]", fontsize=14)

plt.savefig("../figures/variance_c.pdf", bbox_inches = "tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/variance_c.pdf", "../figures/variance_c.pdf"))

plt.show()

variance = EE - E*E
variance_2 = EE_2 - E_2*E_2

print("optimal alpha non interacting:", alpha[np.argmin(E)], "with energy: ", min(E))
print("optimal alpha non interacting:", alpha3[np.argmin(E_2)], "with energy: ", min(E_2))
print("the average deviation ")
print()
print("optimal alpha non interacting:", alpha[np.argmin(variance)], "with variance: ", min(variance))
print("optimal alpha non interacting:", alpha3[np.argmin(variance_2)], "with variance: ", min(variance_2))

omega = [0.05, 0.5, 1.0]


for counter, i in enumerate(omega):
    data_1 = get_data("../data/E_0_expectationvalues_1.00_%3.2f.txt" % i, ["MCC", "E", "EE", "r12"])
    data_2 = get_data("../data/E_1_expectationvalues_0.87_%3.2f.txt" % i, ["MCC", "E", "EE", "r12"])

    start = int(0.05*len(np.array(data_1["r12"])))

    print("omega: ", i)
    print("average_E0: ", np.mean(np.array(data_1["r12"])[start:]), "standrad deviation: ", np.std(np.array(data_1["r12"])[start:]) )
    print("average_E1: ", np.mean(np.array(data_2["r12"])[start:]), "standrad deviation: ", np.std(np.array(data_2["r12"])[start:])  )
    print()
