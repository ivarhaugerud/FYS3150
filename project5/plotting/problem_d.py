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

alpha = np.arange(0.05, 1.5, 0.05)
E  = np.zeros(len(alpha))
EE = np.zeros(len(alpha))

for counter, i in enumerate(alpha):
    data = get_data("../data/find_beta_expectationvalues_0.89_%3.2f.txt" % i, ["MCC", "E", "EE", "r12"])

    E[counter]  = np.mean(np.array(data["E"])[5000:])
    EE[counter] = np.mean(np.array(data["EE"])[5000:])

plt.figure(1)
plt.plot(alpha, E, "o")
plt.legend(loc="best", fontsize=13)
plt.xlabel(r"Variational parameter $\alpha$", fontsize=14)
plt.ylabel(r"Average energy $\langle E \rangle$ [a.u]", fontsize=14)
plt.savefig("../figures/benchmark.pdf", bbox_inches = "tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/benchmark.pdf", "../figures/benchmark.pdf"))

plt.figure(2)
plt.plot(alpha, EE-E*E, "o")
plt.legend(loc="best", fontsize=13)
plt.xlabel(r"Variational parameter $\alpha$", fontsize=14)
plt.ylabel(r"Variance energy $σ_E^2$ [a.u$^2$]", fontsize=14)
plt.savefig("../figures/benchmark_variance.pdf", bbox_inches = "tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/benchmark_variance.pdf", "../figures/benchmark_variance.pdf"))

plt.show()
