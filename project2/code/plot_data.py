import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

plt.style.use("bmh")
sns.color_palette("hls", 1)

#Importing libarys, and make the plotting look better

def get_data(filename, variables):
    df = pd.read_csv(filename,\
                    delim_whitespace=True, \
                    engine='python', \
                    names=variables)
    return df
    #using pandas to read the data files

#linear regression function
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

iterations_data = {} #data for iterations
iterations_data = get_data("data/iterations_jacoby.txt", ["n", "iterations"])

#calculate linear regression
"""
m, c, delta_m, delta_c = linear_regresion(np.log10(iterations_data["n"]), np.log10(iterations_data["iterations"]))

#print the calculated values
print("For best fit line we calculate:")
print("   Slope:    %5.4f pm  %4.3f" % (m, delta_m))
print("Constant:    %5.4f pm  %4.3f" % (c, delta_c))

plt.xscale("log")
plt.yscale("log")
plt.plot(iterations_data["n"], iterations_data["iterations"], 'go', label="datapoints")
plt.plot(iterations_data["n"], 10**(m*np.log10(iterations_data["n"]))*10**c, label="linear best fit, slope: %3.2f" % m)
plt.fill_between(iterations_data["n"], 10**((m+delta_m)*np.log10(iterations_data["n"]))*10**(c+delta_c), 10**((m-delta_m)*np.log10(iterations_data["n"]))*10**(c-delta_c), alpha=0.3, label=r"uncertainty best fit: $\pm$ 0.01")

plt.xlabel(r"Size of matrix $n$")
plt.ylabel(r"Number of iterations to reach $\epsilon$")
plt.legend(loc="best")
#plt.savefig("../figures/iterations_needed.pdf")
plt.show()
"""

plt.figure(1, figsize=(7,4.5))
wavefunc_one_e = get_data("data/wavefunc_one_e.txt", ["x", "wf1", "wf2", "wf3", "wf4", "wf5"])
wavefunc_two_e = get_data("data/wavefunc_two_e_1.00.txt", ["x", "wf1", "wf2", "wf3", "wf4", "wf5"])

plt.plot(wavefunc_one_e["x"], wavefunc_one_e["wf1"]**2/np.trapz(wavefunc_one_e["wf1"]**2, wavefunc_one_e["x"]), label="one electron")
plt.plot(wavefunc_two_e["x"], wavefunc_two_e["wf1"]**2/np.trapz(wavefunc_two_e["wf1"]**2, wavefunc_two_e["x"]), label=r"two electrons $\omega_r=1$")

plt.ylabel(r"Radial solution $\vert u(r) \vert^2 $", fontsize=15)
plt.xlabel(r"Postion $\rho$", fontsize=15)
plt.legend(loc="best", fontsize=13)
plt.savefig("../figures/difference_interaction.pdf")
plt.show()

plt.figure(5, figsize=(7,4.5))
plt.plot(wavefunc_one_e["x"], wavefunc_one_e["wf1"]**2/np.trapz(wavefunc_one_e["wf1"]**2, wavefunc_one_e["x"]), label=r"$\psi_0$")
plt.plot(wavefunc_one_e["x"], wavefunc_one_e["wf2"]**2/np.trapz(wavefunc_one_e["wf2"]**2, wavefunc_one_e["x"]), label=r"$\psi_1$")
plt.plot(wavefunc_one_e["x"], wavefunc_one_e["wf3"]**2/np.trapz(wavefunc_one_e["wf3"]**2, wavefunc_one_e["x"]), label=r"$\psi_2$")

plt.legend(loc="best", fontsize=13)
plt.xlabel(r"position $\rho_i$", fontsize=15)
plt.ylabel(r"radial solution $\vert u(\rho)\vert^2$", fontsize=15)
plt.savefig("../figures/wavefunc_one_e.pdf")

plt.show()
omega = np.array([0.01, 0.50, 1.00, 5.00])

plt.figure(10, figsize=(7,4.5))
plt.xscale("log")
for i in range(len(omega)):
    wavefunc_two_e = get_data("data/wavefunc_two_e_%3.2f.txt" % omega[i], ["x", "wf1", "wf2", "wf3", "wf4", "wf5"])
    plt.plot(wavefunc_two_e["x"], wavefunc_two_e["wf1"]**2/np.trapz(wavefunc_two_e["wf1"]**2, wavefunc_two_e["x"]), label=r"$\omega$=%3.2f" % omega[i])

plt.legend(loc="best", fontsize=13)
plt.xlabel(r"position $\rho_i$", fontsize=15)
plt.ylabel(r"radial solution $\vert u(\rho)\vert^2$", fontsize=15)
#plt.savefig("../figures/wavefunc_many_omega.pdf")

plt.show()
