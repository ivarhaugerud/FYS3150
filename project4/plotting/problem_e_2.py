import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os

plt.style.use("bmh")
sns.color_palette("hls", 1)

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

temp = np.arange(2.31, 2.50, 0.01)
L_vals = np.arange(40, 70, 20)

average_E    = np.zeros(len(temp))
average_EE   = np.zeros(len(temp))
average_M    = np.zeros(len(temp))
average_MM   = np.zeros(len(temp))
average_Mabs = np.zeros(len(temp))

counter = 0

for L in L_vals:
    for T in temp:
        results = get_data("../new_data/%i_%3.2f_expect.txt" % (L, T), ["MCs", "E", "EE", "M", "MM", "absM"])

        half_way_index = int(len(results["E"])*0.5)
        number_MCs = np.array(results["MCs"])[-1]

        average_E[counter]    = np.mean(np.array(results["E"])[half_way_index:])
        average_EE[counter]   = np.mean(np.array(results["EE"])[half_way_index:])
        average_M[counter]    = np.mean(np.array(results["M"])[half_way_index:])
        average_MM[counter]   = np.mean(np.array(results["MM"])[half_way_index:])
        average_Mabs[counter] = np.mean(np.array(results["absM"])[half_way_index:])

        counter += 1

        #plt.plot(results["absM"])
    #plt.show()

    counter = 0
    plt.figure(1)
    plt.plot(temp, average_Mabs, "o", label=r"$L$=%3.0f" % L)


    plt.figure(2)
    plt.plot(temp, (average_EE-average_E**2)/(temp*temp)*L*L, "o", label=r"$T$=%3.0f" % L)

    plt.figure(3)
    plt.plot(temp, average_E, "o", label=r"$L$=%3.0f" % L)

    plt.figure(4)
    plt.plot(temp, (average_MM-average_M**2)/(temp), "o", label=r"$L$=%3.0f" % L)

temp = np.arange(2.31, 2.50, 0.01)
L_vals = np.arange(80, 110, 20)

average_E    = np.zeros(len(temp))
average_EE   = np.zeros(len(temp))
average_M    = np.zeros(len(temp))
average_MM   = np.zeros(len(temp))
average_Mabs = np.zeros(len(temp))

counter = 0

for L in L_vals:
    for T in temp:
        results = get_data("../new_data/%i_%4.3f_expect.txt" % (L, T), ["MCs", "E", "EE", "M", "MM", "absM"])

        half_way_index = int(len(results["E"])*0.5)
        number_MCs = np.array(results["MCs"])[-1]

        average_E[counter]    = np.mean(np.array(results["E"])[half_way_index:])
        average_EE[counter]   = np.mean(np.array(results["EE"])[half_way_index:])
        average_M[counter]    = np.mean(np.array(results["M"])[half_way_index:])
        average_MM[counter]   = np.mean(np.array(results["MM"])[half_way_index:])
        average_Mabs[counter] = np.mean(np.array(results["absM"])[half_way_index:])

        counter += 1

        #plt.plot(results["absM"])
    #plt.show()

    counter = 0
    plt.figure(1)
    plt.plot(temp, average_Mabs, "o", label=r"$L$=%3.0f" % L)


    plt.figure(2)
    plt.plot(temp, (average_EE-average_E**2)/(temp*temp)*L*L, "o", label=r"$L$=%3.0f" % L)

    plt.figure(3)
    plt.plot(temp, average_E, "o", label=r"$L$=%3.0f" % L)

    plt.figure(4)
    plt.plot(temp, (average_MM-average_M**2)/(temp), "o", label=r"$L$=%3.0f" % L)


plt.figure(2)
plt.xlabel(r"temperature $T$")
plt.ylabel(r"heat capacity $C_V$")
plt.legend(loc="best")
plt.savefig("../figures/problem_e_heat_cap.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/problem_e_heat_cap.pdf", "../figures/problem_e_heat_cap.pdf"))

plt.figure(1)
plt.xlabel(r"temperature $T$")
plt.ylabel(r"absolute magnetization $|M|$")
plt.legend(loc="best")
plt.savefig("../figures/problem_e_absM.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/problem_e_absM.pdf", "../figures/problem_e_absM.pdf"))

plt.figure(3)
plt.xlabel(r"temperature $T$")
plt.ylabel(r"energy $E$")
plt.legend(loc="best")
plt.savefig("../figures/problem_e_E.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/problem_e_E.pdf", "../figures/problem_e_E.pdf"))

plt.figure(4)
plt.xlabel(r"temperature $T$")
plt.ylabel(r"magnetic suceptibility $\chi$")
plt.legend(loc="best")
plt.savefig("../figures/problem_e_MS.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/problem_e_MS.pdf", "../figures/problem_e_MS.pdf"))
plt.show()
