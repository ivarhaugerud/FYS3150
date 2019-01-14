import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

plt.style.use("bmh")
sns.color_palette("hls", 1)

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

def expec_value_E(T):
    J = 1
    beta = 1/T
    return -2*J*np.sinh(8*beta*J)/(3+np.cosh(8*beta*J))

def expec_value_EE(T):
    J = 1
    beta = 1/T
    Z = 4*(3+np.cosh(8*beta*J))
    return 64*J*J*np.cosh(8*beta*J)/Z

def expec_value_M(T):
    return 0*T

def expec_value_MM(T):
    J = 1
    beta = 1/T
    Z = 4*(3+np.cosh(8*beta*J))
    return 8*(np.exp(beta*8*J)+1)/Z # I actually find 32, oh well, this is per particle at least

def heat_cap(T):
    J = 1
    k = 1
    beta = 1/T
    Z = 4*(3+np.cosh(8*beta*J))
    return 1/(k*T*T)*(expec_value_EE(T)-expec_value_E(T)**2)


def sucept(T):
    J = 1
    k = 1
    beta = 1/T
    Z = 4*(3+np.cosh(8*beta*J))
    return 1/(k*T)*(expec_value_MM(T)-expec_value_M(T)**2)

temp = np.arange(0.5, 4.00, 0.2)
print(temp)
lent = len(temp)
average_E = np.zeros((lent, 2))
average_M = np.zeros((lent, 2))
average_EE = np.zeros((lent, 2))
average_MM = np.zeros((lent, 2))

counter = 0

for T in temp:
    results = get_data("../data/test_%3.2f_expect.txt" % (T), ["E", "EE", "M", "MM", "absM"])
    halfway_index = int(0.5*len(np.array(results["E"])))

    average_E[counter, 0] = np.mean(np.array(results["E"])[halfway_index:])
    average_M[counter, 0] = np.mean(np.array(results["M"])[halfway_index:])
    average_EE[counter, 0] = np.mean(np.array(results["EE"])[halfway_index:])
    average_MM[counter, 0] = np.mean(np.array(results["MM"])[halfway_index:])

    average_E[counter, 1] = np.std(np.array(results["E"])[halfway_index:])
    average_M[counter, 1] = np.std(np.array(results["M"])[halfway_index:])
    average_EE[counter, 1] = np.std(np.array(results["EE"])[halfway_index:])
    average_MM[counter, 1] = np.std(np.array(results["MM"])[halfway_index:])

    counter += 1

plt.errorbar(temp, average_E[:, 0], fmt='o', yerr=average_E[:, 1], label="numerical E")
plt.plot(temp, expec_value_E(temp), label="analytical E")
plt.show()

plt.errorbar(temp, average_EE[:, 0], fmt='o', yerr=average_EE[:, 1], label="numerical M")
plt.plot(temp, expec_value_EE(temp), label="analytical")
plt.show()

plt.errorbar(temp, average_M[:, 0], fmt='o', yerr=average_M[:, 1], label="numerical E")
plt.plot(temp, expec_value_M(temp), label="analytical E")
plt.show()

plt.errorbar(temp, average_MM[:, 0], fmt='o', yerr=average_MM[:, 1], label="numerical M")
plt.plot(temp, expec_value_MM(temp), label="analytical")
plt.show()

plt.plot(temp, heat_cap(temp))
plt.plot(temp, (average_EE[:, 0]-average_E[:, 0]**2)/(temp*temp), "o")
plt.show()

plt.plot(temp, sucept(temp))
plt.plot(temp, (average_MM[:, 0]-average_M[:, 0]**2)/(temp), "o")
plt.show()
