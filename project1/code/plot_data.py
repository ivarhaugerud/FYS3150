import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
plt.style.use("bmh")
sns.color_palette("hls", 1)
#Importing libarys, and make the plotting look better

def foo(x):
    return  1 - (1 - np.exp(-10))*x - np.exp(-10*x)
    #Analytic expression for the derivative

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

###
#AVERAGE TIME RELATIONSHIP BETWEEN 8N AND 4N ALGORITHM
###

number_of_files_read = 6
n_no_exp = np.linspace(1, number_of_files_read, 2*number_of_files_read-1)
n = np.power(10, n_no_exp)
logn = n_no_exp
max_deviation = np.zeros(number_of_files_read*2-1)
#Defining arrays for plotting and reading

mean_time_general = np.zeros(number_of_files_read*2-1)
mean_time_special = np.zeros(number_of_files_read*2-1)
std_general = np.zeros(number_of_files_read*2-1)
std_special = np.zeros(number_of_files_read*2-1)
relations_ship_times = np.zeros(number_of_files_read*2-1)
std_relationship = np.zeros(number_of_files_read*2-1)
#Defining arrays for plotting and reading

time_data = {}
for i in range(2*number_of_files_read-1):
    #Loop over and collects all of the data, and calculate the interesting values
    time_data[str(n[i])] = get_data("data/time_"+str(int(n[i]))+".txt", ["general", "special"])
    mean_time_general[i] = np.mean(time_data[str(n[i])]["general"])
    mean_time_special[i] = np.mean(time_data[str(n[i])]["special"])
    std_special[i] = np.std(time_data[str(n[i])]["general"])
    std_general[i] = np.std(time_data[str(n[i])]["general"])

    relations_ship_times[i] = np.mean(mean_time_general[i]/mean_time_special[i])
    std_relationship[i] = np.sqrt((std_general[i]/mean_time_general[i])**2 + (std_special[i]/mean_time_special[i])**2)

#Variabels needed for plotting the time relationship between 4n and 8n
number_of_runs = len(time_data["10.0"]["general"])
stabelize = -6
average_over_last = np.mean(relations_ship_times[stabelize:])
std = np.std(relations_ship_times[stabelize:])

#Plotting
ax = plt.figure(1)
plt.xscale("log")
plt.fill_between(n[stabelize:], (average_over_last+np.std(relations_ship_times[stabelize:]))*np.ones(len(n[stabelize:])), (average_over_last-np.std(relations_ship_times[stabelize:]))*np.ones(len(n[stabelize:])), alpha=0.3)
plt.errorbar(n, relations_ship_times, fmt='o', yerr=std_relationship)
plt.xlabel(r"$n$", fontsize=15)
plt.ylabel(r"Runtime relationship $\frac{\langle t_{8n} \rangle }{\langle t_{4n} \rangle }$", fontsize=13)
plt.savefig("../figures/time_difference.pdf")
plt.show()

print("Average time relationship: %5.4f pm %5.4f" % ((average_over_last), (std)))

###
#AVERAGE TIME LU DECOMP
###

#Defining arrays for plotting and reading
LU_files = 7
n = np.power(10, np.linspace(0.5, LU_files/2, LU_files))

time_data_LU = {}
mean_time = np.zeros(LU_files)
std_time  = np.zeros(LU_files)

for i in range(LU_files):
    #Reading data using pandas
    time_data_LU[str(n[i])] = get_data("data/time_LU_"+str(int(n[i]))+".txt", ["time"])
    mean_time[i] = np.mean(time_data_LU[str(n[i])]["time"])
    std_time[i] = np.std(time_data_LU[str(n[i])]["time"])

#Plotting
ax = plt.figure(1)
plt.xscale("log")
plt.yscale("log")
m, c, delta_m, delta_c = linear_regresion(np.log10(n[1:]), np.log10(mean_time[1:]))
plt.errorbar(n, mean_time, fmt='*', yerr=std_time, label="measurments")
plt.plot(n[1:],  10**c*10**(m*np.log10(n[1:])), "--", label="best fit, slop: %3.2f pm %3.2f" % (m, delta_m))
plt.xlabel(r"$n$", fontsize=15)
plt.ylabel(r"Runtime of LU algorithm [s]", fontsize=15)
plt.legend(loc="best")
plt.savefig("../figures/time_difference_LU.pdf")

#printing the interesting values
print("Slope:        %5.4f           Constant term:    %5.4f" % (m, c,))
print("Delta Slope:  %5.4f     Delta Constant term:    %5.4f" % (delta_m, delta_c))
print("If log(N)=5 it would take aprox: %5.4f months to run code" % ( (10**c * 10**(m*5))/(60*60*24*365/12.0)))
plt.show()

###
#DEVIATION ANALYTICAL AND NUMERICAL
###

#Defining arrays for plotting and reading
number_of_files_read = 7
n_no_exp = np.linspace(1, number_of_files_read, 2*number_of_files_read-1)
n = np.power(10, n_no_exp)
logn = n_no_exp
max_deviation = np.zeros(int(number_of_files_read*2-1))

data = {}
for i in range(int(2*number_of_files_read-1)):
    data[str(int(n[i]))] = (get_data("data/data_"+str(int(n[i]))+".txt", ['x', "f", "deviation", "log_deviation"]))

#Defining arrays for plotting and reading
x = np.linspace(0, 1, int(1e6))
analytic = foo(x)

#plotting
plt.figure(1)
plt.plot(x, analytic, '-r', label="analytical")
for i in range(int(number_of_files_read)):
    plt.plot(data[str(int(n[i]))]["x"], data[str(int(n[i]))]["f"], "--", label=r"$n=10^%i$" % (i+1))
plt.legend(loc="best")
plt.xlabel(r"$x$", fontsize=18)
plt.ylabel(r"$u(x)$", fontsize=15)
plt.savefig("../figures/graphs.pdf")

for i in range(int(2*number_of_files_read-1)):
    #reads the maximum deviation
    max_deviation[i] = np.max(data[str(int(n[i]))]["log_deviation"][1:-1])

#using linear regression on our data
m, c, delta_m, delta_c = linear_regresion(logn[0:-3], max_deviation[0:-3])

#plots it
plt.figure(3)
plt.xscale("log")
plt.yscale("log")
plt.plot(n, 10**max_deviation, "k*", label="measurments")
plt.plot(n[0:-3], 10**c * 10**(m*np.log10(n[0:-3])), "--", label=r"best fit, slope: %3.2f $\pm$ %3.2f" % (m, delta_m))
plt.xlabel(r"$n$", fontsize=15)
plt.ylabel(r"$Max \, (\epsilon)$", fontsize=15)
plt.legend(loc="best")
plt.savefig("../figures/log_difference.pdf")
plt.show()
