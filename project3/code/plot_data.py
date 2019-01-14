import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import sys

plt.style.use("bmh")
sns.color_palette("husl", 8)
#sns.hls_palette(8, l=.3, s=.8)
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


#PLOT ESCAPE VELOCITY
"""
respective_T_max = [5000, 10000, 20000, 40000, 50000, 50000]
counter = 0
fig = plt.figure()
ax = fig.add_subplot(111)

for i in [8.878, 8.888, 8.8835, 8.88575, 8.8846, 8.88515]: #8.857, 8.868,
    position_x = get_data("data/escape_velocity/%5.4f" % i +"_VV__x.txt", ["sun", "earth"])
    position_y = get_data("data/escape_velocity/%5.4f" % i +"_VV__y.txt", ["sun", "earth"])
    t = np.linspace(0, respective_T_max[counter], int(len(position_x)-2))
    plt.plot(t, np.sqrt(np.array(position_x["earth"])[2:]**2 + np.array(position_y["earth"])[2:]**2), label=r"$v_0=$%5.4f" % i)
    counter += 1

plt.yscale("log")
plt.xscale("log")

plt.xlabel(r"time [years]")
plt.ylabel(r"radial distance [AU]")
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12),
          ncol=3, fancybox=True, shadow=True)
plt.savefig("../figures/escape_velocity.pdf")#is now perfect
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/escape_velocity.pdf", "../figures/escape_velocity.pdf"))
plt.show()
"""
#PLOT EARTH SUN ORBIT
"""
plt.figure(1)
objects = ["sun", "earth"]

position_x = get_data("data/sun_at_rest_VV__x.txt", objects)
position_y = get_data("data/sun_at_rest_VV__y.txt", objects)
plt.plot(position_x["earth"], position_y["earth"], label="earth")
plt.plot(position_x["sun"], position_y["sun"],"o", label="sun")

plt.axis("equal")
plt.xlabel(r"$x-$position [AU]")
plt.ylabel(r"$y$-position [AU]")
plt.legend(loc="best", fontsize=12)
#plt.savefig("../figures/eart_sun_orbit.pdf", bbox_inches="tight") #this plot is perfect
plt.show()
"""
#CONSERVATION OF ENERGY AND ANGULAR MOMENTUM
"""
plt.figure(2)
types_energy = ["total", "kinetic", "potential", "ang_mom"]
energy_VV = get_data("data/sun_at_rest_VV_energy.txt", types_energy)
energy_FE = get_data("data/sun_at_rest_FE_energy.txt", types_energy)
t = np.linspace(0, 50, len(energy_VV))

plt.plot(t, (energy_VV["total"]-np.array(energy_VV["total"])[0])/np.array(energy_VV["total"])[0], label=r"$E_{vv}$")
plt.plot(t, (energy_VV["ang_mom"]-np.array(energy_VV["ang_mom"])[0])/np.array(energy_VV["ang_mom"])[0], "--", label=r"$L_{vv}$")

plt.plot(t, abs(energy_FE["total"]-np.array(energy_FE["total"])[0])/abs(np.array(energy_FE["total"])[0]), label=r"$E_{fe}$")
plt.plot(t, (energy_FE["ang_mom"]-np.array(energy_FE["ang_mom"])[0])/np.array(energy_FE["ang_mom"])[0], "--", label=r"$L_{fe}$")

plt.legend(loc="best", fontsize=12)
plt.xlabel("time [years]")
plt.ylabel(r"absolute relative deviation from start")
plt.yscale("log")
plt.savefig("../figures/deviation_E_L.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/deviation_E_L.pdf", "../figures/deviation_E_L.pdf"))
plt.show()
"""
#FINAL 3D PLOT
"""
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')

objects = ["sun", "mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune", "pluto", "moon", "io", "europa", "halley"]

#ax.set_xlim(-10, 3)
#ax.set_ylim(-10, 3)
#ax.set_zlim(-3, 5)
position_x = get_data("data/simulator_VV__x.txt", objects)
position_y = get_data("data/simulator_VV__y.txt", objects)
position_z = get_data("data/simulator_VV__z.txt", objects)

for object in objects:
    ax.plot(position_x[object], position_y[object], position_z[object], label=object)

#box = ax.get_position()
#ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
#          ncol=3, fancybox=True, shadow=True)
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

# Put a legend below current axis
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=7, fontsize=8)

ax.set_xlabel(r"$x-$position [AU]")
ax.set_ylabel(r"$y-$position [AU]")
ax.set_zlabel(r"$z-$position [AU]")
plt.savefig("../figures/total_plot.pdf", bbox_to_anchor="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/total_plot.pdf", "../figures/total_plot.pdf"))

plt.show()
"""
#CHANGING MASS OF JUPITER
"""
def highResPoints(x,y,factor=10):
    r = [0]
    for i in range(1,len(x)):
        dx = x[i]-x[i-1]
        dy = y[i]-y[i-1]
        r.append(np.sqrt(dx*dx+dy*dy))
    r = np.array(r)

    # rtot is a cumulative sum of r, it's used to save time
    rtot = []
    for i in range(len(r)):
        rtot.append(r[0:i].sum())
    rtot.append(r.sum())

    dr = rtot[-1]/(NPOINTS*RESFACT-1)
    xmod=[x[0]]
    ymod=[y[0]]
    rPos = 0 # current point on walk along data
    rcount = 1
    while rPos < r.sum():
        x1,x2 = x[rcount-1],x[rcount]
        y1,y2 = y[rcount-1],y[rcount]
        dpos = rPos-rtot[rcount]
        theta = np.arctan2((x2-x1),(y2-y1))
        rx = np.sin(theta)*dpos+x1
        ry = np.cos(theta)*dpos+y1
        xmod.append(rx)
        ymod.append(ry)
        rPos+=dr
        while rPos > rtot[rcount+1]:
            rPos = rtot[rcount+1]
            rcount+=1
            if rcount>rtot[-1]:
                break

    return xmod,ymod

NPOINTS = 10
RESFACT=10
MAP='viridis' # choose carefully, or color transitions will not appear smoooth

fig = plt.figure()
ax3 = fig.add_subplot(111) # high resolution color map

position_x = get_data("data/jupiter_1000_VV__x.txt", ["sun", "earth", "jupiter"])
position_y = get_data("data/jupiter_1000_VV__y.txt", ["sun", "earth", "jupiter"])

xHiRes,yHiRes = highResPoints(position_x["earth"],position_y["earth"],RESFACT)
npointsHiRes = int(len(xHiRes))
cm = plt.get_cmap(MAP)
ax3.set_color_cycle([cm(1.*i/(npointsHiRes-1))
                     for i in range(npointsHiRes-1)])
for i in range(npointsHiRes-1):
    ax3.plot(xHiRes[i:i+2],yHiRes[i:i+2], "-")

xHiRes,yHiRes = highResPoints(position_x["jupiter"],position_y["jupiter"],RESFACT)
npointsHiRes = len(xHiRes)
cm = plt.get_cmap(MAP)

ax3.set_color_cycle([cm(1.*i/(npointsHiRes-1))
                     for i in range(npointsHiRes-1)])

for i in range(npointsHiRes-1):
    ax3.plot(xHiRes[i:i+2],yHiRes[i:i+2], "o")

plt.axis("equal")
plt.plot([0,0],[0,0], "ro")

plt.xlabel(r"$x-$position [AU]")
plt.ylabel(r"$y-$position [AU]")
plt.axis("equal")
plt.savefig("../figures/jupiter_1000_mass.pdf", bbox_to_anchor="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/jupiter_1000_mass.pdf", "../figures/jupiter_1000_mass.pdf"))

plt.show()
"""
#PLOT RADIAL DISTANCE DIFFERENT JUPITER MASS
"""
m = ["1", "10", "100"]
line_style = ["-", "--", "-."]#["-", "--", "--"]
counter = 0

position_x_no_jup = get_data("data/sun_at_rest_VV__x.txt", ["sun", "earth"])
position_y_no_jup = get_data("data/sun_at_rest_VV__y.txt", ["sun", "earth"])

fig = plt.figure(10)
ax = fig.add_subplot(111)

for i in (m):
    position_x = get_data("data/jupiter_"+i+"_VV__x.txt", ["sun", "earth", "jupiter"])
    position_y = get_data("data/jupiter_"+i+"_VV__y.txt", ["sun", "earth", "jupiter"])
    t = np.linspace(0, 24, len(position_x))

    plt.plot(t, np.sqrt(np.array(position_x["earth"])**2 + np.array(position_y["earth"])**2)-np.sqrt(np.array(position_x_no_jup["earth"])**2+np.array(position_y_no_jup["earth"])**2), line_style[counter], label=r"$M_J\,\,*=10^{%1.0f}$" % (np.log10(int(i))))
    counter += 1

ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
          ncol=3, fancybox=True, shadow=True, fontsize=12)
plt.xlabel("time [years]")
plt.ylabel("change in radial distance of earth [AU]")
plt.savefig("../figures/radial_distance_jupiter.pdf", bbox_to_anchor="tight") #the figure is now perfect
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/radial_distance_jupiter.pdf", "../figures/radial_distance_jupiter.pdf"))
plt.show()
"""
#FOURIER
"""
counter = 0

fig = plt.figure(11, dpi=100)
ax = fig.add_subplot (111)

for i in (m):
    position_x = get_data("data/jupiter_"+i+"_VV__x.txt", ["sun", "earth", "jupiter"])
    position_y = get_data("data/jupiter_"+i+"_VV__y.txt", ["sun", "earth", "jupiter"])

    N = len(position_x)
    x = np.linspace(0, 24, N)
    Sf = x[1] - x[0]

    y = (np.sqrt(np.array(position_x["earth"])**2 + np.array(position_y["earth"])**2)-1) #np.sqrt(np.array(position_x_no_jup["earth"])**2+np.array(position_y_no_jup["earth"])**2))#*np.exp(-(x-12)**2/32)

    #plt.plot(y)
    #plt.show()

    freqencies = np.fft.fftfreq(N, Sf)
    amp = np.absolute(np.fft.fft(y))/N

    counter += 1

    ax.plot(freqencies[:int(N/2)], amp[:int(N/2)], label="m="+str(i))
ax.set_xlabel(r"Frequency [1/years]")
ax.set_ylabel(r"Fourier koeffisient [AU]")
#plt.axis([0, 2, -0.00001, 0.0003])
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
          ncol=3, fancybox=True, shadow=True, fontsize=12)
plt.savefig("../figures/fourier.pdf")
plt.show ()
"""
#PLOT DIFFERENCE SUN AT REST
"""
fig = plt.figure(6)
ax = fig.add_subplot(111)
objects = ["sun", "earth", "jupiter"]
position_x_sun_at_rest = get_data("data/jupiter_1_VV__x.txt", objects)
position_y_sun_at_rest = get_data("data/jupiter_1_VV__y.txt", objects)
position_x = get_data("data/simulator_VV__x.txt", objects)
position_y = get_data("data/simulator_VV__y.txt", objects)

objects = ["earth", "sun", "jupiter"]
for i in objects:
    time = np.linspace(0, 24, len(position_x[i]))
    plt.plot(time, 100*(np.sqrt(np.array(position_x[i])**2 + np.array(position_y[i])**2) - np.sqrt(np.array(position_x_sun_at_rest[i])**2 + np.array(position_y_sun_at_rest[i])**2)), label=i)

plt.xlabel(r"time [years]")
plt.ylabel(r"change in radial position[$10^{-2}$AU]")
box = ax.get_position()
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
          ncol=3, fancybox=True, shadow=True, fontsize=12)
plt.savefig("../figures/sun_movement_effect.pdf", bbox_inches="tight") #plot is perfect
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/sun_movement_effect.pdf", "../figures/sun_movement_effect.pdf"))
plt.show()
"""
#PLOT TIME RELATION
"""
stable_values = 27

plt.figure(3)
plt.xscale("log")
FE_time = get_data("data/FE_timer.txt", ["n", "T"])
VV_time = get_data("data/VV_timer.txt", ["n", "T"])
n = np.array(FE_time["n"])
exponents = np.log10(n)


plt.plot(n, np.array(FE_time["T"]), "o")
m, c, delta_m, delta_c = linear_regresion(exponents[-stable_values:], np.log10(FE_time["T"][-stable_values:]))
plt.plot(n[-stable_values:], 10**c*10**(m*exponents[-stable_values:]), label="FE-slope:%3.3f(%1.0f)" % (m, delta_m*1000))

plt.plot(n, np.array(VV_time["T"]), "o")
m, c, delta_m, delta_c = linear_regresion(exponents[-stable_values:], np.log10(VV_time["T"][-stable_values:]))
plt.plot(n[-stable_values:], 10**c*10**(m*exponents[-stable_values:]), label="VV-slope:%3.3f(%1.0f)" % (m, delta_m*1000))

plt.yscale("log")
plt.legend(loc="best")
plt.ylabel(r"time used [s]")
plt.xlabel(r"number of datapoints $N$")
plt.legend(loc="best", fontsize=12)
#plt.savefig("../figures/timer_VV_and_FE.pdf", bbox_inches="tight") #this figure is perfect now
#os.system('pdfcrop %s %s &> /dev/null &'%("../figures/timer_VV_and_FE.pdf", "../figures/timer_VV_and_FE.pdf"))
plt.show()

time_relation = np.array(FE_time["T"])[-stable_values:]/np.array(VV_time["T"])[-stable_values:]
print("Average time relationship FE/VV:", np.mean(time_relation))
print("Standard deviation of time relationship FE/VV:", np.std(time_relation))
"""

#DEVIATION AFTER ONE REVOLUTION
"""
plt.figure(4)

deviation_FE = []
deviation_VV = []
exponents = np.linspace(3, 8, 21)
numbers   = np.power(10, exponents)
for i in exponents:
    x = get_data("data/%3.2f_FE__x.txt" % i, ["sun", "earth"])
    y = get_data("data/%3.2f_FE__y.txt" % i, ["sun", "earth"])
    z = get_data("data/%3.2f_FE__z.txt" % i, ["sun", "earth"])
    delta_x_earth = np.array(x["earth"])[-1] - np.array(x["earth"])[0]
    delta_y_earth = np.array(y["earth"])[-1] - np.array(y["earth"])[0]
    delta_z_earth = np.array(z["earth"])[-1] - np.array(z["earth"])[0]

    deviation_FE.append(np.sqrt(delta_x_earth**2 + delta_y_earth**2 + delta_z_earth**2))

for i in exponents:
    x = get_data("data/%3.2f_VV__x.txt" % i, ["sun", "earth"])
    y = get_data("data/%3.2f_VV__y.txt" % i, ["sun", "earth"])
    z = get_data("data/%3.2f_VV__z.txt" % i, ["sun", "earth"])
    delta_x_earth = np.array(x["earth"])[-1] - np.array(x["earth"])[0]
    delta_y_earth = np.array(y["earth"])[-1] - np.array(y["earth"])[0]
    delta_z_earth = np.array(z["earth"])[-1] - np.array(z["earth"])[0]

    deviation_VV.append(np.sqrt(delta_x_earth**2 + delta_y_earth**2 + delta_z_earth**2))


plt.yscale("log")
plt.xscale("log")

m, c, delta_m, delta_c = linear_regresion(exponents, np.log10(np.array(deviation_FE)))
plt.plot(numbers, 10**c * 10**(m*exponents), label="FE, slope = %4.3f(%d)" % (m, round(delta_m*1000)))
plt.plot(numbers, np.array(deviation_FE), "o")
m, c, delta_m, delta_c = linear_regresion(exponents, np.log10(np.array(deviation_VV)))
plt.plot(numbers, 10**c * 10**(m*exponents), label="VV, slope = %4.3f(%d)" % (m, round(delta_m*1000)))
plt.plot(numbers, np.array(deviation_VV), "o")
plt.legend(loc="best")
plt.xlabel(r"number of datapoints $N$")
plt.ylabel(r"deviation initial vs final position [AU]")
plt.savefig("../figures/deviation_vs_n.pdf", bbox_inches="tight")
plt.show()
"""
#DIFFERENT EXPONENT IN NEWTONS GRAVIY
"""
plt.figure(123)
exponents = np.linspace(2.996, 2.998, 2)
for i in exponents:
    position_x = get_data("data/different_exponent_%5.4f_VV__x.txt" % i, ["sun", "earth"])
    position_y = get_data("data/different_exponent_%5.4f_VV__y.txt" % i, ["sun", "earth"])
    t = np.linspace(0, 50, len(position_x))

    plt.plot(t, np.sqrt(position_x["earth"]**2 + position_y["earth"]**2), label=r"$\beta=$"+str(i))

exponents = np.linspace(2.999, 3, 6)
for i in exponents:
    position_x = get_data("data/different_exponent_%5.4f_VV__x.txt" % i, ["sun", "earth"])
    position_y = get_data("data/different_exponent_%5.4f_VV__y.txt" % i, ["sun", "earth"])
    plt.plot(t, np.sqrt(position_x["earth"]**2 + position_y["earth"]**2), label=r"$\beta=$"+str(i))

plt.ylabel(r"radial distance [AU]")
plt.xlabel(r"time [years]")
plt.legend(loc="best", fontsize=12)
plt.axis([0, 50, 1-0.005, 1+0.09])
plt.savefig("../figures/betta_plot.pdf", bbox_to_anchor="tight")
plt.show()

#os.system('pdfcrop %s %s &> /dev/null &'%("../figures/betta_plot.pdf", "../figures/betta_plot.pdf"))
"""
#STUDYING PERIHELION
"""
perihelion_x = get_data("data/perihelion_relativity_x.txt", ["sun", "mercury"])
perihelion_y = get_data("data/perihelion_relativity_y.txt", ["sun", "mercury"])

perihelion_y = (np.array(perihelion_y["mercury"]))
perihelion_x = (np.array(perihelion_x["mercury"]))

position_x = get_data("data/relativistic_VV__x.txt", ["sun", "mercury"])
position_y = get_data("data/relativistic_VV__y.txt", ["sun", "mercury"])

plt.plot(position_x["mercury"], position_y["mercury"])
plt.plot(perihelion_x, perihelion_y)
plt.show()


thetap = np.array(np.arctan2(perihelion_y, perihelion_x))

plt.plot(thetap, "o")

t = np.linspace(0, 400, len(thetap))
m, c, delta_m, delta_c = linear_regresion(t, thetap)
print(m*100*(60*60)/np.pi*180, c, delta_m*100*(60*60)/np.pi*180, delta_c)
plt.show()

years = 400
cent  = float(years)/100

anal = 42.98/(60*60)*np.pi/180
num  = (thetap[-1]-thetap[0])/cent

print("Analytisk:      %3.4e" % (anal*(60*60)/np.pi*180))
print("Numerisk:       %3.4e" % (num*(60*60)/np.pi*180))
print("Relativt avvik: %3.2e" % (abs(num-anal)/anal))

mercury_years = years*365.25/87.969

print("Number of Mercurian years:   %d" % mercury_years)
print("Number of perihelions:       %d" % len(thetap))

"""
#PLOTTING EFFECTIVE POTENTIAL
"""

def Veff(r, b, c = 1, k = 1):
    return 0.5*k/(r**2) - c/((b-1.0)*r**(b-1.0))


r = np.linspace(0.5, 5, 1000)

c = 4*np.pi**2*3.0024584e-6
k = 4*np.pi**2*3.0024584e-6


norm = np.abs(np.min((Veff(r, 2, c, k))))
betas = [2, 2.9, 2.999, 2.999999, 3, 3.1]


for i in range(6):
    plt.subplot(2,3,i+1)
    plt.plot(r, Veff(r, betas[i], c, k)/norm)
    plt.title(r'$\beta = %.10g$' % betas[i])
    plt.ticklabel_format(axis='y',style='sci',scilimits=(0,1))

plt.subplot(231)
plt.ylabel(r'$E_{\oplus}$')
plt.subplot(234)
plt.xlabel(r'$r$[AU]')
plt.ylabel(r'$E_{\oplus}$')
plt.subplot(235)
plt.xlabel(r'$r$[AU]')
plt.subplot(236)
plt.xlabel(r'$r$[AU]')

plt.tight_layout(pad=0.1, w_pad=0.01, h_pad=1.0)
#plt.savefig("../figures/chaning_beta.pdf", bbox_inches="tight")
plt.show()
"""
#PLOT ENERGY for JUPITER SYSTEM
"""
exponents = np.linspace(3, 8.25, 22)
planets = ["sun", "earth", "jupiter"]
change_in_energy_VV = np.zeros(len(exponents))
change_in_energy_FE = np.zeros(len(exponents))

counter = 0

for i in exponents:
    energy = get_data("data/3_body/%5.4f_VV_energy.txt" % i, ["total", "kinetic", "potential", "ang_mom"])
    change_in_energy_VV[counter] = abs(np.array(energy["total"])[-1] - np.array(energy["total"])[0])/abs(np.array(energy["total"])[0])

    energy = get_data("data/3_body/%5.4f_VV__FE_energy.txt" % i, ["total", "kinetic", "potential", "ang_mom"])
    change_in_energy_FE[counter] = abs(np.array(energy["total"])[-1] - np.array(energy["total"])[0])/abs(np.array(energy["total"])[0])

    counter += 1

m, c, delta_m, delta_c = linear_regresion(exponents, np.log10(change_in_energy_VV))
plt.plot(np.power(10, exponents), 10**(c+m*exponents), label="VV - slope: %3.2f(%1.0f)" % (m, delta_m*100))
m, c, delta_m, delta_c = linear_regresion(exponents, np.log10(change_in_energy_FE))
plt.plot(np.power(10, exponents), 10**(c+m*exponents), label="FE - slope: %3.2f(%1.0f)" % (m, delta_m*100))

plt.plot(np.power(10, exponents), change_in_energy_VV, "o")
plt.plot(np.power(10, exponents), change_in_energy_FE, "o")

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"datapoints $N$")
plt.ylabel(r"relative change in energy")
plt.legend(loc="best", fontsize=12)

plt.savefig("../figures/conservation_energy.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/conservation_energy.pdf", "../figures/conservation_energy.pdf"))
plt.show()
"""
#ENERGY FOR DIFFERENT MJ
"""
exponents = np.linspace(3, 8, 21)
masses = ["1", "10", "100", "1000"]
planets = ["sun", "earth", "jupiter"]

change_in_energy_1    = np.zeros(len(exponents))
change_in_energy_10   = np.zeros(len(exponents))
change_in_energy_100  = np.zeros(len(exponents))
change_in_energy_1000 = np.zeros(len(exponents))

counter = 0

for i in exponents:
    energy = get_data("data/3_body/%5.4f_1_VV_energy.txt" % i, ["total", "kinetic", "potential", "ang_mom"])
    change_in_energy_1[counter] = abs(np.array(energy["total"])[-1] - np.array(energy["total"])[0])/abs(np.array(energy["total"])[0])

    energy = get_data("data/3_body/%5.4f_10_VV_energy.txt" % i, ["total", "kinetic", "potential", "ang_mom"])
    change_in_energy_10[counter] = abs(np.array(energy["total"])[-1] - np.array(energy["total"])[0])/abs(np.array(energy["total"])[0])

    energy = get_data("data/3_body/%5.4f_100_VV_energy.txt" % i, ["total", "kinetic", "potential", "ang_mom"])
    change_in_energy_100[counter] = abs(np.array(energy["total"])[-1] - np.array(energy["total"])[0])/abs(np.array(energy["total"])[0])

    energy = get_data("data/3_body/%5.4f_1000_VV_energy.txt" % i, ["total", "kinetic", "potential", "ang_mom"])
    change_in_energy_1000[counter] = abs(np.array(energy["total"])[-1] - np.array(energy["total"])[0])/abs(np.array(energy["total"])[0])

    counter += 1

plt.plot(np.power(10, exponents), change_in_energy_1,    "o", label=r"$M_j\cdot=10^{0}$")
plt.plot(np.power(10, exponents), change_in_energy_10,   "o", label=r"$M_j\cdot=10^{1}$")
plt.plot(np.power(10, exponents), change_in_energy_100,  "o", label=r"$M_j\cdot=10^{2}$")
plt.plot(np.power(10, exponents), change_in_energy_1000, "o", label=r"$M_j\cdot=10^{3}$")

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"datapoints $N$")
plt.ylabel(r"absolute relative change in energy")
plt.legend(loc="best", fontsize=12)

plt.savefig("../figures/conservation_energy2.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/conservation_energy2.pdf", "../figures/conservation_energy2.pdf"))
plt.show()
"""
#ANG MOM FOR DIFFERENT MJ
"""
exponents = np.linspace(3, 8, 21)
masses = ["1", "10", "100", "1000"]
planets = ["sun", "earth", "jupiter"]

change_in_energy_1    = np.zeros(len(exponents))
change_in_energy_10   = np.zeros(len(exponents))
change_in_energy_100  = np.zeros(len(exponents))
change_in_energy_1000 = np.zeros(len(exponents))

counter = 0

for i in exponents:
    energy = get_data("data/3_body/%5.4f_1_VV_energy.txt" % i, ["total", "kinetic", "potential", "ang_mom"])
    change_in_energy_1[counter] = np.mean(abs(np.array(energy["ang_mom"]) - np.array(energy["ang_mom"])[0])/abs(np.array(energy["ang_mom"])[0]))

    energy = get_data("data/3_body/%5.4f_10_VV_energy.txt" % i, ["total", "kinetic", "potential", "ang_mom"])
    change_in_energy_10[counter] = np.mean(abs(np.array(energy["ang_mom"]) - np.array(energy["ang_mom"])[0])/abs(np.array(energy["ang_mom"])[0]))

    energy = get_data("data/3_body/%5.4f_100_VV_energy.txt" % i, ["total", "kinetic", "potential", "ang_mom"])
    change_in_energy_100[counter] = np.mean(abs(np.array(energy["ang_mom"]) - np.array(energy["ang_mom"])[0])/abs(np.array(energy["ang_mom"])[0]))

    energy = get_data("data/3_body/%5.4f_1000_VV_energy.txt" % i, ["total", "kinetic", "potential", "ang_mom"])
    change_in_energy_1000[counter] = np.mean(abs(np.array(energy["ang_mom"]) - np.array(energy["ang_mom"])[0])/abs(np.array(energy["ang_mom"])[0]))

    counter += 1

plt.plot(np.power(10, exponents), change_in_energy_1,    "o", label=r"$M_j\cdot=10^{0}$")
plt.plot(np.power(10, exponents), change_in_energy_10,   "o", label=r"$M_j\cdot=10^{1}$")
plt.plot(np.power(10, exponents), change_in_energy_100,  "o", label=r"$M_j\cdot=10^{2}$")
plt.plot(np.power(10, exponents), change_in_energy_1000, "o", label=r"$M_j\cdot=10^{3}$")

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"datapoints $N$")
plt.ylabel(r"relative change in angular momentum")
plt.legend(loc="best", fontsize=12)

#plt.savefig("../figures/conservation_ang_mom_2.pdf", bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%("../figures/conservation_ang_mom_2.pdf", "../figures/conservation_ang_mom_2.pdf"))
plt.show()
"""
#ANG MOM FE VS VV
"""
exponents = np.linspace(3, 8.25, 22)
planets = ["sun", "earth", "jupiter"]

change_in_energy_FE = np.zeros(len(exponents))
change_in_energy_VV = np.zeros(len(exponents))

counter = 0
for i in exponents:
    energy = get_data("data/3_body/%5.4f_VV_energy.txt" % i, ["total", "kinetic", "potential", "ang_mom"])
    change_in_energy_VV[counter] = np.mean(abs(np.array(energy["ang_mom"]) - np.array(energy["ang_mom"])[0])/abs(np.array(energy["ang_mom"])[0]))

    energy = get_data("data/3_body/%5.4f_VV__FE_energy.txt" % i, ["total", "kinetic", "potential", "ang_mom"])
    change_in_energy_FE[counter] = np.mean(abs(np.array(energy["ang_mom"]) - np.array(energy["ang_mom"])[0])/abs(np.array(energy["ang_mom"])[0]))

    counter += 1

m, c, delta_m, delta_c = linear_regresion(exponents, np.log10(change_in_energy_VV))
plt.plot(np.power(10, exponents), 10**(c+m*exponents), label="VV - slope: %2.1f(%1.0f)" % (m, delta_m*10))
m, c, delta_m, delta_c = linear_regresion(exponents, np.log10(change_in_energy_FE))
plt.plot(np.power(10, exponents), 10**(c+m*exponents), label="FE - slope: %4.3f(%1.0f)" % (m, delta_m*1000))

plt.plot(np.power(10, exponents), change_in_energy_VV, "o")
plt.plot(np.power(10, exponents), change_in_energy_FE, "o")

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"datapoints $N$")
plt.ylabel(r"relative change in angular momentum")
plt.legend(loc="best", fontsize=12)
#plt.savefig("../figures/conservation_ang_mom.pdf", bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%("../figures/conservation_ang_mom.pdf", "../figures/conservation_ang_mom.pdf"))
plt.show()
"""
#CHANGE IN POSITION VV VS FE
"""
exponents = np.linspace(3.5, 8.5, 21)
planets = ["sun", "earth", "jupiter"]
change_in_position_VV = np.zeros(len(exponents))
change_in_position_FE = np.zeros(len(exponents))
counter = 0

for i in exponents:
    x = get_data("data/3_body/%5.4f_VV__FE__x.txt" % i, planets)
    y = get_data("data/3_body/%5.4f_VV__FE__y.txt" % i, planets)
    change_in_position_FE[counter] = np.sqrt(np.array(x["earth"])[-1]**2 + np.array(y["earth"])[-1]**2)

    x = get_data("data/3_body/%5.4f_VV__x.txt" % i, planets)
    y = get_data("data/3_body/%5.4f_VV__y.txt" % i, planets)
    change_in_position_VV[counter] = np.sqrt(np.array(x["earth"])[-1]**2 + np.array(y["earth"])[-1]**2)

    counter += 1

change_in_position_VV = abs(change_in_position_VV - np.sqrt(np.array(x["earth"])[-1]**2 + np.array(y["earth"])[-1]**2))
change_in_position_FE = abs(change_in_position_FE - np.sqrt(np.array(x["earth"])[-1]**2 + np.array(y["earth"])[-1]**2))

plt.plot(np.power(10, exponents), change_in_position_VV, "o", label="VV")
plt.plot(np.power(10, exponents), change_in_position_FE, "o", label="FE")

plt.legend(loc="best")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"datapoints $N$")
plt.ylabel(r"deviation final position [AU]")
plt.savefig("../figures/change_position.pdf", bbox_inches="tight")

os.system('pdfcrop %s %s &> /dev/null &'%("../figures/change_position.pdf", "../figures/change_position.pdf"))
plt.show()
"""
#PLOT CHANGE IN POSITION DIFFERENT MJ
"""
exponents = np.linspace(3, 8, 21)
planets = ["sun", "earth", "jupiter"]
change_in_position_1 = np.zeros(len(exponents))
change_in_position_10 = np.zeros(len(exponents))
change_in_position_100 = np.zeros(len(exponents))
change_in_position_1000 = np.zeros(len(exponents))

counter = 0

for i in exponents:
    x1 = get_data("data/3_body/%5.4f_1_VV__x.txt" % i, planets)
    y1 = get_data("data/3_body/%5.4f_1_VV__y.txt" % i, planets)
    change_in_position_1[counter] = np.sqrt(np.array(x1["earth"])[-1]**2 + np.array(y1["earth"])[-1]**2)

    x10 = get_data("data/3_body/%5.4f_10_VV__x.txt" % i, planets)
    y10 = get_data("data/3_body/%5.4f_10_VV__y.txt" % i, planets)
    change_in_position_10[counter] = np.sqrt(np.array(x10["earth"])[-1]**2 + np.array(y10["earth"])[-1]**2)

    x100 = get_data("data/3_body/%5.4f_100_VV__x.txt" % i, planets)
    y100 = get_data("data/3_body/%5.4f_100_VV__y.txt" % i, planets)
    change_in_position_100[counter] = np.sqrt(np.array(x100["earth"])[-1]**2 + np.array(y100["earth"])[-1]**2)

    x1000 = get_data("data/3_body/%5.4f_1000_VV__x.txt" % i, planets)
    y1000 = get_data("data/3_body/%5.4f_1000_VV__y.txt" % i, planets)
    change_in_position_1000[counter] = np.sqrt(np.array(x1000["earth"])[-1]**2 + np.array(y1000["earth"])[-1]**2)

    counter += 1

change_in_position_1 = abs(change_in_position_1 - np.sqrt(np.array(x1["earth"])[-1]**2 + np.array(y1["earth"])[-1]**2))
change_in_position_10 = abs(change_in_position_10 - np.sqrt(np.array(x10["earth"])[-1]**2 + np.array(y10["earth"])[-1]**2))
change_in_position_100 = abs(change_in_position_100 - np.sqrt(np.array(x100["earth"])[-1]**2 + np.array(y100["earth"])[-1]**2))
change_in_position_1000 = abs(change_in_position_1000 - np.sqrt(np.array(x1000["earth"])[-1]**2 + np.array(y1000["earth"])[-1]**2))

plt.plot(np.power(10, exponents), change_in_position_1, "o", label=r"$M_j\cdot=10^{0}$")
plt.plot(np.power(10, exponents), change_in_position_10, "o", label=r"$M_j\cdot=10^{1}$")
plt.plot(np.power(10, exponents), change_in_position_100, "o", label=r"$M_j\cdot=10^{2}$")
plt.plot(np.power(10, exponents), change_in_position_1000, "o", label=r"$M_j\cdot=10^{3}$")

plt.legend(loc="best", fontsize=12)
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"datapoints $N$")
plt.ylabel(r"deviation final position [AU]")
#plt.savefig("../figures/change_position_MJ.pdf", bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%("../figures/change_position_MJ.pdf", "../figures/change_position_MJ.pdf"))
plt.show()
"""
