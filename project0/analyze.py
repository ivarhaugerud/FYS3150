import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.style.use("bmh")
sns.color_palette("hls", 1)


infile = open("data_cec.txt")
h = []
diff_2_flt = []
diff_2_dbl = []
diff_3_flt = []
diff_3_dbl = []
for line in infile:
    h.append(float(line.split()[0]))
    diff_2_flt.append(float(line.split()[1]))
    diff_2_dbl.append(float(line.split()[2]))
    diff_3_flt.append(float(line.split()[3]))
    diff_3_dbl.append(float(line.split()[4]))


plt.figure(1)
plt.loglog(h, diff_2_flt, 'r', label="single precision")
plt.loglog(h, diff_2_dbl, 'g', label="double precision")
plt.legend(loc="best")
plt.loglog(h, diff_2_flt, 'r.')
plt.loglog(h, diff_2_dbl, 'g*')
plt.xlabel("h")
plt.ylabel("Deviation between numerical and analytic")

plt.title("Method diff2c")
plt.savefig("method1.pdf")

plt.figure(2)
plt.loglog(h, diff_3_flt, 'b', label="single precision")
plt.loglog(h, diff_3_dbl, 'k', label="double precision")
plt.legend(loc="best")
plt.loglog(h, diff_3_flt, 'b.')
plt.loglog(h, diff_3_dbl, 'k*')
plt.xlabel("h")
plt.ylabel("Deviation between numerical and analytic")

plt.title("Method diff3c")
plt.savefig("method2.pdf")
