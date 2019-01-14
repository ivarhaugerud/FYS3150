import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
plt.style.use("bmh")
sns.color_palette("hls", 1)

def get_data(filename, start, end, variables):
    df = pd.read_csv(filename,\
                    #skiprows=int(start),\
                    #skipfooter=int(end),\
                    delim_whitespace=True, \
                    engine='python', \
                    names=variables)
    return df

data = abs(get_data('data.txt', 1, 26, ['h', 'd_1_single', 'd_2_single', 'd_1_double', 'd_2_double', "richardsson"]))

plt.figure(1)
plt.loglog(data['h'], data['d_1_single'], 'r', label="method 1, single")
plt.loglog(data['h'], data['d_1_double'], 'g', label="method 1, double")
plt.loglog(data['h'], data['d_2_single'], 'b', label="method 2, single")
plt.loglog(data['h'], data['d_2_double'], 'k', label="method 2, double")
plt.loglog(data['h'], data['richardsson'], 'c', label="richardsson, double")


plt.legend(loc="best")
plt.loglog(data['h'], data['d_1_single'], 'r.')
plt.loglog(data['h'], data['d_1_double'], 'g*')
plt.loglog(data['h'], data['richardsson'], 'c*')

#plt.figure(2)

plt.loglog(data['h'], data['d_2_single'], 'b*')
plt.loglog(data['h'], data['d_2_double'], 'k.')

plt.xlabel("h")
plt.ylabel("Deviation between numerical and analytic")

plt.show()
