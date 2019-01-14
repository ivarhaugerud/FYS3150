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

temp = [2.31]#np.arange(0.5, 4.00, 0.2)



results = get_data("../data/large_run_2.31_expect.txt", ["MC", "E", "EE", "M", "MM", "absM"])


plt.plot(results["E"])
plt.show()
