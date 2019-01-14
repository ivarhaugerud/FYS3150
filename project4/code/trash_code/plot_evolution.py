import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def get_data(filename, variables):
    df = pd.read_csv(filename,\
                    delim_whitespace=True, \
                    engine='python', \
                    names=variables)
    return df
    #using pandas to read the data files

positions = np.loadtxt("../data/test_0.50_spins_.txt")
L = int(np.sqrt(len(positions[0, :])))

plt.subplot(111)
plt.axis('off')
for i in range(len(positions[:, 0])):
    spinns = np.reshape(positions[i, :], (L, L))

    plt.imshow(spinns, cmap='Greys',  interpolation='nearest')
    #plt.draw()
    #plt.pause(0.000000001)
plt.show()
