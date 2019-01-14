import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

def get_data(filename, variables):
    df = pd.read_csv(filename,\
                    delim_whitespace=True, \
                    engine='python', \
                    names=variables)
    return df
    #using pandas to read the data files

#positions = np.loadtxt("../final_data/medium_long_spin_1.000_spins.txt")
positions = np.loadtxt("../final_data/tester_1.000_spins.txt")

L = int(np.sqrt(len(positions[0, :])))

counter = 0
plt.figure(figsize=(20, 8))
print(len(positions[:, 0]))
print(np.shape(positions))

for i in [0, 2, 3, 4]:
    plt.subplot(2,5,1+counter)
    spinns = np.reshape(positions[i, :], (L, L))
    plt.axis('off')
    plt.imshow(spinns, cmap='Greys',  interpolation='nearest')

    counter += 1

positions = np.loadtxt("../final_data/large_long_spin_1.000_spins.txt")

for i in [-2]:
    plt.subplot(2,5,1+counter)
    spinns = np.reshape(positions[i, :], (L, L))
    plt.axis('off')
    plt.imshow(spinns, cmap='Greys',  interpolation='nearest')
    counter += 1

positions = np.loadtxt("../final_data/task_d_2.400_spins_2.txt")
L = int(np.sqrt(len(positions[0, :])))

for i in range(int(len(positions[:, 0])-1)):
    plt.subplot(2,5,1+counter)
    spinns = np.reshape(positions[i, :], (L, L))
    plt.axis('off')
    plt.imshow(spinns, cmap='Greys',  interpolation='nearest')

    counter += 1

plt.tight_layout(pad=0.1, w_pad=0, h_pad=0)
plt.savefig("../figures/evolution.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/evolution.pdf", "../figures/evolution.pdf"))
plt.show()
