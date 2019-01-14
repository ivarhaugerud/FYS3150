import os, sys
import numpy as np

"""
Program to run a C++ program for finding the ideal value of rho_max
This is done by looping over different values for n and rho and comparing
them to the analytical eigenvalues for the 5 lowest energy staes
"""

q_min = 2.5 #lowest  value of rho_max
q_max = 7.5 #largest value of rho_max
n_min = 100 #minimum n value
n_max = 200 #maximum n value
number_of_q = 20                #number of iterations
q = np.zeros(number_of_q + 6)   #add some extra elements
q[:-6] = np.linspace(q_min, q_max, number_of_q)

#add some extra values
q[-6] = 4.45
q[-5] = 4.50
q[-4] = 4.55
q[-3] = 4.60
q[-2] = 4.65
q[-1] = 4.75

n = np.linspace(n_min, n_max, number_of_q+6) #the n values we will use
epsilon = 5 #the value we will use to epsilon

#Making file for resulting times
outfile = open('data/data.txt', 'w') #the data file we will read from

#loop over the n values
for k in range(len(n)):
    outfile = open('data/data_%4.3f.txt' % n[k], 'w') #current text file to write to

    #loop over q values
    for i in range(len(q)):
        cmdline  = './program ' + str(int(n[k])) + " " + str(q[i]) + " " + str(epsilon) #what to write in command line
        outfile.write(str(q[i]))
        print(cmdline) #to check how far the code has ran

        failure = os.system(cmdline)

        if failure: #checks if the cpp code returned zero
            print ('running time_diff failed')
            sys.exit(1)

        infile = open('data/temporary_data.txt', 'r')

        for line in infile: # write the read data to a new file
            outfile.write('   %g    ' % (float(line.split()[0])))
        outfile.write("\n")

    #close files
    outfile.close()
    infile.close()
