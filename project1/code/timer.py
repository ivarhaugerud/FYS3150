import sys, os

"""
Program runs time_diff.cpp w times for n = 10**exponent. Remember to make the
executable time_diff before running this program. The results (time used for
the two algorithms) are written to a file named time_n.txt.
"""

try:
    exponent = float(sys.argv[1]) #exponent of n
    w = int(sys.argv[2])          #number of times to repeat experiment

except:
    print ("Usage of this script", sys.argv[0], "Integration points", sys.argv[1],\
           "Times repeating experiment", sys.argv[2])
    sys.exit(1)

n = int(10**exponent)


#Making file for resulting times
outfile = open('data/time_' + str(n) + '.txt', 'w')

# Define command line text string
cmdline  = './time_diff ' + str(exponent)
# Now run code, here c++ code  which has been compiled and linked
cmd = cmdline

for i in range(0, w):

    failure = os.system(cmd)
    if failure: #checks if the cpp code returned zero
        print ('running time_diff failed')
        sys.exit(1)

    infile = open('data/storefile.txt', 'r')

    for line in infile: # write the read data to a new file
        outfile.write('%g   %g\n' % (float(line.split()[0]), float(line.split()[1])))

#close files
outfile.close()
infile.close()
