# Project 1: Solving the 1D Poisson equation

This folder has 3 subfolders, [code](https://github.uio.no/cecilgl/FYS4150/tree/master/project1/code), [latex](https://github.uio.no/cecilgl/FYS4150/tree/master/project1/latex) and [figures](https://github.uio.no/cecilgl/FYS4150/tree/master/project1/figures). What the folder contains is quite self explanatory. A description of what the different files in each folder does is explained below.

## Code
All of the code for this project was run with the following specs:

Computer  : 2,7 GHz Intel Core i5Processor <br />
Version   : MacOS Sierra 10.12.6 <br />
Processor : 2,7 GHz Intel Core i5 <br />
Memory    : 8 GB 1867 MHz DDR3 <br />
Graphics  : Intel Iris Graphics 6100 1536 MB <br />

The main file of this project is the file ```generate_data.cpp```. This file contains a function for writing results to file, named ```results2file```, a function giving the right hand side of our differential equation,```f ``` and a function giving the analytic solution of the problem, named ```u_analytic```. In this code we run the ```4n``` algorithm, since it is the fastest. The data it writes to file is the x-values, u'', the absolute difference between analytical and numerical solution, and the logarithm of the relative deviation between numerical and analytical solution.

The ```main``` function allocates memory dynamically. A number n must be provided as a command line argument when running the executable of ```generate_data.cpp```, where n+2 is the number of points the problem is solved for.

The file ```generate_data.cpp``` includes the headerfile ```rref_special.cpp``` which includes the functions for the special Gaussian elimination.

The file ```timer.py``` runs the 4n and 8n algorithm and reads the time used by each function, and writes the data to the file. In this file we can choose the number for ```n```, and the number of runs for each run using sys.argv in python. This function uses the C++ file ```time_diff.cpp```.

The file ```timer_LU.py``` runs LU decomposition algorithm and reads the time used by the algorithm, and writes the data to a file. In this file we can choose the number for ```n```, and the number of runs for each run using sys.argv in python. This function uses the C++ file ```time_diff_LU.cpp```.

The file ```test_rref_tri.cpp``` checks if our algorithm is implemented correctly, and still works, by importing the Gaussian elimination function, and comparing with a know calculated results, for a 4x4 matrix.

The files ```rref_tri.hpp``` and ```rref_special.hpp``` are the functions we import in the other programs to do the Gaussian elimination. The two files are respectivly the 8n and the 4n algorithm.

The file ```plot_data.py``` plots the data, which is stored in the folder /data, and save the figures in the folder ../figures.

The file ```LU.cpp``` is used to generate data using the LU decomposition. This is only  used to test the results from our algorithm to the results we find using LU decomposition.

The two files ```lib.h``` and ```lib.hpp``` are not written by us, and is just used to conduct the LU decomposition.

The data is not added to GitHub, to reduce the amount of storage used by GitHub, since some of the data files are quite large. The data is saved in the subdirectory ```data```. The code we did not end up needing is saved in the subdirectory ```trash_code```


## Figures
In the directory ```figures``` we save the figures which is plotted by the program ```plot_data.py``` in the directory ```code```. Most of the figures are not pushed to GitHub to save space.

## LaTeX
In this directory we have all the code used for LaTeX. For our references we use BibTeX.
