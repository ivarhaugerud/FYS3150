# Project 2: Eigenvalue problems in pyhsics: Solving the Schrödinger equation

This folder has 3 subfolders, [code](https://github.uio.no/cecilgl/FYS4150/tree/master/project2/code), [figures](https://github.uio.no/cecilgl/FYS4150/tree/master/project2/figures) and [latex](https://github.uio.no/cecilgl/FYS4150/tree/master/project2/latex). What the folder contains is quite self explanatory. A description of what the different files in each folder does is explained below. All of the code is commented to explain what each step does.

## Code :computer:
All of the code for this project was run with the following specs:

Computer  : 2,7 GHz Intel Core i5Processor <br />
Version   : MacOS Sierra 10.12.6 <br />
Processor : 2,7 GHz Intel Core i5 <br />
Memory    : 8 GB 1867 MHz DDR3 <br />
Graphics  : Intel Iris Graphics 6100 1536 MB <br />


The data is not added to GitHub:octocat:, to reduce the amount of storage used by GitHub:octocat:, since some of the data files are quite large. The data is saved in the subdirectory ```data```. The code we did not end up needing is saved in the subdirectory [code/trash_code](https://github.uio.no/cecilgl/FYS4150/tree/master/project2/code/trash_code).

The file ```functions.cpp``` include all of the functions we use in this project. The functions does the following:
write results to file, calculate the largest non diagonal element, rotate in jacobis method, find the analytical eigenvalues for the test function, create a tridiagonal matrix with scalar input, create tridiagonal matrix with vector input, to compute the whole of jacobis method, and to find the eigenvectors of the input matrix using jacobis method.
These functions are included in other programs through the headerfile ```functions.hpp```.

The file ```count_iterations.cpp``` runs jacobis method for different matrix sizes and calculated the iterations needed to reach a given accuracy, and writes the results to file. The file takes three input arguments; the exponent of the lowest size matrix, the exponent of the highest size matrix, and the step factor for each run.

The file ```size_for_leading_digits.cpp``` calculates the needed matrix size to reach a specific accuracy in the lowest eigenvalues. It takes three command line arguments: the start size of the matrix, the maximal position value, and tolerance.

The file ```test_all.cpp``` runs all of the four test functions we have implemented for our system. The code takes no command line arguments. It tests the following: are the analytical eigenvalues the same as the numerical ones using the armadillo eigenvalue finder and the algorithm using jacobis method which we have implemented, tests eigenvectors for a 3x3 matrix versus analytical solution, test jacobis method eigenvalues vs analytical, and lastly test if the off function chooses the largest non-diag element.

The file ```2_d.cpp``` find the eigenvalues and eigenvectors, and writes them to file, for the single electron in a 3D harmonic oscillator potential for the lowest energy eigenstates.

The file ```2_e.cpp``` find the eigenvalues and eigenvectors, and writes them to file, for the two interacting electrons in a 3D harmonic oscillator potential for the lowest energy eigenstates.

The data generated in C++ is plotted in python through the two files ```plot_data.py``` and ```plot_3D_data.py```. The 3D plot for finding ideal value of rho_max is made in ```plot_3D_data.py```, while all other plots are made from ```plot_data.py```. The date for the 3D plot is calculated by ```try_3D_plot.cpp```, which takes three inputarguments from the command line: the size of the matrix, rho_max and the tolerance. This program is run for different command line arguments through the python file ```data_3D_plot.py```.

## Figures :chart_with_downwards_trend:
In the directory [figures](https://github.uio.no/cecilgl/FYS4150/tree/master/project2/figures) we save the figures which is plotted by the programs ```plot_data.py``` and ```plot_3D_data.py``` in the directory [code](https://github.uio.no/cecilgl/FYS4150/tree/master/project2/code). Most of the figures are not pushed to GitHub:octocat: to save space.

## LaTeX :memo:
In the directory [latex](https://github.uio.no/cecilgl/FYS4150/tree/master/project2/latex) we have all the code used for LaTeX. For our references we use BibTeX.
