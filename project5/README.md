# Project 5: Variational Monte Carlo of electrons in a quantum dot

This folder has 4 subfolders, [code](https://github.uio.no/cecilgl/FYS4150/tree/master/project5/code), [plotting](https://github.uio.no/cecilgl/FYS4150/tree/master/project5/plotting), [figures](https://github.uio.no/cecilgl/FYS4150/tree/master/project5/figures), [latex](https://github.uio.no/cecilgl/FYS4150/tree/master/project5/latex). What the folder contains is quite self explanatory. A description of what the different files in each folder does is explained below. All of the code is commented to explain what each step does.

## Code :computer:
Most of the code for this project was run with the following specs:

Computer  : 2,7 GHz Intel Core i5Processor <br />
Version   : MacOS Sierra 10.12.6 <br />
Processor : 2,7 GHz Intel Core i5 <br />
Memory    : 8 GB 1867 MHz DDR3 <br />
Graphics  : Intel Iris Graphics 6100 1536 MB <br />

The parallelized code was run on a 64-core computer.

We have a main class ```MainClass``` which does the Monte Carlo simulation for a given trial wave function and corresponding local energy. In this class, we have a method ```Run```which runs the simulation with number of Monte Carlo cycles as arguments.

The program ```run.cpp``` does the simulation with the optimal step length. Command line arguments are alpha, beta, the exponent of the number of MC cycles and omega. The trial wave function is currently set to psiT2.

```different_alpha.cpp``` is used to do the MC simulations for different values of alpha, from alpha = 0.6 to alpha = 1.6 with steps of alpha on 0.2. The trial wave function is currently set to pisT1.

```different_beta.cpp``` is used to do the MC simulations for different values of beta with the trial wave function psiT2, from beta = 0.05 to beta = 1.5 with steps of beta on 0.05. This program is not used when we study the optimal combination of alpha and beta, as we use a parallelized program, see parallelization. 

```different_omega.cpp``` is used to do the MC simulations for different values of omega from 0.01 to 1 with a step of 0.01.

```different_step.cpp``` does the VMC simulation for different alphas and different step lengths and is used to find the optimal step length.

The file ```WFandE.hpp``` contains our trial wave functions and local energies.

The data is not added to GitHub:octocat:, to reduce the amount of storage used by GitHub:octocat:, since some of the data files are quite large. The code we did not end up needing is saved in the subdirectory [code/trash_code](https://github.uio.no/cecilgl/FYS4150/tree/master/project5/code/trash_code).


## Parallelization

The code for finding the energy for different alpha's and beta's is parallelized using MPI. The program is called ```parallel.cpp```.

When using MPI, remember to write ```module load mpi``` before making and ```CXX=mpicxx make``` instead of ```make``` and ``` mpirun ./program``` to run.

## Plotting

We have chosen to have all the code for plotting in a separate directory, [plotting](https://github.uio.no/cecilgl/FYS4150/tree/master/project5/plotting). Here, we have one plotting program for each problem where plotting is needed. All data is analyzed in these programs.

## Figures :chart_with_downwards_trend:
In the directory [figures](https://github.uio.no/cecilgl/FYS4150/tree/master/project5/figures) we save all the figures produced by the code in [plotting](https://github.uio.no/cecilgl/FYS4150/tree/master/project5/plotting).

## LaTeX :memo:
In the directory [latex](https://github.uio.no/cecilgl/FYS4150/tree/master/project5/latex) we have all the code used for LaTeX. For our references we use BibTeX.
