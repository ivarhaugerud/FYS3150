# Project 4: Solving the 2D Ising model

This folder has 6 subfolders, [code](https://github.uio.no/cecilgl/FYS4150/tree/master/project4/code), [parallel_run](https://github.uio.no/cecilgl/FYS4150/tree/master/project4/parallel_run), [plotting](https://github.uio.no/cecilgl/FYS4150/tree/master/project4/plotting), [figures](https://github.uio.no/cecilgl/FYS4150/tree/master/project4/figures), [latex_ivar](https://github.uio.no/cecilgl/FYS4150/tree/master/project4/latex_ivar) and [latex_cecilie](https://github.uio.no/cecilgl/FYS4150/tree/master/project4/latex_cecilie). What the folder contains is quite self explanatory. A description of what the different files in each folder does is explained below. All of the code is commented to explain what each step does.

## Code :computer:
All of the code for this project was run with the following specs:

Computer  : 2,7 GHz Intel Core i5Processor <br />
Version   : MacOS Sierra 10.12.6 <br />
Processor : 2,7 GHz Intel Core i5 <br />
Memory    : 8 GB 1867 MHz DDR3 <br />
Graphics  : Intel Iris Graphics 6100 1536 MB <br />

We have a main class ```MainClass``` which does the Monte Carlo simulation for a given size of the lattice. In this class, we have a method ```Run```which runs the simulation with temperature and number of Monte Carlo cycles as arguments. The header file is ```MainClass.hpp``` and the functions are defined in ```MainClass.cpp```.

The program ```run_b.cpp``` does the simulations for a lattice of size L for temperatures from 0.5J/kB to 4J/kB. Command line arguments are first the lattice size then the logarithm of the number of Monte Carlo cycles.

The program ```run_c.cpp``` is used to do simulations for T = 1.0J/kB and T = 2.4J/kB.

The program ```run_timer.cpp``` runs a simulation without parallelization for a lattice size and logarithm of Monte Carlo cycles as command line arguments and caluclates the time the program uses to do the simuations for eigth different temperatures.

```problem_e.cpp``` is an unparallelized version for solving problem e, where Monte Carlo simulations are done for different temperatures and lattice sizes.

The data is not added to GitHub:octocat:, to reduce the amount of storage used by GitHub:octocat:, since some of the data files are quite large. The code we did not end up needing is saved in the subdirectory [code/trash_code](https://github.uio.no/cecilgl/FYS4150/tree/master/project4/code/trash_code).

## Parallelization

The code for running problem e parallelized is found in the directory [parallel_run](https://github.uio.no/cecilgl/FYS4150/tree/master/project4/parallel_run).

The program ```problem_e.cpp``` solves problem e using parallelization.

The program ```timer.cpp``` runs a simulation with parallelization for a lattice size and logarithm of Monte Carlo cycles as command line arguments and caluclates the time the program uses to do the simuations for eigth different temperatures.

When using MPI, remember to write ```module load mpi``` before making and ```CXX=mpicxx make``` instead of ```make``` and ``` mpirun ./program``` to run.

## Plotting

We have chosen to have all the code for plotting in a separate directory, [plotting](https://github.uio.no/cecilgl/FYS4150/tree/master/project4/plotting). Here, we have one plotting program for each problem where plotting is needed.

## Figures :chart_with_downwards_trend:
In the directory [figures](https://github.uio.no/cecilgl/FYS4150/tree/master/project4/figures) we save . Most of the figures are not pushed to GitHub:octocat: to save space.

## LaTeX :memo:
In the directories [latex_ivar](https://github.uio.no/cecilgl/FYS4150/tree/master/project4/latex_ivar) and [latex_cecilie](https://github.uio.no/cecilgl/FYS4150/tree/master/project4/latex_cecilie) we have all the code used for LaTeX. For our references we use BibTeX.
