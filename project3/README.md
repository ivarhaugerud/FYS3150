
# Project 3: The Solar System :waning_gibbous_moon: :earth_africa: :sunny: :sparkles:

This folder has 3 subfolders, [code](https://github.uio.no/cecilgl/FYS4150/tree/master/project3/code), [figures](https://github.uio.no/cecilgl/FYS4150/tree/master/project3/figures) and [latex](https://github.uio.no/cecilgl/FYS4150/tree/master/project3/latex). What the folder contains is quite self explanatory. A description of what the different files in each folder does is explained below.

The simulations use an included data file containing the positions and velocities of a wide range of objects in the Solar System. The positions and velocities are dated 2018-Sep-28 00:00:00.0000 TDB and extracted from [NASAs JPL HORIZONS](http://ssd.jpl.nasa.gov/horizons.cgi#top).

## Code :computer:
All of the code for this project was run with the following specs:

Computer  : 2,7 GHz Intel Core i5Processor <br />
Version   : MacOS Sierra 10.12.6 <br />
Processor : 2,7 GHz Intel Core i5 <br />
Memory    : 8 GB 1867 MHz DDR3 <br />
Graphics  : Intel Iris Graphics 6100 1536 MB <br />


The data is not added to GitHub:octocat:, to reduce the amount of storage used by GitHub:octocat:, since some of the data files are quite large. The data is saved in the subdirectory [code/data](https://github.uio.no/cecilgl/FYS4150/tree/master/project3/code/data). The code we did not end up needing is saved in the subdirectory [code/trash_code](https://github.uio.no/cecilgl/FYS4150/tree/master/project3/code/trash_code).

We have one main class called SolarSystemSimulator. This class is used to simulate the orbit of any objects in the Solar System around the Solar Systems. The header file is called ```SolarSystemSimulator.hpp``` and the corresponding functions can be found in ```SolarSystemSimulator.cpp```.

The subclass SunAtRest simulates the Solar System with the Sun at rest in the origin. The header file is called ```SunAtRest.hpp``` and the corresponding functions can be found in ```SunAtRest.cpp```.

The subclass RelativisticCorrection adds a relativistic correction to Newtons gravitational law. This class is specialized to the case of two objects, and is used for studying the perihelion precession of Mercury. The header file is called ```RelativisticCorrection.hpp``` and the corresponding functions can be found in ```RelativisticCorrection.cpp```.

The subclass ChangeExponent gives us the possibility of changing the exponent of Newton's gravitational law. The headerfile is called ```ChangeExponent.hpp```, and the corresponding functions can be found in ```ChangeExponent.cpp```. The file ```different_exponent.cpp``` is using this subclass to calculate the position for different values of the exponent.

The file ```stability_3body.cpp``` increases the mass of Jupiter and runs the simulation for different number of datapoints.

The file ```main.cpp``` is the general file which we can run for different command line arguments to produce all the results we could want. The arguments are: max time[years], exponent of the number of data points, sun at rest[true/false], relativity[true/false], VV/FE, save energy files[true/false], name of planets.
Example for running a program with relativistic corrections for 100 years and 10^10 data points for the Sun and Mercury using Velocity Verlet:
```./program 100 10 false true VV false sun mercury```.

To read information about the masses, initial position, and initial velocities we use ```json.hpp``` and the data file ```planet_data.json```

The file ```find_escape_velocity.cpp``` runs the simulation for many different escape velocities.

The file ```plot_data.py``` is where we plot all of the data stored in the folder ```data```, and store the figures in the directory [figures](https://github.uio.no/cecilgl/FYS4150/tree/master/project3/figures).

## Figures :chart_with_downwards_trend:
In the directory [figures](https://github.uio.no/cecilgl/FYS4150/tree/master/project3/figures) we save the figures which is plotted by the program ```plot_data.py``` in the directory [code](https://github.uio.no/cecilgl/FYS4150/tree/master/project3/code).

## LaTeX :memo:
In the directory [latex](https://github.uio.no/cecilgl/FYS4150/tree/master/project3/latex) we have all the code used for LaTeX. For our references we use BibTeX.
