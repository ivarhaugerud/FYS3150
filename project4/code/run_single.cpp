#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <random>
#include <ctime>
#include <string>
#include <typeinfo>

#include "MainClass.hpp"

using namespace std;
using namespace arma;

//runs for a single configuration which can be specified in command line
//used for example for making plot of spins with black down and white up

int main(int argc, char const *argv[])
  {
    int L = atoi(argv[1]);
    double exp_MC_cyc = atof(argv[2]);
    double T = atof(argv[3]);
    //variabels we can change from command line

    MainClass ising_model(L, "test", 4); //(size, file name, number of data lines)
    ising_model.initialize_random();                //initialize
    ising_model.equilibrate(T, pow(10, exp_MC_cyc));//equilibrate
    ising_model.Run(T, pow(10, exp_MC_cyc));        //run

    return 0;
  }
