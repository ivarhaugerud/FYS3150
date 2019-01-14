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

//times the code with parallelization and writes results to file

int main(int argc, char const *argv[])
  {
    double exp_MC_cyc = atof(argv[1]);
    int L = atoi(argv[2]);
    //what we can change from the command line

    double T = 1.0;

    MainClass ising_model(L, "timer", pow(10, exp_MC_cyc));
    //initialize

    ising_model.initialize_random(); //random start
    ising_model.Timer("start", "timer");  //starts timer

    //begins calculating
    for (double T = 1.0; T < 2.75+0.001; T += 0.25)
    {
      ising_model.Run(T, pow(10, exp_MC_cyc));
    }
    ising_model.Timer("stop", "timer"); //stops timer and writes time to file

    return 0;
  }
