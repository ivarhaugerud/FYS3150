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

//this program solves task b - generates data for many temperatures

int main(int argc, char const *argv[])
  {
    //can change these variabels easily
    int L = atoi(argv[1]);
    double exp_MC_cyc = atof(argv[2]);

    MainClass ising_model(L, "task_b", 10000); //(size, filename, number of times to write to file)
    ising_model.initialize_random();           //intialize random

    for (double T = 0.5; T < 4.0; T += 0.05)
    {
      ising_model.equilibrate(T, exp_MC_cyc);   //equilibrate
      ising_model.Run(T, pow(10, exp_MC_cyc));  //run
      ising_model.reset();                      //reset expectation values (not spin positions, as these are used)
    }
    return 0;
  }
