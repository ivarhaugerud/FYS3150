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


int main(int argc, char const *argv[])
  {
    int L = atoi(argv[1]);
    double exp_MC_cyc = atof(argv[2]);
    double T = atof(argv[3]);

    MainClass ising_model(L, "task_c", pow(10, 4));
    ising_model.initialize_random();
    ising_model.Run(T, pow(10, exp_MC_cyc));

    return 0;
  }
