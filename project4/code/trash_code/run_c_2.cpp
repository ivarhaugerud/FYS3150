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

    MainClass ising_model(L, "task_c_2", pow(10, 1));
    ising_model.initialize_random();


    for (double T = 4.3; T < 4.800001; T += 0.1)
    {
      ising_model.Run(T, pow(10, exp_MC_cyc));
      ising_model.reset();
    }
    return 0;
  }
