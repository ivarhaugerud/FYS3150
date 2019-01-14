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

    MainClass ising_model(L, "task_b", 10000);
    ising_model.initialize_random();

    for (double T = 0.5; T < 4.0; T += 0.01)
    {
      ising_model.equilibrate(T, exp_MC_cyc);
      ising_model.Run(T, pow(10, exp_MC_cyc));
      ising_model.reset();
    }
    return 0;
  }
