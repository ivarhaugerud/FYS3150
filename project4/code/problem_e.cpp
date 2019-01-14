#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <random>
#include <ctime>
#include <string>

#include "MainClass.hpp"

using namespace std;
using namespace arma;

//used to solve problem e, this was used before we got the parallelization to work

int main(int argc, char const *argv[])
{
  int exp_MC_cyc = 6;

  for (int L = 20; L < 110; L += 20)
  {
    MainClass ising_model(L, to_string(L), 1000);
    ising_model.initialize_random();
    //initialize for each lattice size, random configuration each time

    for (double T = 1.95; T <= 2.10; T += 0.01)
    {
      ising_model.equilibrate(T, pow(10, exp_MC_cyc))
      ising_model.Run(T, pow(10, exp_MC_cyc));
      ising_model.reset();
    }
  }

  return 0;
}
