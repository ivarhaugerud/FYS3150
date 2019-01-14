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
  double T_1 = 1;
  double T_2 = 2.4;

  MainClass ising_unorder(L, "c_unorder", pow(10, exp_MC_cyc-2));

  ising_unorder.initialize_random();
  ising_unorder.equilibrate(T_1, exp_MC_cyc);
  ising_unorder.Run(T_1, pow(10, exp_MC_cyc));
  ising_unorder.reset();

  ising_unorder.initialize_random();
  ising_unorder.equilibrate(T_2, exp_MC_cyc);
  ising_unorder.Run(T_2, pow(10, exp_MC_cyc));
  ising_unorder.reset();

  /*
  MainClass ising_order(L, "c_order", pow(10, exp_MC_cyc-1));

  ising_order.initialize_up();
  ising_order.Run(T_1, pow(10, exp_MC_cyc));
  ising_order.reset();
  ising_order.initialize_up();
  ising_order.Run(T_2, pow(10, exp_MC_cyc));
  ising_order.reset();
  */

  return 0;
}
