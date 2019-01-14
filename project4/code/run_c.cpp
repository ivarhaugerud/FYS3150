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

//used to solve problem c

int main(int argc, char const *argv[])
{
  int L = atoi(argv[1]);
  double exp_MC_cyc = atof(argv[2]);
  //from command line

  double T_1 = 1;
  double T_2 = 2.4;
  //temperatures we are interested in

  MainClass ising_unorder(L, "c_unorder", pow(10, exp_MC_cyc-2));
  //initialize random start

  //run for T_1
  ising_unorder.initialize_random();
  ising_unorder.equilibrate(T_1, exp_MC_cyc);
  ising_unorder.Run(T_1, pow(10, exp_MC_cyc));
  ising_unorder.reset();

  //run for T_2
  ising_unorder.initialize_random();
  ising_unorder.equilibrate(T_2, exp_MC_cyc);
  ising_unorder.Run(T_2, pow(10, exp_MC_cyc));
  ising_unorder.reset();


  MainClass ising_order(L, "c_order", pow(10, exp_MC_cyc-1));
  //initialze orderd start

  //run for T_1
  ising_order.initialize_up();
  ising_order.Run(T_1, pow(10, exp_MC_cyc));
  ising_order.reset();

  //run for T_2
  ising_order.initialize_up();
  ising_order.Run(T_2, pow(10, exp_MC_cyc));
  ising_order.reset();

  return 0;
}
