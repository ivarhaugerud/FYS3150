#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <random>
#include <ctime>
#include <string>
#include <typeinfo>

#include "MainClass.hpp" //the class
#include "wfandE.hpp"    //the functions for energy and prob-density

using namespace std;
using namespace arma;

//run for different step sizes for different alpha values
int main(int argc, char const *argv[])
  {
    long long int MCC  = pow(10, atoi(argv[1]));
    double omega = 1;
    int number_of_particles = 2;
    double alpha = 1.0;
    double beta = 1.0;
    double step = 1.0;

    MainClass VMC("find_step", 100, omega, number_of_particles, alpha, beta, step, psi_T1, E_T1);
    //           (file name,   number of data lines, omega, nr_of_particles, alpha, beta, step, psi, E)

    for (double alpha = 0.3; alpha < 2.5; alpha += 0.1)
    {
      VMC.change_alpha(alpha);
      VMC.DeleteDatafile();

      for (double step = 0.1; step < 3; step += 0.1)
      {
          VMC.change_step(step);
          VMC.Run(MCC);
          VMC.reset();
      }
    }
  return 0;
  }
