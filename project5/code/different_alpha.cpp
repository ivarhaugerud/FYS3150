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

//code for running with many different values of alpha
int main(int argc, char const *argv[])
  {
    long long int MCC  = pow(10, atoi(argv[1]));
    double alpha = 1.00;
    double beta = 1.00;
    int number_of_particles = 2;
    double omega = 1;
    double step = 1;

    MainClass VMC("final_data_2", 1000, omega, number_of_particles, alpha, beta, step, psi_T1, E_T1);
    //(file name, number of data lines, omega, nr_of_particles, alpha, beta, step)

    for (double alpha = 0.6; alpha < 1.6+0.01; alpha += 0.2)
    {
      VMC.change_alpha(alpha);
      step = pow(10,  0.143191866131  -0.5039370357278281*log10(alpha)); //best fit function for ideal step

      VMC.DeleteDatafile();
      //VMC.Equilibrate(pow(10, 5));
      VMC.Run(MCC);
      VMC.reset();
    }
  return 0;
  }
