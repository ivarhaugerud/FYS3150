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

//run for different values of omega
int main(int argc, char const *argv[])
  {
    double alpha       = atof(argv[1]);
    double beta        = atof(argv[2]);
    long long int MCC  = pow(10, atof(argv[3]));
    double omega       = atof(argv[4]);
    double step = pow(10,  0.143191866131  -0.5039370357278281*log10(alpha)); //best fit function for ideal step

    int number_of_particles = 2;


    MainClass VMC("virial_E0", 100000, omega, number_of_particles, alpha, beta, step, psi_T1, E_T0);
    //           (file name,   number of data lines, omega, nr_of_particles, alpha, beta, step, psi, E)

    for (double omega = 0.01; omega < 1+0.001; omega += 0.01)
    {
      VMC.change_omega(omega);
      VMC.Equilibrate( pow(10, 5) );
      VMC.Run(MCC);
      VMC.reset();
    }

    return 0;
  }
