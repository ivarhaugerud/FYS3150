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

//run code for many different beta valuess
int main(int argc, char const *argv[])
  {
    double alpha       = atof(argv[1]);
    long long int MCC  = pow(10, atoi(argv[2]));
    double omega       = atof(argv[3]);

    double step = pow(10,  0.143191866131  -0.5039370357278281*log10(alpha)); //best fit function for ideal step

    MainClass VMC("different_beta", 10000, omega, 2, alpha, 1, step, psi_T2, E_T2);
    //(file name, number of data lines, omega, nr_of_particles, alpha, beta, step)

    for (double beta = 0.05; beta < 1.5; beta += 0.05)
    {
      VMC.change_beta(beta);
      VMC.DeleteDatafile();
      VMC.Run(MCC);
      VMC.reset();
    }


    //VMC.DeleteDatafile();
    //VMC.Run(MCC);
    return 0;
  }
