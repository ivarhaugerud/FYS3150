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

//simple run function for general input
int main(int argc, char const *argv[])
  {
    double alpha       = atof(argv[1]);
    double beta        = atof(argv[2]);
    long long int MCC  = pow(10, atof(argv[3]));
    double omega       = atof(argv[4]);
    double step = pow(10,  0.143191866131  -0.5039370357278281*log10(alpha)); //best fit function for ideal step

    int number_of_particles = 2;
    cout << step << endl;

    MainClass VMC("omegas_E2", 100000, omega, number_of_particles, alpha, beta, step, psi_T2, E_T2);
    //           (file name,   number of data lines, omega, nr_of_particles, alpha, beta, step, psi, E)

    VMC.DeleteDatafile();
    VMC.Equilibrate( pow(10, 5) );
    VMC.Run(MCC);

    return 0;
  }
