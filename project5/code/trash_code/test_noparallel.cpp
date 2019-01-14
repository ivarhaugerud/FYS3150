#include <mpi.h>
#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <random>
#include <ctime>
#include <string>
#include <typeinfo>

#include "MainClass.hpp"
#include "WFandE.hpp"

using namespace std;
using namespace arma;

int main(int argc, char *argv[])
  {
    //variabels used
    int my_rank, numprocs;

    numprocs = 1;
    my_rank = 0;


    long long int MCC  = pow(10, 6);
    double omega       = 1;

    long double step;
    double alpha, beta;

    //Alpha interval
    double Amin = 0.5;
    double Amax = 1.5;
    double deltaA = 0.01;
    /*
    Determine number of interval which are used by all processes
    myloop_begin gives the starting point on process my_rank
    myloop_end gives the end point for summation on process my_rank
    */
    int no_alphas = ((Amax - Amin)/deltaA+1)/numprocs;
    double myloop_begin = my_rank*no_alphas*deltaA + Amin;
    double myloop_end = (my_rank+1)*no_alphas*deltaA + Amin - deltaA;
    if ( (my_rank == numprocs-1) && ( myloop_end < Amax) ) myloop_end = Amax;

    alpha = 0.5;
    beta = 0.5;
    step = pow(10,  0.143191866131 - 0.5039370357278281*log10(alpha));

    MainClass VMC("noparallel", 1000, omega, 2, alpha, beta, step, psi_T2, E_T2); //(file name, number of data lines, omega, nr_of_particles, alpha, beta, step)

    for (alpha = myloop_begin; alpha <= myloop_end + deltaA/2.0; alpha += deltaA)
    {
      step = pow(10,  0.143191866131 - 0.5039370357278281*log10(alpha));
      VMC.change_step(step);
      VMC.change_alpha(alpha);
      for (beta = 0.05; beta < 1.5; beta += 0.01)
      {
        VMC.change_beta(beta);
        VMC.DeleteDatafile();
        VMC.Run(MCC);
        VMC.reset();
      }
    }

    //VMC.DeleteDatafile();
    //VMC.Run(MCC);

    return 0;
  }
