#include <mpi.h>
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

//run for both alpha and beta in parallel
int main(int argc, char *argv[])
{
  //variabels used
  int my_rank, numprocs;

  //  MPI initializations
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  cout << my_rank << "/" << numprocs << endl;

  long long int MCC  = pow(10, 8);
  double omega       = 1;

  double step, alpha, beta;

  //Alpha interval
  double Amin = 0.78;
  double Amax = 1.415;
  double deltaA = 0.005;
  /*
  Determine number of interval which are used by all processes
  myloop_begin gives the starting point on process my_rank
  myloop_end gives the end point for summation on process my_rank
  */
  int no_alphas = ((Amax - Amin)/deltaA+1)/numprocs;
  no_alphas = 4;
  double myloop_begin = my_rank*no_alphas*deltaA + Amin;
  double myloop_end = (my_rank+1)*no_alphas*deltaA + Amin - deltaA;
  if ( (my_rank == numprocs-1) && ( myloop_end < Amax) ) myloop_end = Amax;

  alpha = 0.1;
  beta = 0.1;
  step = pow(10,  0.143191866131 - 0.5039370357278281*log10(alpha));

  MainClass VMC("parallel", 1000, omega, 2, alpha, beta, step, psi_T2, E_T2);
  //(file name, number of data lines, omega, nr_of_particles, alpha, beta, step)

  for (alpha = myloop_begin; alpha <= myloop_end + deltaA/2.0; alpha += deltaA)
  {
    step = pow(10,  0.143191866131  -0.5039370357278281*log10(alpha)); //best fit function for ideal step
    VMC.change_step(step);
    VMC.change_alpha(alpha);
    VMC.DeleteDatafile();
    for (beta = 0; beta <= 0.9; beta += 0.005)
    {
      VMC.change_beta(beta);
      VMC.Equilibrate(pow(10,5));
      VMC.Run(MCC);
      VMC.reset();
      cout << "Rank " << my_rank << ": " << alpha << "  " << beta << endl;
    }
  }
  // End MPI
  MPI_Finalize ();

  return 0;
}
