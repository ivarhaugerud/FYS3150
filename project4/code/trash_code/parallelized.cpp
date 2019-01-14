#include "mpi.h"
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


int main(int argc, char const *argv[])
{
  int n_spins, mcs, my_rank, numprocs;
  //  MPI initializations
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  if (my_rank == 0 && argc <= 1) {
    cout << "Bad Usage: " << argv[0] <<
      " read output file" << endl;
    exit(1);
  }

  double exp_MC_cyc = 8;
  double mcs = pow(10, 8)

  int no_intervalls = mcs/numprocs;
  int myloop_begin = my_rank*no_intervalls + 1;
  int myloop_end = (my_rank+1)*no_intervalls;
  if ( (my_rank == numprocs-1) &&( myloop_end < mcs) ) myloop_end = mcs;

  double initial_temp = 2.1;
  double final_temp = 2.3;
  double temp_step = 0.02;

  for (double L = 40; L <= 100; L += 20)
    double n_spins = L*L


    // broadcast to all nodes common variables
    MPI_Bcast (&n_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&initial_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MainClass ising_model(L, "test", 1000);

    ising_model.initialize_random();
    //ising_model.Run(2.4, pow(10, 5));

    double  TimeStart, TimeEnd, TotalTime;
    TimeStart = MPI_Wtime();

    for (double T = initial_temp; T < final_temp; T += temp_step)
    {
      cout << T << endl;
      ising_model.Run(T, pow(10, exp_MC_cyc));
      ising_model.reset();
    }

  TimeEnd = MPI_Wtime();
    TotalTime = TimeEnd-TimeStart;
    if ( my_rank == 0) {
      cout << "Time = " <<  TotalTime  << " on number of processors: "  << numprocs  << endl;
    }

  // End MPI
  MPI_Finalize();

  return 0;
}
