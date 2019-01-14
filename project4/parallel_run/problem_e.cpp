#include <mpi.h>
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


int main(int argc, char *argv[])
{
  //variabels used
  int my_rank, numprocs;
  int exp_MC_cyc = 6;

  //  MPI initializations
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  //temperature interval used for one run
  double Tmin = 2.2090;
  double Tmax = 2.3265;
  double deltaT = 0.0025;
  int nr_data_lines = 1000;
  /*
  Determine number of interval which are used by all processes
  myloop_begin gives the starting point on process my_rank
  myloop_end gives the end point for summation on process my_rank
  */
  int no_temps = ((Tmax - Tmin)/deltaT+1)/numprocs;
  double myloop_begin = my_rank*no_temps*deltaT + Tmin;
  double myloop_end = (my_rank+1)*no_temps*deltaT + Tmin - deltaT;
  if ( (my_rank == numprocs-1) && ( myloop_end < Tmax) ) myloop_end = Tmax;

  //initialize all classes
  MainClass ising_model_40(40,   to_string(40),  nr_data_lines);
  MainClass ising_model_60(60,   to_string(60),  nr_data_lines);
  MainClass ising_model_80(80,   to_string(80),  nr_data_lines);
  MainClass ising_model_100(200, to_string(200), nr_data_lines);

  //give random configuration
  ising_model_40.initialize_random();
  ising_model_60.initialize_random();
  ising_model_80.initialize_random();
  ising_model_100.initialize_random();

  //write that it has begun
  ising_model_40.equilibrate(myloop_begin,  pow(10, exp_MC_cyc));
  cout << "Rank: " << my_rank << " done equlibriating 40." << endl;
  ising_model_60.equilibrate(myloop_begin,  pow(10, exp_MC_cyc));
  cout << "Rank: " << my_rank << " done equlibriating 60." << endl;
  ising_model_80.equilibrate(myloop_begin,  pow(10, exp_MC_cyc));
  cout << "Rank: " << my_rank << " done equlibriating 80." << endl;
  ising_model_100.equilibrate(myloop_begin, pow(10, exp_MC_cyc));
  cout << "Rank: " << my_rank << " done equlibriating 100." << endl;

for (double T = myloop_begin; T <= myloop_end + 0.000001; T += deltaT)
 {
  cout << T << "  " << my_rank <<  endl; //write how far it has come

  //calibrate
  ising_model_40.equilibrate(T,  pow(10, exp_MC_cyc));
  ising_model_60.equilibrate(T,  pow(10, exp_MC_cyc));
  ising_model_80.equilibrate(T,  pow(10, exp_MC_cyc));
  ising_model_100.equilibrate(T, pow(10, exp_MC_cyc));

  //run and reset for all latice sizes
  ising_model_40.Run(T, pow(10, exp_MC_cyc));
  ising_model_40.reset();

  ising_model_60.Run(T, pow(10, exp_MC_cyc));
  ising_model_60.reset();

  ising_model_80.Run(T, pow(10, exp_MC_cyc));
  ising_model_80.reset();

  ising_model_100.Run(T, pow(10, exp_MC_cyc));
  ising_model_100.reset();
 }

  // End MPI
  MPI_Finalize ();

  return 0;
}
