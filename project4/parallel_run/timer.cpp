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
  //variabels
  int my_rank, numprocs;
  int exp_MC_cyc = 5;
  int L = atoi(argv[1]);

  //variabels for timer
  double start;
  double stop;
  double time_used;

  //intializing the MPI
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  //  MPI initializations

  //the temperature values we will use
  double Tmin = 1;
  double Tmax = 2.75;
  double deltaT = 0.25;
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

   //initialize class
   MainClass ising_model(L, to_string(L),  nr_data_lines);
   ising_model.initialize_random(); //random configuration

   start = clock(); //start clock

   //start looping over temperature values
   for (double T = myloop_begin; T <= myloop_end + 0.000001; T += deltaT)
     {
      ising_model.Run(T, pow(10, exp_MC_cyc));
     }

  MPI_Finalize (); //finalize the MPI

  stop = clock();  //finish timer
  time_used = (stop-start)/CLOCKS_PER_SEC; //calculate time used

  //write time used to file
  std::ofstream outfile;
  outfile.open("../final_data/timer_" + to_string(my_rank) + ".txt", std::ios_base::app);
  if (!outfile.is_open())
    cout << "Could not open file" << endl;
  else
    outfile << L << "      " << my_rank << "    " << time_used << endl;
  return 0;
}
