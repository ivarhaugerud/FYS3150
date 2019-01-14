#ifndef MAINCLASS
#define MAINCLASS

//import what we need
#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <string>
#include <iomanip>
#include <fstream>

//make it easier
using namespace std;
using namespace arma;

class MainClass
{
//the public variables
public:

  //writing to file
  string filename;
  int data_lines;

  //scalers
  double omega;
  double prob_at_R;
  double E;
  int nr_accepted;
  double step;
  double alpha;
  double beta;

  //functions
  double (*psi)(mat, double, double, double);
  double (*calcE)(mat, double, double, double);

  //matrices and vectors
  mat R;
  mat R_trial;
  vec expectationvalues;
  vec average_expectations;

  //random number generators
  std::mt19937 generator;
  std::uniform_real_distribution<double> zero_to_one_distribution;

  //the class
  MainClass();
  MainClass(string save_name, int amount_of_data, double OMEGA,
            int nr_particles, double ALPHA, double BETA, double STEP,
            double wavefunc(mat, double, double, double),
            double energyfunc(mat, double, double, double));

  //running and equilibrating
  void Metropolis();
  void Run(int nr_cycles);
  void Equilibrate(int equilibration_cycles);

  //calculate potential energy for virial theorem
  double calc_potential(mat r);

  //change varaiable values
  void  change_step(double  new_step);
  void  change_alpha(double new_alpha);
  void  change_beta(double  new_beta);
  void  change_omega(double new_omega);

  // reset and writing/deleteing files
  void reset();
  void write_expectation_values(int cycle_nr, ofstream& OutputFile);
  void write_accepts(int cycle_nr, ofstream& OutputFile);
  void write(ofstream& OutputFile);
  void DeleteDatafile();
};

#endif
