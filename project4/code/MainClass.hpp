#ifndef MAINCLASS
#define MAINCLASS

#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <string>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace arma;

class MainClass
{
public:
  int L;
  double E;
  double M;
  //save as double to avoid overflow

  //for writing to file
  int count_accepts;
  string filename;
  int data_lines;

  //vectors
  vec boltzmann_factor;
  vec expectationvalues;
  Mat<int> spins;
  Col<int> index_vector;

  //for timing
  double start;
  double stop;
  double time_used;

  //random number generators
  std::mt19937 generator;
  std::uniform_real_distribution<double> zero_to_one_distribution;
  std::uniform_int_distribution<int> zero_to_L_distribution;

  MainClass();
  MainClass(int size, string save_name, int amount_of_data);

  //the function names are self explanetory
  void initialize_random();
  void initialize_up();

  //calcualting initial energy and magnetization
  double calc_energy();
  double calc_magnetization();

  //run metropolis, two different ways, with and without counting accepts
  //metroplis is called in the Run function
  void Metropolis();
  void Metropolis_Count_Accept();
  void Run(double T, int nr_cycles);

  //time and equilibrate data, as well as reseting expectation values
  void equilibrate(double T, int nr_cycles);
  void equilibrate_magnetization(double T);
  void Timer(string start_or_stop, string filename);
  void reset();

  //write data to file
  void write_spinns(ofstream& OutputFile);
  void write_expectation_values(int cycle_nr, ofstream& inputFile);
};

#endif
