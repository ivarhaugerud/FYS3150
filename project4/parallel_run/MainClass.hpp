#ifndef MAINCLASS
#define MAINCLASS

//all functions are identical to the ones in ../code/MainClass.hpp, see there for comments
//only some functions are removed

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
  int count_accepts;
  string filename;

  double start;
  double stop;
  double time_used;

  int data_lines;
  vec boltzmann_factor;
  vec expectationvalues;
  Mat<int> spins;
  Col<int> index_vector;

  std::mt19937 generator;
  std::uniform_real_distribution<double> zero_to_one_distribution;
  std::uniform_int_distribution<int> zero_to_L_distribution;

  MainClass();
  MainClass(int size, string save_name, int amount_of_data);

  void initialize_random();
  double calc_energy();
  double calc_magnetization();
  void Metropolis();

  void write_expectation_values(int cycle_nr, ofstream& inputFile);
  void equilibrate(double T, int nr_cycles);
  void reset();
  void Run(double T, int nr_cycles);
};

#endif
