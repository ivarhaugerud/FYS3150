#include "MainClass.hpp"
//all functions are identical to the ones in ../code/MainClass.cpp, see there for comments
//only some functions are removed
MainClass::MainClass()
{}

MainClass::MainClass(int size, string save_name, int amount_of_data)
{
  L = size;
  M = 0;
  E = 0;
  count_accepts = 0;

  spins = Mat<int>(L, L);
  index_vector = Col<int>(L+2);
  expectationvalues = vec(5, fill::zeros);
  boltzmann_factor  = vec(17, fill::zeros);
  filename = save_name;
  data_lines = amount_of_data;

  std::mt19937 generator (std::clock());
  std::uniform_real_distribution<double> zero_to_one_distribution(0.0, 1.0);
  std::uniform_int_distribution<int> zero_to_L_distribution(0, L-1);

  index_vector = Col<int>(L+2);
  for (int i = 1; i < L+2; i++)
  {
    index_vector(i) = i-1;
  }
  index_vector(0) = L-1;
  index_vector(L+1) = 0;
}

void MainClass::initialize_random()
  {
    for (int i = 0; i < L; i++)
    {
      for (int j = 0; j < L; j++)
      {
        spins(i, j) = (zero_to_one_distribution(generator) > 0.5) ? 1 : -1;
      }
    }
  }

double MainClass::calc_energy()
  {
    int energy = 0;        //in units of J

    //loops over everything except bottom row
    for (int i = 0; i < L; i++)
    {
      for (int j= 0; j<L; j++)
      {
        energy -= spins(i, j)*( spins( index_vector(i+1), index_vector(j) ) + spins(index_vector(i), index_vector(j+1)));
      }
    }
    return energy;
  }

double MainClass::calc_magnetization()
  {
    int mag = 0;

    for (int i = 0; i < L; i++)
    {
      for (int j = 0; j < L; j++)
      {
        mag += spins(i, j);
      }
    }
    return mag;
  }


void MainClass::Metropolis()
    {
      int delta_E;
      int x_i;
      int y_i;
      int neighbour_sum;

      for (int i = 0; i < L*L; i++)
      {
        //choose random spin
        x_i = zero_to_L_distribution(generator)%L+1; //we plus with one since we will use them on "index_vector", which is shifted by +1
        y_i = zero_to_L_distribution(generator)%L+1; //we plus with one since we will use them on "index_vector", which is shifted by +1

        neighbour_sum = spins( index_vector(y_i), index_vector(x_i+1) ) + spins( index_vector(y_i), index_vector(x_i-1) ) + spins( index_vector(y_i+1), index_vector(x_i) ) + spins( index_vector(y_i-1), index_vector(x_i) );

        delta_E = 2*neighbour_sum * spins( index_vector(y_i), index_vector(x_i) ); //J=1

        if (delta_E < 0 || zero_to_one_distribution(generator) <= boltzmann_factor(delta_E + 8))
        {
          E += delta_E;
          M -= 2*spins( index_vector(y_i), index_vector(x_i));
          spins( index_vector(y_i) , index_vector(x_i) ) *= -1;
        }
      }
    }

void MainClass::Run(double T, int nr_cycles)
{
    E = calc_energy();
    M = calc_magnetization();
    int denominator = nr_cycles/data_lines;

    stringstream stream;
    stream << fixed << setprecision(4) << T;

    ofstream outfile_expect("../final_data/" + filename + "_" + stream.str() + "_expect.txt");
    if (!outfile_expect.is_open())
      cout<<"Could not open file" << endl;

    int LL = L*L;
    int LLLL = LL*LL;

    for (int cycle = 0; cycle <= nr_cycles; cycle++)
    {

      Metropolis();

      expectationvalues(0) += E/LL;
      expectationvalues(1) += E*E/LLLL;
      expectationvalues(2) += M/LL;
      expectationvalues(3) += M*M/LLLL;
      expectationvalues(4) += fabs(M)/LL;

      if ( (cycle) % denominator == 0)
      {
        write_expectation_values(cycle+1, outfile_expect);
      }
    }
 }

void MainClass::equilibrate(double T, int nr_cycles)
{
  double beta = 1.0/(T); //k=J=1

  boltzmann_factor(-8+8) = exp(8.0*beta);
  boltzmann_factor(-4+8) = exp(4.0*beta);
  boltzmann_factor(0+8)  = exp(0);
  boltzmann_factor(4+8)  = exp(-4.0*beta);
  boltzmann_factor(8+8)  = exp(-8.0*beta);

  for (int i = 0; i < nr_cycles/10; i++)
    {
      Metropolis();
    }
  }

void MainClass::reset()
{
  expectationvalues(0) = 0;
  expectationvalues(1) = 0;
  expectationvalues(2) = 0;
  expectationvalues(3) = 0;
  expectationvalues(4) = 0;
}

void MainClass::write_expectation_values(int cycle_nr, ofstream& OutputFile)
  {
    vec dummy_expectationvalues = expectationvalues/(cycle_nr);
    OutputFile << cycle_nr << " " << dummy_expectationvalues(0) << " " << dummy_expectationvalues(1) << " " << dummy_expectationvalues(2) << " " << dummy_expectationvalues(3) << " " << dummy_expectationvalues(4) << "\n";
  }
