#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <random>
#include <ctime>
#include <string>
#include <typeinfo>

using namespace std;
using namespace arma;

Mat<int> initialize(int L)
  {
    std::mt19937 generator (std::time(NULL));
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    Mat<int> spins = Mat<int>(L, L);
    for (int i = 0; i < L; i++)
    {
      for (int j = 0; j < L; j++)
      {
        spins(i, j) = (dis(generator) > 0.5) ? 1 : -1;
      }
    }
    return spins;
  }

  Mat<int> initialize_up(int L)
    {
      std::mt19937 generator (std::time(NULL));
      std::uniform_real_distribution<double> dis(0.0, 1.0);

      cout << typeid(generator).name() << endl;

      Mat<int> spins = Mat<int>(L, L);
      for (int i = 0; i < L; i++)
      {
        for (int j = 0; j < L; j++)
        {
          spins(i, j) = 1;
        }
      }
      return spins;
    }

double calc_energy(Mat<int> matrix)
  {
    int energy = 0;        //in units of J
    int L = matrix.n_rows; //always symetric, i think
    Col<int> index_vector(L+2);//= vec(L+2);

    for (int i = 1; i < L+2; i++)
    {
      index_vector(i) = i-1;
    }
    index_vector(0) = L-1;
    index_vector(L+1) = 0;

    //loops over everything except bottom row
    for (int i = 0; i < L; i++)
    {
      for (int j= 0; j<L; j++)
      {
        energy -= matrix(i, j)*( matrix( index_vector(i+1), index_vector(j) ) + matrix(index_vector(i), index_vector(j+1)));
      }
    }

    return energy;
  }

double calc_magnetization(Mat<int> matrix)
  {
    int magnetization = 0;        //in units of J
    int L = matrix.n_rows;        //always symetric, i think

    for (int i = 0; i < L; i++)
    {
      for (int j = 0; j < L; j++)
      {
        magnetization += matrix(i, j);
      }
    }
    return magnetization;
  }

void Metropolis(Mat<int> & matrix, int & E, int & M, vec boltzmann_factor, Col<int> index_vector)
  {
    int L = matrix.n_rows;             //always symetric, i think
    int delta_E;
    int x_i;
    int y_i;
    int neighbour_sum;
    double r;

    std::mt19937 generator (std::time(NULL));
    std::uniform_int_distribution<int> dis1(0, L-1);
    std::uniform_real_distribution<double> dis2(0.0, 1.0);

    for (int i = 0; i < L*L; i++)
    {
      //choose random spinn
      x_i = dis1(generator)+1; //we plus with one since we will use them on "index_vector", which is shifted by +1
      y_i = dis1(generator)+1; //we plus with one since we will use them on "index_vector", which is shifted by +1

      neighbour_sum = matrix( index_vector(y_i), index_vector(x_i+1) ) + matrix( index_vector(y_i), index_vector(x_i-1)) + matrix( index_vector(y_i+1), index_vector(x_i) ) + matrix(index_vector(y_i-1), index_vector(x_i));

      r = dis2(generator);
      delta_E = -2*neighbour_sum * matrix( index_vector(y_i), index_vector(x_i) ); //J=1
      if (delta_E < 0 || r <= boltzmann_factor(delta_E + 8))
      {
        E += delta_E;
        M -= 2*matrix( index_vector(y_i), index_vector(x_i));

        matrix( index_vector(y_i) , index_vector(x_i) ) *= -1;

        /*
        cout << matrix << endl;
        cout << delta_E << endl;
        cout << boltzmann_factor(delta_E + 8) << "  " << r << endl;
        cout << endl;
        cout << endl;
        */
      }
    }
  }

void run_all(int L, int nr_cycles, double T_min, double T_max, double T_step, string filename)
  {
    //make datafile ready
    ofstream outfile("../data/" + filename + ".txt");
    if (!outfile.is_open())
      cout<<"Could not open file" << endl;

    //Mat<int> spinn = initialize(L);
    Mat<int> spinn = initialize_up(L);

    vec expectationvalues = vec(5);   //energy, energy^2, mangetization, mangetization^2, abs(mangetization)
    vec boltzmann_factor  = vec(17, fill::zeros);

    int E = calc_energy(spinn);
    int M = calc_magnetization(spinn);

    //create index vector
    Col<int> index_vector(L+2);
    for (int i = 1; i < L+2; i++)
      {
        index_vector(i) = i-1;
      }
    index_vector(0) = L-1;
    index_vector(L+1) = 0;

    for (double T = T_min; T < T_max; T+=T_step)
    {
      double beta = 1.0/(T);
      boltzmann_factor(-8+8) = exp(8.0*beta);
      boltzmann_factor(-4+8) = exp(4.0*beta);
      boltzmann_factor(0+8)  = exp(0);
      boltzmann_factor(4+8)  = exp(-4.0*beta);
      boltzmann_factor(8+8)  = exp(-8.0*beta);

      for (int cycle = 0; cycle <= nr_cycles; cycle++)
      {
        Metropolis(spinn, E, M, boltzmann_factor, index_vector);

        expectationvalues(0) += E;
        expectationvalues(1) += E*E;
        expectationvalues(2) += M;
        expectationvalues(3) += M*M;
        expectationvalues(4) += fabs(M);


        if ( cycle % 100000 ==0)
        {
          //
        }
        //outfile << E << "\n";
      }
      //expectationvalues /= nr_cycles*L*L;

      //outfile << T << " " << expectationvalues(0) << " " << expectationvalues(1) << " " << expectationvalues(2) << " " << expectationvalues(3) << " " << expectationvalues(4) << "\n";

      expectationvalues(0) = 0;
      expectationvalues(1) = 0;
      expectationvalues(2) = 0;
      expectationvalues(3) = 0;
      expectationvalues(4) = 0;

      spinn = initialize(L);
      E = calc_energy(spinn);
      M = calc_magnetization(spinn);
    }
  }

int main(int argc, char const *argv[])
  {

    run_all(10, pow(10, 6), 2.5, 4.5, 50, "positions");
    return 0;
  }
