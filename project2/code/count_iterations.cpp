#include <iostream>
#include <cmath>
#include <armadillo>

//Code to count the number of iterations as a function of the size of the matrix and write the results to file

#include "functions.hpp" //headerfile over functions we want to use

using namespace arma;
using namespace std;

int main(int argc, char *argv[])
{
  //dynamical memory allocation for variabels we want to alter
  int n_low = atoi(argv[1]);          //the size of the smallest matrix we want to use
  int n_hig = atoi(argv[2]);          //the size of the largest  matrix we want to use
  double step_factor = atof(argv[3]); //the step factor, use this instaed off adding a constant term since this looks better in a log-log plot

  double a = 1;   //constant values along diagonal and secondary diagonal
  double d = 2;   //constant values along diagonal and secondary diagonal
  double epsilon = pow(10, -7); //tolerance
  int number_of_runs = 0;       //define number of iterations needed
  mat matrix;                   //define a matrix to be calculated, size of this matrix will change each run

  //open a text file to write the data to
  ofstream outfile("data/iterations_jacoby.txt");
      if (!outfile.is_open())
      {
        cout<<"Could not open file" << endl;
      }

  //loop over until n > n_hig
  for (int i = n_low; i <= n_hig;)
  {
    matrix = create_tridiagonal(i, a, d);                             //create matrix
    vec jac_meth = jacobys_method(matrix, epsilon, number_of_runs);   //calc eigenvalues, and find number of iterations

    outfile << i << "   " << number_of_runs << endl;    //write the size of the matrix and the number of iterations to the text file
    i *= step_factor; //multiply by the stepfactor
    number_of_runs = 0; //set the number of iterations to zero
  }
  return 0;
}
