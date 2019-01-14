#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
//Importing libaries

#include "lib.hpp"
//lib.hpp contains the LU decomposition algorithm, which we need

using namespace std;
//makes printing easier

//Create the tri-diagonal matrix with 2's along the diagonal, and 1 along the second largest diagonals
double** tridiagonal(int number)
{
  double** matrix;
  matrix = new double*[number];
  for(int i = 0; i < number; i++)
    matrix[i] = new double[number];
  //Create a double pointer for our nxn matrix

  //Loops over all elements in the matrix, and check if we want to make the element 1 or 2 dependent on the difference between 1 and 2
  //If not we fill the element with zero

  for (int i = 0; i < number; i++)
    for (int j = 0; j < number; j++)
    {
      if (i == j)
      {
        matrix[i][j] = 2;
      }
      else if (i == j+1)
      {
        matrix[i][j] = -1;
      }
      else if (i == j-1)
      {
        matrix[i][j] = -1;
      }
      else
      {
        matrix[i][j] = 0;
      }
    }
  return matrix;
}

// The function which we want to find the numerical solution for
double func(double var)
{
  return 100*exp(-10*var);
}

//Function which writes the data to file
void results2file(string filename, int argc, double x[], double f[])
{
  ofstream outfile("data/" + filename + to_string(argc) + ".txt");
  if (!outfile.is_open())
    cout<<"Could not open file" << endl;
  else
  {
    for (int i = 0; i < argc; i++)
      outfile << setw(15) << x[i] << setw(15) << f[i] << endl;
  }
}

int main(int argc, char *argv[])
{
  //Define the size of our matrix and arrays by dynaic memory allocation
  int n = atoi(argv[1]);

  double **A = tridiagonal(n); //Creates the matrix

  //Defining the vectors we need, using double for all
  double *b = new double[n];
  int *indx = new int[n];
  double *d = new double[n];
  double *x = new double[n];

  double h = 1.0/(n+1); //definition of h
  double hh = h*h; //to reduce #FLOPS
  for (int i = 0; i < n; i++)
    {
    x[i] = h*(i+1);       //x values
    b[i] = func(x[i])*hh; //calc f_star
    }

  ludcmp(A, n, indx, d);  //LU decomposition takes two steps, step 1
  lubksb(A, n, indx, b);  //Step 2

  results2file("LU_", n, x, b); //save data to file, so we can compare the results from our algoithm with the results from LU decomposition

  return 0;
}
