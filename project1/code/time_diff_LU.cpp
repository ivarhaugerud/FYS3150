#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <time.h>
//importing libaries

//importing our own functions
#include "lib.hpp"

//make priting easier
using namespace std;

//function for defining the A^ matrix
// also used in LU.cpp
double** tridiagonal(int number)
{
  double** matrix;
  matrix = new double*[number];
  for(int i = 0; i < number; i++)
    matrix[i] = new double[number];

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

// analytical expression for f
double func(double var)
{
  return 100*exp(-10*var);
}

//Function for writing time used by the algorithms to file
void results2file(string filename, double time_used)
{
  ofstream outfile("data/" + filename);
  if (!outfile.is_open())
    cout<<"Could not open file" << endl;
  else
  {
    outfile << time_used << endl;
  }
}

int main(int argc, char *argv[])
{
  int n = pow(10, atof(argv[1]));
  //send in the exponent of 10, for n

  //define matrixs and arrays
  double **A = tridiagonal(n);
  double *b = new double[n];
  int *indx = new int[n];
  double *d = new double[n];
  double *x = new double[n];

  double h = 1.0/(n+1); //definition of h
  double hh = h*h; //reduces number of flops
  for (int i = 0; i < n; i++)
    {
    x[i] = h*(i+1);
    b[i] = func(x[i])*hh;
    }

  double start = clock(); //start clock
  ludcmp(A, n, indx, d);  //the two steps for LU decomposition, first step
  lubksb(A, n, indx, b);  // second step
  double finish = clock();//finish timer
  double time_used = (finish-start)/CLOCKS_PER_SEC; //print time used
  results2file("storefile_LU.txt", time_used);      //store data in text file

  return 0;
}
