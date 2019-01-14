#include <iostream>
#include <cmath>
#include <armadillo>
//importing libaries

//Our own functions:
#include "functions.hpp"

//make printing better
using namespace std;
using namespace arma;

/*
*/

//function for writing results to file
void results2file(vec analytic, vec numeric, int len_eig)
  {
    ofstream outfile("data/temporary_data.txt");
    if (!outfile.is_open())
      cout<<"Could not open file" << endl;
    else
    {
      for (int i = 0; i < len_eig; i++)
      {
      outfile << abs(analytic(i)-numeric(i+1)) << endl;
      }
    }
  }

/*
*/

/*
*/

//This program works together with data_3D_plot.py

int main(int argc, char *argv[])
{
  //use dynamic memory allocation for the value of n, send in expoent of 10
  int n          = atoi(argv[1]);
  double rho_max = atof(argv[2]);
  double epsilon = pow(10, -atof(argv[3]));
  double rho_nul = 0;
  int len_eig = 5;
  int counter = 0;
  double h = (rho_max-rho_nul)/(n); // stepsize
  double hh = h*h;                  // reduce FLOPS

  vec d   = vec(n+1);       // diagonal elements
  vec e   = vec(n+1);       // secondary-diagonal elements
  vec rho = vec(n+1);       // position variable
  vec analytical_eigenvalues = vec(len_eig);

  for (int t = 0; t < len_eig; t++)
  {
    analytical_eigenvalues(t) = 3+t*4;
  }

  rho(n) = 0;
  d(n) = 0;
  e(n) = 0;


  for (int i = 1; i < n; i++)
  {
    rho(i) = i*h; //by definition
    d(i)   = 2.0/hh + rho(i)*rho(i);  //by definition
    e(i)   = -1.0/hh; //by definition
  }


  mat matrix = create_tridiagonal_vector(n, e, d);              //create matrix
  vec jac_meth = sort(jacobys_method(matrix, epsilon, counter));//create eigenvalues
  results2file(analytical_eigenvalues, jac_meth, len_eig);      //write the results to file, this file will be read by the python program data_3D_plot.py

  return 0;
}
