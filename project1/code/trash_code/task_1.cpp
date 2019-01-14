#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
//Including ibaries

//Our own methods:
#include "rref_special.hpp"

using namespace std;

/*
*/

//Function for writing the data to file
void results2file(string filename, int n, double x[], double f[], double analytic[])
  {
    ofstream outfile("data/" + filename + to_string(n) + ".txt");
    if (!outfile.is_open())
      cout<<"Could not open file" << endl;
    else
    {
      for (int i = 0; i < n+2; i++)
        outfile << setw(15) << x[i] << setw(15) << f[i] << setw(15) << abs(f[i]-analytic[i]) << setw(15) << log10(abs((f[i]-analytic[i])/(analytic[i]))) << endl;
    }
  }

/*
*/

//Make inline functions to make the code run faster
inline double func(double var){ return 100*exp(-10*var);}
inline double u_analytic(double var){ return 1-(1-exp(-10))*var-exp(-10*var);}
/*
*/

int main(int argc, char *argv[])
  {
  // send in the size of the arrays as the exponent to make life easier
  int n = pow(10, atof(argv[1]));

  //defining arrays, make them n+2 to include the end points
  double *x = new double[n+2];
  double *u = new double[n+2];
  double *f = new double[n+2];
  double *analytic = new double[n+2];

  //these elements will stay constant due to the boundary conditions, and are therfor defined here
  u[0] = 0;
  u[n+1] = 0;

  double h = 1.0/(n+1); //from definition
  double hh = h*h; //reduce number of FLOPS

  for (int i = 0; i < n+2; i++)
    {
    x[i] = h*i;
    f[i] = func(x[i])*hh;
    analytic[i] = u_analytic(x[i]);
    //arrays we need for the algorithm, and writing to file
    }

  //More arrays we need for our algorithm, and solution
  double *dtilde = new double[n+1];
  double *ftilde = new double[n+1];

  //these elements will not get indexed in the following for-Loop, do it here instead
  dtilde[n] = (n+1.0)/n;

  for (int i = 1; i < n; i++)
  {
    dtilde[i] = (i+1.0)/i;
    //define all the known array elements outside the algorithm
  }

  rref_special(u, f, n, dtilde, ftilde); //runs algorithml
  results2file("data_", n, x, u, analytic); //saves the data

  //delte the arrays
  delete [] x; delete [] u; delete [] f;
  delete [] analytic; delete [] dtilde; delete [] ftilde;

  return 0;
  }
