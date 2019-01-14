#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <time.h>
//importing libaries

//Our own methods:
#include "rref_tri.hpp"
#include "rref_special.hpp"

//make printing better
using namespace std;

/*
*/

//function for writing results to file
void results2file(double time_1, double time_2)
  {
    ofstream outfile("data/storefile.txt");
    if (!outfile.is_open())
      cout<<"Could not open file" << endl;
    else
    {
      outfile << setw(15) << time_1 << setw(15) << time_2 << endl;
    }
  }

/*
*/

//make inline functions to reduce computation time
inline double func(double var){return 100*exp(-10*var);}
inline double u_analytic(double var){return 1-(1-exp(-10))*var-exp(-10*var);}

/*
*/

int main(int argc, char *argv[])
{
  //use dynamic memory allocation for the value of n, send in expoent of 10
  int n = pow(10, atof(argv[1]));

  //define arrays
  double *x = new double[n+2];
  double *u = new double[n+2];
  double *f = new double[n+2];
  double *analytic = new double[n+2];

  //since these are the boundary conditions, we define them here
  u[0] = 0;
  u[n+1] = 0;

  double h = 1.0/(n+1); //from definition
  double hh = h*h; //reduce number of FLOATS
  for (int i = 0; i < n+2; i++)
    {
    //calculates the analytical solution
    x[i] = h*i;
    f[i] = func(x[i])*hh;
    analytic[i] = u_analytic(x[i]);
    }

  //define more arrays
  double *a = new double[n];
  double *b = new double[n+1];
  double *c = new double[n];

  double *dtilde = new double[n+1];
  double *btilde = new double[n+1];
  double *ftilde = new double[n+1];

  //these will not be included in the for-loop, will therfor define them here
  b[n] = 2;
  dtilde[n] = (n+1.0)/n;

  for (int i = 1; i < n; i++)
  {
    //loops over, and define, known values
    a[i] = -1;
    b[i] = 2;
    c[i] = -1;

    dtilde[i] = (i+1.0)/i; //for the special algorithm we have an analytical expression for dtilde, which is calculated here
  }

  //Timing algorithm

  //time of the general algorithm
  double start1 = clock();
  rref_tri(a, b, c, u, f, n, btilde, ftilde);
  double finish1 = clock();

  //time for the special algorithm
  double start2 = clock();
  rref_special(u, f, n, dtilde, ftilde);
  double finish2 = clock();

  //calculate how many seconds used by each algorithm
  double time1 = (finish1-start1)/CLOCKS_PER_SEC;
  double time2 = (finish2-start2)/CLOCKS_PER_SEC;

  //write this to file
  results2file(time1, time2);


  //Delete arrays
  delete [] x; delete [] u; delete [] f;
  delete [] a; delete [] b; delete [] c;
  delete [] btilde; delete [] dtilde; delete [] ftilde;

  return 0;
}
