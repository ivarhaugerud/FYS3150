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

//function for writing data to file
void results2file(string filename, int n, int k, double times_1[], double times_2[])
  {
    ofstream outfile("data/" + filename + to_string(n) + ".txt");
    if (!outfile.is_open())
      cout<<"Could not open file" << endl;
    else
    {
      for (int i = 0; i < k; i++)
        outfile << setw(15) << times_1[i] << setw(15) << times_2[i] << endl;
    }
  }

/*
*/

//analytical expression for f
double func(double var)
  {
    return 100*exp(-10*var);
  }

/*
*/

//analytical expression for the double derivative of f
  double u_analytic(double var)
    {
      return 1-(1-exp(-10))*var-exp(-10*var);
    }
/*
*/

int main(int argc, char *argv[])
  {
  int n = pow(10, atof(argv[1])); //exponent of 10
  int w = atoi(argv[2]);          //number of runs

  //define arrays
  double *x = new double[n+2];
  double *u = new double[n+2];
  double *f = new double[n+2];
  double *times_1 = new double[w];
  double *times_2 = new double[w];

  //the boundary conditions are set
  u[0] = 0;
  u[n+1] = 0;

  double h = 1.0/(n+1); //from definition
  double hh = h*h;      //reduce number of FLOPS
  for (int i = 0; i < n+2; i++)
    {
    x[i] = h*i;
    f[i] = func(x[i])*hh;
    }

  //define new arrays
  double *a = new double[n];
  double *b = new double[n+1];
  double *c = new double[n];

  double *dtilde = new double[n+1];
  double *btilde = new double[n+1];
  double *ftilde = new double[n+1];

  //these are outside the for-loop, and is therefor defined here
  b[n] = 2;
  dtilde[n] = (n+1.0)/n;

  for (int i = 1; i < n; i++)
  {
    a[i] = -1;
    b[i] = 2;
    c[i] = -1;

    dtilde[i] = (i+1.0)/i; //the analytical expression for dtilde
  }

  //defining the time variables as doubles
  double start1;
  double finish1;
  double start2;
  double finish2;

  //calculates the time
  for (int k = 0; k < w; k++)
  {
    start1 = clock();
    rref_tri(a, b, c, u, f, n, btilde, ftilde);
    finish1 = clock();

    start2 = clock();
    rref_special(u, f, n, dtilde, ftilde);
    finish2 = clock();

    times_1[k] = (finish1-start1)/CLOCKS_PER_SEC;
    times_2[k] = (finish2-start2)/CLOCKS_PER_SEC;
  }
  results2file("time_", n, w, times_1, times_2);


  delete [] x;
  delete [] u;
  delete [] f;

  delete [] a;
  delete [] b;
  delete [] c;

  delete [] btilde;
  delete [] dtilde;
  delete [] ftilde;

  return 0;
  }
