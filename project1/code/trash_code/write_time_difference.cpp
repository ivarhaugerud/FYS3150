#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <time.h>

//Our own methods:
#include "rref_tri.hpp"
#include "rref_special.hpp"

using namespace std;

/*
*/

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

double func(double var)
  {
    return 100*exp(-10*var);
  }

/*
*/

  double u_analytic(double var)
    {
      return 1-(1-exp(-10))*var-exp(-10*var);
    }
/*
*/

int main(int argc, char *argv[])
  {
  int n = pow(10, atof(argv[1]));
  int w = atoi(argv[2]);

  double *x = new double[n+2];
  double *u = new double[n+2];
  double *f = new double[n+2];
  double *times_1 = new double[w];
  double *times_2 = new double[w];

  u[0] = 0;
  u[n+1] = 0;

  double h = 1.0/(n+1);
  double hh = h*h;
  for (int i = 0; i < n+2; i++)
    {
    x[i] = h*i;
    f[i] = func(x[i])*hh;
    }

  double *a = new double[n];
  double *b = new double[n+1];
  double *c = new double[n];

  double *dtilde = new double[n+1];
  double *btilde = new double[n+1];
  double *ftilde = new double[n+1];

  b[n] = 2;
  dtilde[n] = (n+1.0)/n;

  for (int i = 1; i < n; i++)
  {
    a[i] = -1;
    b[i] = 2;
    c[i] = -1;

    dtilde[i] = (i+1.0)/i;
  }
  double start1;
  double finish1;
  double start2;
  double finish2;

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
