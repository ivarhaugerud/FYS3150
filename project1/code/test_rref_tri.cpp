#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
//importing libaries

//Our own methods:
#include "rref_tri.hpp"

//makes printing better
using namespace std;


/*
Doing a test for a 4x4 non-symmetric tridiagonal matrix:
[1  2  0  0
 5  3  6  0
 0  1  2  3
 0  0  4  1]

 so a = [5  1  4]
    b = [1  3  2  1]
    c = [2  6  3]

We can i.e. use f = [1  2  3  4]

Gaussian elimination by hand then gives u = [-13/8  21/16  33/32  -1/8]
*/

int main()
{

  /*
  Remembering that our indexation in the rref_tri starts at 1.
  */
  int n = 4;

  //defining arrays
  double a[] = {0, 5, 1, 4};
  double b[] = {0, 1, 3, 2, 1};
  double c[] = {0, 2, 6, 3};
  double f[] = {0, 1, 2, 3, 4};
  double *btilde = new double[n+1];
  double *ftilde = new double[n+1];
  double u[n + 1];

  rref_tri(a, b, c, u, f, n, btilde, ftilde); //runs the algorithm

  //analytic solution
  double sol[] = {0, -13.0/8, 21.0/16, 33.0/32, -1.0/8};

  double tol = 1e-13; //define a tolerance when testing the numerical solution to analytical

  for (int i = 1; i < n + 1; i++)
  //loops over all values and test if the are correct
  {
    if (abs(u[i]- sol[i]) > tol)
    {
      cout << "WARNING for i = " << i << " ; u = " << u[i] << endl;
    }
    else
    {
      cout << "i = " << i << " works fine." << endl;
    }
  }

  /*
  Everything seems to work fine.
  */

  return 0;
}
