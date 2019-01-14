#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>

#include "functions.hpp" //importing functions

using namespace std;
using namespace arma;

// Functions to calculate the needed size of the matrix to achieve three leading digits

int main(int argc, char *argv[])
{
  //dynamic memory allocation
  int n = atoi(argv[1]);
  double rho_max = atof(argv[2]);
  double epsilon = pow(10, -atof(argv[3]));

  double tolerance = pow(10, -4);
  int counter = 0;
  int counter_loop = 0;

  //the number of eigenvalues we are interested in requiering to be less than the tolerance
  int number_of_eigen_values = 4;
  bool achived;


  //loop over all of the eigenvalues
  for (int l = 0; l < number_of_eigen_values; l++)
  {
    achived = false;

    while (!achived) //will continue to increase the size of the matrix by 10 until this is satisfied
    {
      n = n + 10;//pow(10, atof(argv[1])+0.05*counter_loop);
      cout << n << endl;

      double h = rho_max/n; // stepsize
      double hh = h*h;                  // reduce FLOPS

      vec d   = vec(n+1);       // diagonal elements
      vec e   = vec(n+1);       // secondary-diagonal elements
      vec rho = vec(n+1);       // position variable

      rho(n) = 0;
      d(n) = 0;
      e(n) = 0;

      for (int i = 1; i < n; i++)
      {
        rho(i) = i*h;                   //by definition
        d(i)   = 2.0/hh + rho(i)*rho(i);//diagonal elements by definition
        e(i)   = -1.0/hh;               //secondary diagonal elements
      }

      mat matrix = create_tridiagonal_vector(n, e, d);  //create matrix
      vec jac_meth = sort(jacobys_method(matrix, epsilon, counter));  //create eigenvalues

      if (abs((3+l*4)-jac_meth(l+1)) < tolerance) //the analytical eigenvalues
      {
        cout << "Size of matrix needed for eigenvalue  " << l << "  is  " << n << endl; //print the result
        cout << "Found value: " << jac_meth(l+1) << endl;
        cout << endl;
        achived = true;
      }
      counter_loop = counter_loop + 1;

    }
  }
  return 0;
}
