#include <iostream>
#include <cmath>
#include <armadillo>
//#include <cstdio>
#include "functions.hpp"

using namespace std;
using namespace arma;

int main()
{
  int n = 50;
  double epsilon = pow(10, -8);
  double d = 2;
  double a = -1;
  int counter = 0;

  mat A = create_tridiagonal(n, a, d);
  vec A_eig = sort(eig_sym(A));
  vec analytic_eigenvalues = sort(eigenvalues_finder(n, a, d));
  vec jac_meth = sort(jacobys_method(A, epsilon, counter));

  cout << "ARMADILLO EIG_SYM METHOD VS ANALYTICAL:" << endl;
  cout << endl;
  for (int t = 0; t < n; t++)
  {
    if (abs(analytic_eigenvalues(t)-A_eig(t)) < epsilon)
    {
      cout << "eigenvalue " << t << " correct:   " << A_eig(t) << endl;
    }
    else
    {
      cout << "ERROR:" << endl;
      cout << "Wrong calculation for eigenvalue " << t << endl;
      cout << "      Our calculation : " << analytic_eigenvalues(t) << endl;
      cout << "armadillo calculation : " << A_eig(t) << endl;
    }
  }

  cout << endl;
  cout << "ARMADILLO EIG_SYM METHOD VS JACOBYS METHOD:" << endl;
  cout << endl;

  for (int t = 0; t < n; t++)
  {
    if (abs(jac_meth(t)-A_eig(t)) < epsilon)
    {
      cout << "eigenvalue " << t << " correct:   " << A_eig(t) << endl;
    }
    else
    {
      cout << "ERROR:" << endl;
      cout << "Wrong calculation for eigenvalue " << t << endl;
      cout << "      Our calculation : " << jac_meth(t) << endl;
      cout << "armadillo calculation : " << A_eig(t) << endl;
    }
  }
  return 0;
}
