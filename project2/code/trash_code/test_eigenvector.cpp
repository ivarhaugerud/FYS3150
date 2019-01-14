#include <iostream>
#include <cmath>
#include <armadillo>
#include "functions.hpp"

using namespace std;
using namespace arma;

int main()
{
  mat A = mat(3, 3, fill::eye)*2;
  A(0,1) = 1;
  A(1,0) = 1;
  A(1,2) = 1;
  A(2,1) = 1;

  mat eig = eigenvector(A);

  double tol = pow(10, -8);

  //Check eigenvectors:

  bool check00 = (abs(eig(0,0) - 0.5) < tol);
  bool check10 = (abs(eig(1,0) - (-1.0/sqrt(2))) < tol);
  bool check20 = (abs(eig(2,0) - 0.5) < tol);

  bool check01 = (abs(eig(0,1) - 0.5) < tol);
  bool check11 = (abs(eig(1,1) - (1.0/sqrt(2))) < tol);
  bool check21 = (abs(eig(2,1) - 0.5) < tol);

  bool check02 = (abs(eig(0,2) - (-1.0/sqrt(2))) < tol);
  bool check12 = (abs(eig(1,2) - 0) < tol);
  bool check22 = (abs(eig(2,2) - (1.0/sqrt(2))) < tol);

  if  (check00 && check10 && check20\
    && check01 && check11 && check21\
    && check02 && check12 && check22)
    cout << "Correct eigenvectors." << endl;
  else
    cout << "Eigenvectors are wrong." << endl;

  //Check that eigenvectors correspond to correct eigenvalue:
  double val1 = 2.0 - sqrt(2.0);
  double val2 = 2.0 + sqrt(2.0);
  double val3 = 2.0;


  int n = 0;
  vec jac = jacobys_method(A, tol, n);

  tol = pow(10, -2);

  if (abs(jac(0) - val1) < tol && abs(jac(1) - val2) < tol\
   && abs(jac(2) - val3) < tol)
    cout << "Eigenvalues correspond." << endl;
  else
    cout << "Eigenvalues are wrong." << endl;

  return 0;
}
