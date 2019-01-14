#include <iostream>
#include <cmath>
#include <armadillo>

#include "functions.cpp"

using namespace arma;
using namespace std;


int main()
{
  int n = 3;
  double a = 1;
  double d = 2;
  mat M = create_tridiagonal(n, a, d);

  cout << M << endl;

  rotation(M, 0, 1, 0, 1);

  cout << M << endl;

  return 0;
}
