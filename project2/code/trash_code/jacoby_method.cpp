#include <iostream>
#include <cmath>
#include <armadillo>

#include "functions.cpp"

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{
  int n = atoi(argv[1]);
  int a = -1;
  int d = 2;
  double epsilon = pow(10, -8);
  int achieved = 1;

  mat matrix = create_tridiagonal(n, a, d);
  //vec indexes = off(n, matrix, epsilon);

  while (achieved)
  {
    vec indexes = off(n, matrix, epsilon);

    int k = indexes(0);
    int l = indexes(1);

    if ( (k == 0) && (l == 0) )
      {
        achieved = 0;
      }

    else
    {
      double tau = (matrix(l, l)-matrix(k, k))/(2*matrix(k, l));
      double t = tau - sqrt(1+tau*tau);
      double c = 1/sqrt(1+t*t);
      double s = c*t;

      //cout << c << " " << s << endl;
      rotation(matrix, c, s, k, l);
    }
  }
  cout << sort(matrix.diag()) << endl;
  return 0;
}
