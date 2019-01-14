#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;


mat similarity(int n, int k, int l, double theta)
{
  mat A = mat(n,n,fill::eye);
  double c = cos(theta);
  double s = sin(theta);
  A(k,k) = c;
  A(l,l) = c;
  A(k,l) = s;
  A(l,k) = -s;

  return A;
}


int main(int argc, char** argv)
{
  //mat S = similarity(4, 1, 3, 3.14159265358979/2);

  //cout << S << endl;
  mat A = mat(2, 2, fill::zeros);
  A(0, 1) = 2;
  A(1, 0) = 2;
  cout << A << endl;
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, A);

  cout << eigval << endl;
  //cout << eig_sym(S) << endl;
  return 0;
}
