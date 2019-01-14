#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>


using namespace std;
using namespace arma;

int main(int argc, char *argv[])
{
  cube a = cube(3, 10, 2, fill::zeros);
  //cout << a << endl;
  //mat B = a.row(0);
  a(0, 0, 0) = 3;
  a(0, 0, 1) = 4;
  double c = a(0, 0, 1);
  cout << c << endl;
  vec b = vec(2);
  b(0) = 5;
  b(1) = 10;

  b += a.subcube(0, 0, 0, 0, 0, 1);
  cout << b << endl;
  return 0;
}
