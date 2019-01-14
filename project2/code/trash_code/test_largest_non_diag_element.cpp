#include <iostream>
#include <cmath>
#include <armadillo>
#include "functions.hpp"

using namespace std;
using namespace arma;

int main()
{
  int n = 5;
  double epsilon = 0;
  double largest_value = -72.9;

  mat A = mat(n, n, fill::eye);

  A(2, 3) = 5;
  A(3, 2) = 5;
  A(1, 4) = 72.8;
  A(4, 1) = 72.8;
  A(0, 3) = largest_value;
  A(3, 0) = largest_value;


  int k = 0;
  int l = 0;

  bool achieved = false;

  off(A, k, l, epsilon, achieved);

  if  (A(k, l) == largest_value)
  {
    cout << "Correct value found." << endl;
  }
  else
  {
    cout << "Error:" << endl;
    cout << "Largest absolute value given in program: " << largest_value << endl;
    cout << "                    Largest value found: " << A(k, l) << endl;
  }
  return 0;
}
