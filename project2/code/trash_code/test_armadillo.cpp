#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <armadillo>

//Including ibaries

//Our own methods:
using namespace std;
using namespace arma;

/*
*/

vec off(int n, mat matrix, double epsilon)
  {
      double max_value = 0;
      vec index = vec(2);
      for (int i = 0; i < n; i++)
      {
          for (int j = 0; j < n; j++)
          {
            if (abs(matrix(i, j)) > abs(max_value))
            {
              max_value = matrix(i, j);
              index(0) = i;
              index(1) = j;
            }
          }
      }
    if (max_value*max_value > epsilon)
      {
        return index;
      }
    else
      {
        return vec(2, fill::zeros) ;
      }
  }
