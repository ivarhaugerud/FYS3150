#include "functions.hpp"
#include <cmath>
#include <armadillo>
#include <iomanip>

/*
This file includes all the functions we use in this project
what each functions does overfladisk is shown in functions.hpp
*/
//using namespace arma;
using namespace std;

//function for writing results to file

//write results to file, given a matrix as input
void results2file(string filename, mat matrix)
  {
    int n = matrix.n_rows;

    ofstream outfile("data/" + filename + ".txt");
    if (!outfile.is_open())
      cout<<"Could not open file" << endl;
    else
    {
      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < matrix.n_cols; j++)
        {
        outfile << setw(15) << matrix(i, j);
        }
        outfile << endl;
      }
    }
  }


//create tridiagonal matrix
//loops over all tri diagonal elements
//given scaler inputs for diagonal and secondary diagonal elements
mat create_tridiagonal(int n, double a, double d)
{
  mat A = mat(n, n, fill::eye)*d;

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (i == j+1)
        A(i,j) = a;
      else if (i == j-1)
        A(i,j) = a;
    }
  }
  return A;
}

//create tridiagonal matrix
//loops over all tri diagonal elements
//given vector unputs for diagonal and secondary diagonal
mat create_tridiagonal_vector(int n, vec a, vec d)
{
  mat A = mat(n, n);

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (i==j)
        A(i, j) = d(i);
      if (i == j+1)
        A(i,j) = a(i);
      else if (i == j-1)
        A(i,j) = a(i);
    }
  }
  return A;
}

//uses the analytical expression from the task to check eigenvalues
vec eigenvalues_finder(int n, double a, double d)
{
  double pi = atan(1)*4.0;
  double cos_argument = pi/(n+1.0);
  vec lambda_j = vec(n);

  for (int j=0; j < n; j++)
  {
    lambda_j(j) = d + 2*a*cos(cos_argument*(j+1));
  }
  return lambda_j;
}

//loops over all the upper triangle part of the matrix to find largest value
void off(mat matrix, int &k, int &l, double epsilon, bool &achieved)
{
  int n = matrix.n_rows;
  double max_value = 0;
  for (int i = 0; i < n-1; i++)
  {
    for (int j = 1+i; j < n; j++)
      {
        if ((abs(matrix(i, j)) > abs(max_value)))
        {
          max_value = matrix(i, j);
          k = i;
          l = j;
        }
      }
    }
  if (abs(max_value) < epsilon)
  {
    achieved = true;
  }
}

//computes the rotation for the jacobi's method
void rotation(mat& A, double c, double s, int k, int l)
{
  double Aki;
  double Ali;
  double Aik;
  double Ail;

  for (int i = 0; i < A.n_rows; i++)
  {
    //Multiplying by S transposed
    Aki = A(k,i);
    Ali = A(l,i);
    A(k,i) = Aki*c - Ali*s;
    A(l,i) = Aki*s + Ali*c;
  }

  for (int i = 0; i < A.n_rows; i++)
  {
    //Multiplying by S
    Aik = A(i,k);
    Ail = A(i,l);
    A(i,k) = Aik*c - Ail*s;
    A(i,l) = Aik*s + Ail*c;
  }
}

//computes all of the needed things given the input, and returns the eignvalues in a vector
vec jacobys_method(mat matrix, double epsilon, int &counter)
{
  bool achieved = false;
  int k = 0;
  int l = 0;

  double t;
  double tau;
  double c;
  double s;

  while (!achieved)
  {
    //finds largest element, with indexes
    off(matrix, k, l, epsilon, achieved);
    //calculates the new angels and element sizes

    tau = (matrix(l, l)-matrix(k, k))/(2*matrix(k, l));

    if ( tau >= 0 )
       t =  1.0/(tau + sqrt(1.0 + tau*tau));
    else
       t = -1.0/(-tau +sqrt(1.0 + tau*tau));

    c = 1.0/sqrt(1+t*t);
    s = c*t;

    rotation(matrix, c, s, k, l); //does the rotation
    counter += 1;                 //counter for number of iteration

  }
  return matrix.diag();
}


//computes the eigenvector using jacobis method
mat eigenvector(mat A)
{
  int n = A.n_rows;
  mat R = mat(n, n, fill::eye);

  int k = 0;
  int l = 0;

  bool achieved = false;

  double t;
  double tau;
  double c;
  double s;

  double Rik;
  double Ril;

  double epsilon = pow(10, -8);
  //the algorithm is quite similar to jacobis method, but includes now the eigenvectors as well
  while (!achieved)
  {
    off(A, k, l, epsilon, achieved);

    tau = (A(l, l)- A(k, k))/(2*A(k, l));

    if ( tau >= 0 )
       t = 1.0/(tau + sqrt(1.0 + tau*tau));
    else
       t = -1.0/(-tau +sqrt(1.0 + tau*tau));

    c = 1/sqrt(1+t*t);
    s = c*t;

    rotation(A, c, s, k, l);

    for (int i = 0; i < R.n_rows; i++)
    {
      Rik = R(i,k);
      Ril = R(i,l);
      R(i,k) = Rik*c - Ril*s;
      R(i,l) = Rik*s + Ril*c;
    }
  }

  return R;
}
