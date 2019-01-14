#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>


using namespace std;

/*
Imported files
*/

double** tridiagonal(int number, double* a, double* b, double* c)
{
  double** matrix;
  matrix = new double*[number];
  for(int i = 0; i < number; i++)
    matrix[i] = new double[number];

  for (int i = 0; i < number; i++)
    for (int j = 0; j < number; j++)
    {
      if (i == j)
      {
        matrix[i][j] = b[i];
      }
      else if (i == j+1)
      {
        matrix[i][j] = a[j];
      }
      else if (i == j-1)
      {
        matrix[i][j] = c[i];
      }
      else
      {
        matrix[i][j] = 0;
      }
    }
  return matrix;
}

double func(double var)
{
  return 100*exp(-10*var);
}

void rref_tri(double** matrix, double* b, int n)
{
  for (int i = 0; i < n-1; i++)
  {
    double coeff = matrix[i][i+1]/matrix[i][i];
    matrix[i][i+1] -= coeff*matrix[i][i];
    matrix[i+1][i+1] -= coeff*matrix[i+1][i];
    b[i+1] -= coeff*b[i];

    if (i > 0)
    {
      for (int j = 1; j <= i; j++)
      {
        coeff = matrix[i][i-j]/matrix[i][i];
        matrix[i][i-j] -= coeff*matrix[i][i];
        matrix[i+1][i-j] -= coeff*matrix[i+1][i];
        b[i-j] -= coeff*b[i];
      }
    }
    for (int j = 1; j < n; j++)
      {
        coeff = matrix[n-1][n-1-j]/matrix[n-1][n-1];
        matrix[n-1][n-1-j] -= coeff*matrix[n-1][n-1];
        b[n-1-j] -= coeff*b[n-1];
      }
    }

    for (int i = 0; i < n; i++)
    {
      matrix[i][i] /= matrix[i][i];
      b[i] /= matrix[i][i];
    }
  }


  void results2file(string filename, int argc, double x[], double f[])
  {
    ofstream outfile(filename + to_string(argc) + ".txt");
    if (!outfile.is_open())
      cout<<"Could not open file" << endl;
    else
    {
      for (int i = 0; i < argc; i++)
        outfile << setw(15) << x[i] << setw(15) << f[i] << endl;
    }
  }

int main(int argc, char *argv[])
{
  int n = atoi(argv[1]);

  double x[n];
  double u[n];
  double f[n];

  double h = 1.0/(n-1);
  for (int i = 0; i < n; i++)
    {
    x[i] = h*i;
    f[i] = func(x[i])*h*h;
    }

  double a[n-1];
  double b[n];
  double c[n-1];

  for (int i = 0; i < n; i++)
  {
    a[i] = 1;
    c[i] = 1;
    b[i] = (i+1)/i;
  }
  b[n] = -2;

  double** A = tridiagonal(n, a, b, c);

  rref_tri(A, f, n);

  /*
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
    cout << A[j][i] << " ";
    }
    cout << endl;
  }
  */
  results2file("data_task_1b_", n, x, f);

  for(int i = 0; i < n; i++)
  {
    delete [] A[i];
  }
  delete [] A;

  return 0;
}
