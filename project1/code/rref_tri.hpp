/*
Algorithm for general tridiagonal matrix
with b on diagonal and a and c off-diagonal
*/

using namespace std;

void rref_tri(double *a, double *b, double *c, double *u, double *f, int n,\
              double *btilde, double *ftilde)
  {

  btilde[1] = b[1];
  ftilde[1] = f[1];

  double quotient;

  for (int i = 2; i < n+1; i++)
    {
    quotient = a[i-1]/btilde[i-1];
    btilde[i] = b[i]-c[i-1]*quotient;
    ftilde[i] = f[i]-ftilde[i-1]*quotient;
    }

  u[n] = ftilde[n]/btilde[n];
  for (int i = n-1; i > 0; i--)
    {
    u[i] = (ftilde[i]-c[i]*u[i+1])/btilde[i];
    }
  }
