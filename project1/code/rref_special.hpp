/*
Algorithm for special tridiagonal matrix
with identical elements d on diagonal and identical elements
e off-diagonal.
*/


using namespace std;


void rref_special(double *u, double *f, int n, double *btilde, double *ftilde)
{
  ftilde[1] = f[1];

  for (int i = 2; i < n+1; i++)
    {
    ftilde[i] = f[i]+ftilde[i-1]/btilde[i-1];
    }

  u[n] = ftilde[n]/btilde[n];
  for (int i = n-1; i > 0; i--)
    {
    u[i] = (ftilde[i]+u[i+1])/btilde[i];
    }
}
