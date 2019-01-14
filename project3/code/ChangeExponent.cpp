#include "ChangeExponent.hpp"
#include <armadillo>

ChangeExponent::ChangeExponent\
                            (vec mass, mat pos, mat vel, long long n, double T_max, double exponent, string make_filename, bool E)\
: SolarSystemSimulator(mass, pos, vel, n, T_max, make_filename, E)
{
  betta = exponent;
}

mat ChangeExponent::gravity(mat x) //j = planet numner, i = time index
{
  mat acc = mat(m.n_elem, 3, fill::zeros);
  double G = 4*pi*pi;
  double force;

  for (int k = 0; k < m.n_elem; k++)
  {
    for (int j = 0; j < k; j++)
    {
      //everything is the same to the standard gravity, except for the exponent = betta
      force = m(k)*m(j)*pow((x(j, 0)-x(k, 0))*(x(j, 0)-x(k, 0))\
                          + (x(j, 1)-x(k, 1))*(x(j, 1)-x(k, 1))\
                          + (x(j, 2)-x(k, 2))*(x(j, 2)-x(k, 2)), -(1.0+betta)/2); //+1 because we multiply by \vec{r}

      acc(k, 0) += force/m(k)*(x(j, 0)-x(k, 0));
      acc(k, 1) += force/m(k)*(x(j, 1)-x(k, 1));
      acc(k, 2) += force/m(k)*(x(j, 2)-x(k, 2));

      acc(j, 0) -= force/m(j)*(x(j, 0)-x(k, 0));
      acc(j, 1) -= force/m(j)*(x(j, 1)-x(k, 1));
      acc(j, 2) -= force/m(j)*(x(j, 2)-x(k, 2));
    }
  }

  return acc*G;
}
