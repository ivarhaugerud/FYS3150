#include "SunAtRest.hpp"
#include <armadillo>


SunAtRest::SunAtRest():SolarSystemSimulator(){}


SunAtRest::SunAtRest\
                     (vec mass, mat pos, mat vel, long long n, double T_max, string make_filename, bool E)\
: SolarSystemSimulator(mass, pos, vel, n, T_max, make_filename, E)
{
  //paralell shifts the position of all plant in accordance to the change in position of the sun
  for (int i = 1; i < m.n_elem; i++)
  {
    x1.row(i) -= x1.row(0);
  }
  //sets the velocity and position of the sun to the origin
  x1(0,0) = 0;
  x1(0,1) = 0;
  x1(0,2) = 0;
  v1(0,0) = 0;
  v1(0,1) = 0;
  v1(0,2) = 0;
}

mat SunAtRest::gravity(mat x)
{
  mat acc = mat(m.n_elem, 3, fill::zeros);
  double G = 4*pi*pi;
  double force;
  for (int k = 0; k < m.n_elem; k++)
  {
    for (int j = 0; j < k; j++)
    {
      force = m(k)*m(j)*pow((x(j, 0)-x(k, 0))*(x(j, 0)-x(k, 0))\
                          + (x(j, 1)-x(k, 1))*(x(j, 1)-x(k, 1))\
                          + (x(j, 2)-x(k, 2))*(x(j, 2)-x(k, 2)), -3.0/2);

      //does not include  the sun in calculating the acceleration
      if (k != 0)
      {
        acc(k, 0) += force/m(k)*(x(j, 0)-x(k, 0));
        acc(k, 1) += force/m(k)*(x(j, 1)-x(k, 1));
        acc(k, 2) += force/m(k)*(x(j, 2)-x(k, 2));
      }

      if (j != 0)
      {
        acc(j, 0) -= force/m(j)*(x(j, 0)-x(k, 0));
        acc(j, 1) -= force/m(j)*(x(j, 1)-x(k, 1));
        acc(j, 2) -= force/m(j)*(x(j, 2)-x(k, 2));
      }
    }
  }
  return acc*G;
}

void SunAtRest::ResetInitialConditions(vec masses, mat pos, mat vel)
{
  //resets initial condition, but still keeps the sun at the origin
  //and performs the pararllel shift in positions
  m = masses;
  x1 = pos;
  v1 = vel;
  for (int i = 1; i < m.n_elem; i++)
  {
    x1.row(i) -= x1.row(0);
  }
  x1(0,0) = 0;
  x1(0,1) = 0;
  x1(0,2) = 0;
  v1(0,0) = 0;
  v1(0,1) = 0;
  v1(0,2) = 0;
}
