#include "RelativisticCorrection.hpp"
#include <armadillo>


RelativisticCorrection::RelativisticCorrection():SolarSystemSimulator(){}

RelativisticCorrection::RelativisticCorrection\
                            (vec mass, mat pos, mat vel, long long n, double T_max, string make_filename, bool E)\
: SolarSystemSimulator(mass, pos, vel, n, T_max, make_filename, E)
{
  //since relativistic will only be run for two planets, this warning is printed if not
  if (m.n_elem != 2)
  {
    cout << "WRONG NUMBER OF PLANETS!" << endl;
  }

  x1.row(1) -= x1.row(0);
  x1(0,0) = 0;
  x1(0,1) = 0;
  x1(0,2) = 0;
  v1(0,0) = 0;
  v1(0,1) = 0;
  v1(0,2) = 0;
}

mat RelativisticCorrection::gravity(mat x) //0 = planet numner, i = time index
{
  //calcualting the relativistic gravitational force
  mat acc = mat(m.n_elem, 3, fill::zeros);
  double cc = 63239.7263*63239.7263; //c^2 in AU/year
  double rr;
  double acceleration;
  vec l;

  l = vec(3);
  l(0) = (x1(0,1) - x1(1,1))*v1(1,2) - (x1(0,2) - x1(1,2))*v1(1,1);
  l(1) = (x1(0,2) - x1(1,2))*v1(1,0) - (x1(0,0) - x1(1,0))*v1(1,2);
  l(2) = (x1(0,0) - x1(1,0))*v1(1,1) - (x1(0,1) - x1(1,1))*v1(1,0);

  double ll = l(0)*l(0) + l(1)*l(1) + l(2)*l(2);

  rr = (x(0, 0)-x(1, 0))*(x(0, 0)-x(1, 0))\
     + (x(0, 1)-x(1, 1))*(x(0, 1)-x(1, 1))\
     + (x(0, 2)-x(1, 2))*(x(0, 2)-x(1, 2));

  acceleration = m(0)*pow(rr, -3.0/2)*(1 + 3*ll/(rr*cc));

  acc(1, 0) += acceleration*(x(0, 0)-x(1, 0));
  acc(1, 1) += acceleration*(x(0, 1)-x(1, 1));
  acc(1, 2) += acceleration*(x(0, 2)-x(1, 2));

  return acc*G;
}


void RelativisticCorrection::VelocityVerlet()
{
  mat ai;
  mat aii;
  filename += "_VV_";

  //make files ready
  DeleteDatafile(filename, "coordinates");
  DeleteDatafile("perihelion_relativity", "coordinates");
  //if (energy) DeleteDatafile(filename + "energy", "energy");

  for (long long i = 0; i < integrationpoints; i++)
  {
    ai = gravity(x1);
    x2 = x1 + dt*v1 + dt*dt/2*ai;

    aii = gravity(x2);
    v2 = v1 + dt/2*(aii + ai);

    if (i%lines_of_data == 0)
    //{
      WriteData(filename, x1, "coordinates");
      if (energy)
        TotalEnergy(filename + "energy");
    //}

    Perihelion();

    x1 = x2;
    v1 = v2;
  }
  //always write last timestep
  WriteData(filename, x1, "coordinates");
  if (energy)
    TotalEnergy(filename + "energy");
}

double RelativisticCorrection::PotentialEnergy()
{
  //calculates potential energy
  double V = 0;
  double distance;

  for (int i = 0; i < m.n_elem; i++)
  {
    for (int j = 0; j < i; j++)
    {
      distance = 0;
      for (int k = 0; k < 3; k++)
      {
        distance += (x1(i, k) - x1(j, k))*(x1(i, k) - x1(j, k));
      }
      distance = sqrt(distance);

      vec l = vec(3);
      l(0) = (x1(0,1) - x1(1,1))*v1(1,2) - (x1(0,2) - x1(1,2))*v1(1,1);
      l(1) = (x1(0,2) - x1(1,2))*v1(1,0) - (x1(0,0) - x1(1,0))*v1(1,2);
      l(2) = (x1(0,0) - x1(1,0))*v1(1,1) - (x1(0,1) - x1(1,1))*v1(1,0);

      double ll = l(0)*l(0) + l(1)*l(1) + l(2)*l(2);

      double rr = (x1(0, 0)-x1(1, 0))*(x1(0, 0)-x1(1, 0))\
                + (x1(0, 1)-x1(1, 1))*(x1(0, 1)-x1(1, 1))\
                + (x1(0, 2)-x1(1, 2))*(x1(0, 2)-x1(1, 2));

      double cc = 63239.7263*63239.7263;

      V += -G*m(i)*m(j)/distance*(1 + 3*ll/(rr*cc));
    }
  }
  return V;
}

void RelativisticCorrection::Perihelion()
{
  //this tests if we are at the perihelion point, and if true writes position to file
  pp = length(x1.row(1).t());
  p = length(x2.row(1).t());
  perihelion = pp < ppp && pp < p;
  if (perihelion)
  {
    WriteData("perihelion_relativity", x1, "coordinates");
  }
  ppp = pp;
}
