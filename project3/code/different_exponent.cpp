#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <string>
#include <iomanip>
#include <fstream>

#include "SolarSystemSimulator.hpp"
#include "ChangeExponent.hpp"

int main(int argc, char const *argv[])
{
  // structure of arguments from commande line
  // max time [years], min timestep for loop, max timestep for loop

  int dimensions = 3;
  int number_of_planets = 2;
  double T_max  = atof(argv[1]);

  double minimal_betta = atof(argv[2]);
  double maximal_betta = atof(argv[3]);
  double step          = atof(argv[4]);
  int n                = pow(10, atof(argv[5]));

  mat pos = mat(number_of_planets, dimensions);
  mat vel = mat(number_of_planets, dimensions);
  vec m   = vec(number_of_planets);

  //mass of planets
  m(0) = 1.0;
  m(1) = 0.000003003;

  // initial position of planets
  pos(1, 0) = 1;
  pos(1, 1) = 0;
  pos(1, 2) = 0;

  pos(0, 0) = 0;
  pos(0, 1) = 0;
  pos(0, 2) = 0;

  // initial velocities
  vel(1, 1) = 2*pi;
  vel(1, 0) = 0;
  vel(1, 0) = 0;

  vel(0, 0) = 0;
  vel(0, 1) = -vel(1, 1)*m(1)/m(0);
  vel(0, 2) = 0;

  for (double exponent_argument = minimal_betta; exponent_argument <= maximal_betta; exponent_argument += step)
    {
      stringstream stream;
      stream << fixed << setprecision(4) << exponent_argument;
      string expe = stream.str();      //precision of file name should only be two digits

      //initialize class
      cout << exponent_argument << endl;
      ChangeExponent solve(m, pos, vel, n, T_max, exponent_argument, "different_exponent_"+ expe, false);

      //run VelocityVerlet and save position data
      solve.VelocityVerlet();
    }
  return 0;
}
