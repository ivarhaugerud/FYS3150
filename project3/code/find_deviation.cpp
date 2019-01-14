#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <string>
#include <iomanip>
#include <fstream>

#include "SolarSystemSimulator.hpp"
#include "SunAtRest.hpp"
#include "RelativisticCorrection.hpp"


int main(int argc, char const *argv[])
{
  // structure of arguments from commande line
  // max time [years], min timestep for loop, max timestep for loop

  int dimensions = 3;
  int number_of_planets = 2;

  double T_max  = atof(argv[1]);
  int minimal_exponent = atoi(argv[2]);
  int maximal_exponent = atoi(argv[3]);
  int n;


  mat pos = mat(number_of_planets, dimensions);
  mat vel = mat(number_of_planets, dimensions);
  vec m   = vec(number_of_planets);

  // initial position of planets
  pos(1, 0) = 1;
  pos(1, 1) = 0;
  pos(1, 2) = 0;

  pos(0, 0) = 0;
  pos(0, 1) = 0;
  pos(0, 2) = 0;

  // initial velocities
  vel(1,1) = 2*pi;
  vel(1, 0) = 0;
  vel(1, 0) = 0;

  vel(0, 0) = 0;
  vel(0, 1) = 0;
  vel(0, 2) = 0;

  m(0) = 1;
  m(1) = 0.000003003;

  for (double exponent = minimal_exponent; exponent <= maximal_exponent; exponent += 0.1)
  {
    //values for class
    n = pow(10, exponent);

    stringstream stream;
    stream << fixed << setprecision(2) << exponent;
    string expe = stream.str();      //precision of file name should only be two digits

    //initialize class
    SunAtRest solve(m, pos, vel, n, T_max, expe, false);

    //run Forward Euler and save position data and time data
    solve.Timer("start", "FE");
    solve.ForwardEuler();
    solve.Timer("stop", "FE");

    //initialize class
    SunAtRest solve2(m, pos, vel, n, T_max, expe, false);

    //run Forward Euler and save position data and time data
    solve2.Timer("start", "VV");
    solve2.VelocityVerlet();
    solve2.Timer("stop", "VV");
  }

  return 0;
}
