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
  long int n = 5*pow(10, 9);
  int T_max = 50000;

  double minimal_velocity = atof(argv[1]);
  double maximal_velocity = atof(argv[2]);
  double step             = atof(argv[3]);

  mat pos = mat(number_of_planets, dimensions);
  mat vel = mat(number_of_planets, dimensions);
  vec m   = vec(number_of_planets);

  // initial position of planets
  pos(1,0) = 1;
  pos(1, 1) = 0;
  pos(1, 2) = 0;

  pos(0, 0) = 0;
  pos(0, 1) = 0;
  pos(0, 2) = 0;

  // initial velocities
  vel(1, 0) = 0;
  vel(1, 0) = 0;

  vel(0, 0) = 0;
  vel(0, 1) = 0;
  vel(0, 2) = 0;

  m(0) = 1;
  m(1) = 0.000003003;

  for (double initial_velocity = minimal_velocity; initial_velocity <= maximal_velocity; initial_velocity += step)
  {
    //values for class
    vel(1,1) = initial_velocity;

    stringstream stream;
    stream << fixed << setprecision(4) << vel(1, 1);
    string velocity_string = "/escape_velocity/"+stream.str();      //precision of file name should only be two digits

    //initialize class
    SunAtRest solve(m, pos, vel, n, T_max, velocity_string, false);

    //run Forward Euler and save position data and time data
    solve.VelocityVerlet();
  }

  return 0;
}
