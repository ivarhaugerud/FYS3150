#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <string>
#include <iomanip>
#include <fstream>

#include "SolarSystemSimulator.hpp"
#include "SunAtRest.hpp"
#include "json.hpp"

using json = nlohmann::json;

int main(int argc, char const *argv[])
{
  // structure of arguments from commande line
  // max time [years], min timestep for loop, max timestep for loop

  int dimensions = 3;
  int number_of_planets = 3;
  double T_max  = 12;

  double minimal_N = atof(argv[1]);
  double maximal_N = atof(argv[2]);
  double step      = atof(argv[3]);

  mat pos = mat(number_of_planets, dimensions);
  mat vel = mat(number_of_planets, dimensions);
  vec m   = vec(number_of_planets);

  string planet_names[3] = {"sun", "earth", "jupiter"};

  //We use data from 2018-Sep-28 00:00:00.0000 TDB
  std::ifstream i("planet_data.json");
  json planet_data;
  i >> planet_data;

  double sun_mass = planet_data["sun"][0]; //mass is index 0

  for (int i = 0; i < number_of_planets; i++)
  {
    m(i) = (double)planet_data[planet_names[i]][0]/sun_mass; //mass is index 0

    for (int j = 0; j < dimensions; j++)
    {
      pos(i, j) = (double)planet_data[planet_names[i]][j+1];          //already in AU
      vel(i, j) = (double)planet_data[planet_names[i]][j+4]*365.2422; //converts to AU/year
    }
  }
  /*
  for (double n = minimal_N; n <= maximal_N; n += step)
    {
      stringstream stream;
      stream << fixed << setprecision(4) << n;
      string n_string = stream.str();

      SolarSystemSimulator solve(m, pos, vel, pow(10, n), T_max, "/3_body/"+n_string, true);
      solve.VelocityVerlet();
      solve.ForwardEuler();
    }
    */

    //Increase mass of Jupiter.
    vec new_mass = m;
    for (double n = minimal_N; n <= maximal_N; n += step)
      {
        stringstream stream;
        stream << fixed << setprecision(4) << n;
        string n_string = stream.str();

        for (int mass = 1; mass <= 1000; mass *= 10)
        {

          stringstream stream;
          stream << fixed << setprecision(4) << mass;
          string mass_string = stream.str();

          new_mass(2) = m(2)*mass;
          SunAtRest solve(new_mass, pos, vel, pow(10, n), T_max, "/3_body/"+n_string+"_"+mass_string, true);
          solve.VelocityVerlet();
        }
      }

  return 0;
}
