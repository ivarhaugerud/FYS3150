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

#include "json.hpp"
using json = nlohmann::json;
//include json

int main(int argc, char const *argv[])
{
  //structure of arguments from commande line
  // max time [years], number_of_timesteps, sun_at_rest, relativity, FE/VV, energy, plantes...

  int non_planet_arguments      = 6+1; //+1 due to the filename being argument zero
  int number_of_planets         = argc - non_planet_arguments;
  int dimensions                = 3;
  double T_max                  = atof(argv[1]);
  long long number_of_timesteps = pow(10, atof(argv[2]));

  //defining variabels from command line
  string sunrest      = argv[3];
  string relativistic = argv[4];
  string algorithm    = argv[5];
  string energy_      = argv[6];
  bool energy;

  //do we want to write energy data to file?
  if (energy_ == "true")
  {
    energy = true;
  }
  else
  {
    energy = false;
  }

  //do we want to run three unit tests?
  bool run_tests = true;

  //initialize position, veloicty and mass
  mat pos = mat(number_of_planets, dimensions);
  mat vel = mat(number_of_planets, dimensions);
  vec m   = vec(number_of_planets);

  //We use data from 2018-Sep-28 00:00:00.0000 TDB
  std::ifstream i("planet_data.json");
  json planet_data;
  i >> planet_data;

  //structure of values in planet_data.json
  //mass - initial position (x, y, z) - initial velocity(x, y, z)

  double sun_mass = planet_data["sun"][0]; //mass is index 0 and use solar mass as units

  //collect data from json file
  for (int i = 0; i < number_of_planets; i++)
  {
    m(i) = (double)planet_data[argv[i+non_planet_arguments]][0]/sun_mass; //mass is index 0
    for (int j = 0; j < dimensions; j++)
    {
      pos(i, j) = (double)planet_data[argv[i+non_planet_arguments]][j+1];          //already in AU
      vel(i, j) = (double)planet_data[argv[i+non_planet_arguments]][j+4]*365.2422; //converts to AU/year
    }
  }


  //run for sun at rest?
  if (sunrest == "true")
  {
    string planet3 = "planet";
    if (number_of_planets ==3)
      planet3 = argv[non_planet_arguments + 2];
    if (number_of_planets == 3 && planet3 == "jupiter")
    {
      //Increase mass of Jupiter.
      SunAtRest solve(m, pos, vel, number_of_timesteps, T_max, "jupiter_1", energy);
      solve.VelocityVerlet();

      vec m10 = m;
      m10(2) = m(2)*10; //Multiply Jupiter's mass by 10.
      solve.ResetInitialConditions(m10, pos, vel);
      solve.ChangeFileName("jupiter_10");
      solve.VelocityVerlet();

      vec m100 = m;
      m100(2) = m(2)*100; //Multiply Jupiter's mass by 100.
      solve.ResetInitialConditions(m100, pos, vel);
      solve.ChangeFileName("jupiter_100");
      solve.VelocityVerlet();


      vec m1000 = m;
      m1000(2) = m(2)*1000; //Multiply Jupiter's mass by 1000.
      solve.ResetInitialConditions(m1000, pos, vel);
      solve.ChangeFileName("jupiter_1000");
      solve.VelocityVerlet();
    }
    else
    {
      SunAtRest solve(m, pos, vel, number_of_timesteps, T_max, "sun_at_rest", energy);
      if (algorithm == "FE")
      {
        solve.ForwardEuler();
      }
      else if (algorithm == "VV")
      {
        solve.VelocityVerlet();
      }
      else
      {
        cout << "Method " << algorithm << " is not implemented." << endl;
      }
    }
  }

  //want to run relativistic?
  else if (relativistic == "true")
  {
    cout << "RUNNING RELATIVISTICCORRECTION" << endl;
    RelativisticCorrection solve(m, pos, vel, number_of_timesteps, T_max, "relativistic", energy);
    if (algorithm == "FE")
    {
      solve.ForwardEuler();
    }
    else if (algorithm == "VV")
    {
      solve.VelocityVerlet();
    }
    else
    {
      cout << "Method " << algorithm << " is not implemented." << endl;
    }
  }

  //want to run with sun movement and non relativistic
  else
  {
    if (number_of_planets ==3)
    {
      vel(0, 0) = 0;
      vel(0, 1) = 0;
      for (int i = 1; i < number_of_planets; i++)
      {
        vel(0, 0) -= vel(i, 0)*m(i);
        vel(0, 1) -= vel(i, 1)*m(i); //dont have to divide by one since it's in solar masses
      }
    }

    SolarSystemSimulator solve(m, pos, vel, number_of_timesteps, T_max, "simulator", energy);
    if (algorithm == "FE")
    {
      solve.ForwardEuler();
    }
    else if (algorithm == "VV")
    {
      solve.VelocityVerlet();
    }
    //incase of wrong input
    else
    {
      cout << "Method " << algorithm << " is not implemented." << endl;
    }

    //if we want to run unit test, if so run all 3
    if (run_tests)
    {
      SolarSystemSimulator solve(m, pos, vel, number_of_timesteps, T_max, "simulator", energy);

      double energy_start = solve.PotentialEnergy()+solve.KineticEnergy();
      double ang_mom_start = solve.AngularMomentum();

      solve.VelocityVerlet();
      double energy_stop = solve.PotentialEnergy()+solve.KineticEnergy();
      double ang_mom_stop = solve.AngularMomentum();

      cout << endl;

      //test conservation of energy
      if (abs(energy_start-energy_stop) < pow(10, -10))
      {
        cout << "PASSED:" << endl << "energy is conserved to an accuracy less than 10^(-10)" << endl;
      }
      else
      {
        cout << "FAILED:" << endl << "energy is NOT conserved to an accuracy less than 10^(-10)" << endl;
      }
      cout << endl;

      //test conservation angular momentum
      if (abs(ang_mom_start-ang_mom_stop) < pow(10, -12))
      {
        cout << "PASSED:" << endl << "Angular momentum is conserved to an accuracy less than 10^(-12)" << endl;
      }
      else
      {
        cout << "FAILED:" << endl << "Angular momentum is NOT conserved to an accuracy less than 10^(-12)" << endl;
      }
      cout << endl;

      //test calculation of length
      vec test_vector = vec(3);
      test_vector[0] = 1;
      test_vector[1] = 5;
      test_vector[2] = 7;
      double analytical_length = 8.6602540378443873;
      double numerical_length = solve.length(test_vector);

      if (abs(analytical_length-numerical_length) < pow(10, -8))
      {
        cout << "PASSED:" << endl << "Length of vector is calculated to an accuracy less than 10^(-8)" << endl;
      }
      else
      {
        cout << "FAILED:" << endl << "Length of vector is NOT calculated to an accuracy less than 10^(-8)" << endl;
      }
    }
  }
  return 0;
}
