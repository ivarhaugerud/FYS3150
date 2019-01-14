#ifndef RELATIVISTICCORRECTION
#define RELATIVISTICCORRECTION

#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <string>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace arma;

 //inherits from this
#include "SolarSystemSimulator.hpp"

class RelativisticCorrection : public SolarSystemSimulator
{
public:
  bool perihelion;

  //Distances to sun to check for perihelion:
  double ppp = 0;
  double pp  = 0;
  double p   = 0;

  //Initializers
  RelativisticCorrection();
  RelativisticCorrection(vec mass, mat pos, mat vel, long long n, double T_max, string make_filename, bool E);

  //Member functions
  mat gravity(mat x);       //gravitational force in GR
  void VelocityVerlet();    //runs velocity verlet with perihelion test
  double PotentialEnergy(); //calculates potential energy
  void Perihelion();        //this tests if we are at the perihelion point, and if true writes position to file


};

#endif
