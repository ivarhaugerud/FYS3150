#ifndef SUNATREST
#define SUNATREST

#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <string>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace arma;

#include "SolarSystemSimulator.hpp"

class SunAtRest : public SolarSystemSimulator
{
public:
  //Initializers
  SunAtRest();
  SunAtRest(vec mass, mat pos, mat vel, long long n, double T_max, string make_filename, bool E);


  mat gravity(mat x); //calculates all forces except the ones acting on the sun
  void ResetInitialConditions(vec masses, mat pos, mat vel); //resets initial condition
};

#endif
