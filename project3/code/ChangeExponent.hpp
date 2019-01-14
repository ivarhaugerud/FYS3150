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

#include "SolarSystemSimulator.hpp"

class ChangeExponent : public SolarSystemSimulator
{
public:
  //calculates the gravity with another exponent of Newtons gravitational law
  ChangeExponent(vec mass, mat pos, mat vel, long long n, double T_max, double exponent, string make_filename, bool E);
  mat gravity(mat x);
  double betta;
};

#endif
