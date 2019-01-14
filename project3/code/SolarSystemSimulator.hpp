#ifndef SOLARSYSTEMSIMULATOR
#define SOLARSYSTEMSIMULATOR

#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <string>
#include <iomanip>
#include <fstream>
#include <time.h>

using namespace std;
using namespace arma;

const double pi = 3.1415926535897932384; //will be used globally

class SolarSystemSimulator
{
public:
  //variabels which all functions will use
  double G = 4*pi*pi; //always the same
  vec m;
  mat x1;
  mat x2;
  mat v1;
  mat v2;

  long long integrationpoints;
  double dt;
  double T_max;

  int lines_of_data;

  double start;
  double stop;
  double time_used;
  string filename;

  bool energy;


  //Initializers:
  SolarSystemSimulator();
  SolarSystemSimulator(vec mass, mat pos, mat vel, long long n, double T_max, string make_filename, bool E);


  //Member functions:
  virtual mat gravity(mat x); //virtual so that it can be redefined by subclasses. Calculates forces

  double KineticEnergy();             //returns total kinetic energy
  virtual double PotentialEnergy();   //return total potential energy
  double AngularMomentum();           //return total angular momentum
  void TotalEnergy(string filename);  //writes total energy and angular momentum to file

  void ChangeFileName(string new_filename);                           //change file name?
  virtual void ResetInitialConditions(vec masses, mat pos, mat vel);  //reset inital conditions?

  virtual void ForwardEuler();      //run forward Euler for intialized system
  virtual void VelocityVerlet();    //run velocity Verlet for intialized system

  double length(vec r);             //calculates length of vectors

  void WriteData(string filename, mat matrix, string type); //write position data to file
  void DeleteDatafile(string filename, string type);        //delte data files
  void Timer(string start_or_stop, string filename);        //run timer of algorithms
};

#endif
