#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <string>
#include <iomanip>
#include <fstream>

#include "json.hpp"


using namespace std;
using namespace arma;
using json = nlohmann::json;

const double pi = 3.1415926535897932384;


class ODESolver
{
protected:
  vec m;
  mat x1;
  mat x2;
  mat v1;
  mat v2;
  vec t;
  int lines_of_data;
  double dt;

public:
  ODESolver(vec mass, mat pos, mat vel, int n, double deltat)
  {
    m = mass;
    t = vec(n);
    lines_of_data = n/100+2;
    // create time array
    dt = deltat;
    for (int i = 0; i < n; i++)
    {
      t(i) = dt*i;
    }

    //for current time step
    x1 = mat(m.n_elem, 3); //Planets, time, space-coordinate
    v1 = mat(m.n_elem, 3); //Planets, time, space-coordinate

    //for next time step
    x2 = mat(m.n_elem, 3); //Planets, time, space-coordinate
    v2 = mat(m.n_elem, 3); //Planets, time, space-coordinate

    x1 = pos;     //Initial position
    v1 = vel;     //Initial velocity
  }

  vec gravity(int j, mat x) //j = planet numner, i = time index
  {
    vec acc = vec(3, fill::zeros);
    double G = 4*pi*pi;
    for (int k = 0; k < m.n_elem; k++)
    {
      if (k != j)
      {
        double acceleration = m(k)*pow((x(j, 0)-x(k, 0))*(x(j, 0)-x(k, 0))\
                                     + (x(j, 1)-x(k, 1))*(x(j, 1)-x(k, 1))\
                                     + (x(j, 2)-x(k, 2))*(x(j, 2)-x(k, 2)), -3.0/2);
        acc(0) += acceleration*(x(j, 0)-x(k, 0));
        acc(1) += acceleration*(x(j, 1)-x(k, 1));
        acc(2) += acceleration*(x(j, 2)-x(k, 2));
      }
    }

    return -acc*G;
  }


  void DeleteDatafile(string filename, string position)
  {
    ofstream outfile("data/"+ filename +"_" + position + ".txt");
    if (!outfile.is_open())
      cout<<"Could not open file" << endl;
  }

  void WriteData(string filename, string position, mat matrix)
    {
      std::ofstream outfile;
      outfile.open("data/"+ filename +"_" + position + ".txt", std::ios_base::app);

      if (!outfile.is_open())
        cout << "Could not open file" << endl;
      else

      if (position == "x")
      {
        for (int l = 0; l < matrix.n_rows; l++)
        {
          outfile << setw(15) << x1(l, 0);
        }
        outfile << endl;
      }

      if (position == "y")
      {
        for (int l = 0; l < matrix.n_rows; l++)
        {
          outfile << setw(15) << x1(l, 1);
        }
        outfile << endl;
      }

      if (position == "z")
      {
        for (int l = 0; l < matrix.n_rows; l++)
        {
          outfile << setw(15) << x1(l, 2);
        }
        outfile << endl;
      }
    }

  void ForwardEuler(string filename)
  {
    vec a;
    DeleteDatafile(filename, "x");
    DeleteDatafile(filename, "y");
    DeleteDatafile(filename, "z");

    for (int i = 0; i < t.n_elem; i++)
    {
      for (int j = 0; j < m.n_elem; j++)
      {
        a = gravity(j, x1);
        v2.row(j) = v1.row(j) + dt*a.t();
        x2.row(j) = x1.row(j) + dt*v1.row(j);
      }

      if (i%lines_of_data == 0)
      {
        WriteData(filename, "x", x1);
        WriteData(filename, "y", x1);
        WriteData(filename, "z", x1);
      }
      // update new and old positions and velocities
      x1 = x2;
      v1 = v2;
    }
    WriteData(filename, "x", x1);
    WriteData(filename, "y", x1);
    WriteData(filename, "z", x1);
  }

  void VelocityVerlet(string filename)
  {
    vec ai;
    vec aii;
    DeleteDatafile(filename, "x");
    DeleteDatafile(filename, "y");
    DeleteDatafile(filename, "z");

    for (int i = 0; i < t.n_elem; i++)
    {
      for (int j = 0; j < m.n_elem; j++)
      {
        ai = gravity(j, x1);
        x2.row(j) = x1.row(j) + dt*v1.row(j) + dt*dt/2*ai.t();

        aii = gravity(j, x2);
        v2.row(j) + dt/2*(aii.t() + ai.t());
      }
      if (i%lines_of_data == 0)
      {
        WriteData(filename, "x", x1);
        WriteData(filename, "y", x1);
        WriteData(filename, "z", x1);
      }
      x1 = x2;
      v1 = v2;
    }
    WriteData(filename, "x", x1);
    WriteData(filename, "y", x1);
    WriteData(filename, "z", x1);
  }
};

int main(int argc, char const *argv[])
{
  int number_of_planets = 2;
  int number_of_dimensions = 3;

  int minimal_exponent = atoi(argv[1]);
  int maximal_exponent = atoi(argv[2]);

  // mass of planets
  vec mass = vec(number_of_planets);
  mass(0) = 1;
  mass(1) = 0.000003003;

  // initial position of planets
  mat pos = mat(number_of_planets, number_of_dimensions, fill::zeros);
  pos(1,0) = 1;

  // initial velocities
  mat vel = mat(number_of_planets, number_of_dimensions, fill::zeros);
  vel(1,1) = 2*pi;

  int n;
  double dt;

  std::ifstream i("json_file.json");
  json j;
  i >> j;

  double a = j["sun"][0];
  //cout << j["sun"][0] << endl;

  for (double exponent = minimal_exponent; exponent <= maximal_exponent; exponent += 0.25)
  {
    //values for class
    n = pow(10, exponent);
    dt = 1.0/(n-1);

    stringstream stream;
    stream << fixed << setprecision(2) << exponent;
    string expe = stream.str();      //precision of file name should only be two digits

    //initialize class
    vel(1,1) = 2*pi;
    pos(1,0) = 1;
    ODESolver sun_earth1(mass, pos, vel, n, dt);

    //run Forward Euler and save data
    sun_earth1.ForwardEuler("FE_"+ expe);

    //initialize class
    vel(1,1) = 2*pi;
    vel(0,1) = -2*pi*mass(1)/mass(0);
    pos(1,0) = 1;
    ODESolver sun_earth_2(mass, pos, vel, n, dt);
    sun_earth_2.VelocityVerlet("VV_"+ expe);
  }
  return 0;
}
