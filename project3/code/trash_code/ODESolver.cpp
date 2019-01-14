#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <string>
#include <iomanip>
#include <fstream>


using namespace std;
using namespace arma;

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
  double dt;

public:
  ODESolver(vec mass, mat pos, mat vel, int n, double deltat)
  {
    m = mass;
    t = vec(n);
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

    for (int i = 0; i < t.n_elem; i++)
    {
      for (int j = 0; j < m.n_elem; j++)
      {
        a = gravity(j, x1);
        v2.row(j) = v1.row(j) + dt*a.t();
        x2.row(j) = x1.row(j) + dt*v1.row(j);
      }

      if (i%10000 == 0)
      {
        WriteData(filename, "x", x1);
        WriteData(filename, "y", x1);
        WriteData(filename, "z", x1);
      }
      // update new and old positions and velocities
      x1 = x2;
      v1 = v2;
    }
  }

  void VelocityVerlet(string filename)
  {
    vec ai;
    vec aii;

    for (int i = 0; i < t.n_elem; i++)
    {
      for (int j = 0; j < m.n_elem; j++)
      {
        ai = gravity(j, x1);
        x2.row(j) = x1.row(j) + dt*v1.row(j) + dt*dt/2*ai.t();

        aii = gravity(j, x2);
        v2.row(j) + dt/2*(aii.t() + ai.t());
      }
      if (i%100 == 0)
      {
        WriteData(filename, "x", x1);
        WriteData(filename, "y", x1);
        WriteData(filename, "z", x1);
      }
      x1 = x2;
      v1 = v2;
    }
  }
};

int main(int argc, char const *argv[])
{
  int number_of_planets = 2;

  // mass of planets
  vec mass = vec(number_of_planets);
  mass(0) = 1;
  mass(1) = 0.000003003;


  // initial position of planets
  mat pos = mat(number_of_planets, 3, fill::zeros);
  pos(1,0) = 1;

  // initial velocities
  mat vel = mat(number_of_planets, 3, fill::zeros);
  vel(1,1) = 2*pi;


  int n = pow(10, 6);
  double dt = 1.0/(n-1);

  //initialize class
  ODESolver sun_earth(mass, pos, vel, n, dt);

  //run Forward Euler and save data
  sun_earth.ForwardEuler("new");
  //sun_earth.VelocityVerlet("new");

  return 0;
}
