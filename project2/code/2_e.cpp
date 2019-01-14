#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <iomanip>
//libaries we need

#include "functions.hpp" //header file with all the functions

using namespace std;
using namespace arma;

// This program find both the eigenvectors and eigenvalues for the 3D harmonic oscillator potential with two electrons and stores the data for the 5 lowest energy eigenstates in a .txt file

int main(int argc, char *argv[])
{
  //dynamic memory allocation for the variable we want to alter
  int n          = atoi(argv[1]);
  double rho_max = atof(argv[2]);
  double epsilon = pow(10, -atof(argv[3]));
  double omega   = atof(argv[4]);

  double omega_omega = omega*omega; // reduce FLOPS
  double h = rho_max/(n+1);         // stepsize
  double hh = h*h;                  // reduce FLOPS
  double hh2 = 2.0/hh;              // reduce FLOPS
  double hh1 = -1.0/hh;             // reduce FLOPS
  int counter = 0;

  vec d   = vec(n+1);       // diagonal elements
  vec e   = vec(n+1);       // secondary-diagonal elements
  vec rho = vec(n+1);       // position variable

  for (int i = 0; i < n+1; i++)
  {
    rho(i) = (i+1)*h;                                       //definition of rho
    d(i)   = hh2 + omega_omega*rho(i)*rho(i) + 1./rho(i);   //definition of diagonal elements
    e(i)   = hh1;                                           //definition of secondary diagonals
  }

  mat matrix = create_tridiagonal_vector(n, e, d);          //create the tridiagonal matrix
  vec jac_meth = jacobys_method(matrix, epsilon, counter);  //calculate eigenvalues and store in vector
  mat eigvecs = eigenvector(matrix);                        //calculate eigenvectors

  //Want the five lowest states

  int number = 5;             //the number of lowest states we are interested in
  double energy = pow(10, 5); //randomly chosen maximum energy for searching
  double smaller_energy = 0;  //randomly chosen minimum energy for searching
  int m = 0;                  //dummy index
  vec indices = vec(number);  //vector to store indecies

  for (int i = 0; i < number; i++) //loop over number of states we are intererested in
  {
    for (int j = 1; j < jac_meth.n_elem; j++)    //loop over all eigenvalues
    {
      if (jac_meth(j) < energy && jac_meth(j) > smaller_energy) //check if the value is smaller or not
      {
        m = j;                //if so this is the new smallest eigenvalue, and we need to store it
        energy = jac_meth(j); //and save it as the current lowest energy
      }
    }
    smaller_energy = jac_meth(m);   //after the for-loop this is the smallest one
    energy = pow(10, 5);            //new dummy-energy
    indices(i) = m;                 //save the index
    cout << "Eigenvalue " << i << ":  " << smaller_energy << endl; //tell us what the lowest energies were
  }

  mat states = mat(eigvecs.n_rows, number+1); //use the indexes just stored to find the corresponding eigenvectors
  for (int i = 0; i < number; i++)            //loop over states and store them in matrix for writing
  {
    states.col(i+1) = eigvecs.col(indices(i));
  }

  for (int i = 0; i < n; i++) //want to store the position vector for plotting
  {
    states(i, 0) = rho(i);
  }

  stringstream stream;
  stream << fixed << setprecision(2) << omega;
  string s = stream.str();      //precision of file name should only be two digits

  results2file("wavefunc_two_e_"+s, states);  //save data
  return 0;
}
