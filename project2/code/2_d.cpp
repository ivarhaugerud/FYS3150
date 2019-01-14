#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>

#include "functions.hpp" //headerfile with functions we want to use

using namespace std;
using namespace arma;

//This program finds the eigenvalues and eigenvectors for the single electron in 3D harmonic oscillator, and stores the data for the 5 lowest energy eigenstates in a .txt file

int main(int argc, char *argv[])
{
  //Dynamical memory allocation of variables we want to alter
  int n          = atoi(argv[1]);
  double rho_max = atof(argv[2]);
  double epsilon = pow(10, -atof(argv[3]));

  double h = rho_max/(n+1);         // stepsize
  double hh = h*h;                  // reduce FLOPS
  double hh2 = 2.0/hh;              // reduce FLOPS
  double hh1 = -1.0/hh;             // reduce FLOPS
  int counter = 0;                  // counts number of iternations

  mat j   = mat(3, 4, 5);
  vec d   = vec(n+1);       // diagonal elements
  vec e   = vec(n+1);       // secondary-diagonal elements
  vec rho = vec(n+1);       // position variable

  //Define the last elements in the array
  rho(n) = 0;
  d(n) = 0;
  e(n) = 0;

  for (int i = 1; i < n; i++)
  {
    rho(i) = i*h;                         //definition of rho
    d(i)   = hh2 + rho(i)*rho(i);         //diagonal elements definition
    e(i)   = hh1;                         //secondary-diagonal definition
  }

  mat matrix = create_tridiagonal_vector(n, e, d); //Function to create the tridiagonal matrix
  vec jac_meth = jacobys_method(matrix, epsilon, counter);  //calculates the eigenvvalues and return them as a vector
  mat eigvecs = eigenvector(matrix); //computes the eigenvectors of the given matrix input


  //Want the five lowest states
  int number = 5; //the number of states we are interested in
  double energy = pow(10, 5); //define something as the highest energy for first iteration for for-loop
  double smaller_energy = 0;  //define something as the lowest  energy for first iteration for for-loop
  int m = 0;                  //define the index which will be used
  vec indices = vec(number);  //the vector where we will store the indeces

  //loop over all eigenvalues and find the the lowest one, and it's corresponding index
  //we need the index to find the corresponding eigenvector
  for (int i = 0; i < number; i++)
  {
    for (int j = 1; j < jac_meth.n_elem; j++)
    {
      if (jac_meth(j) < energy && jac_meth(j) > smaller_energy)
      {
        m = j;
        energy = jac_meth(j);
      }
    }
    smaller_energy = jac_meth(m); //the smaller enery is now the smallest one found
    energy = pow(10, 5);          //and the largest one has to be redefined
    indices(i) = m;               //save the index
    cout << smaller_energy << endl;
  }

  mat states = mat(eigvecs.n_rows, number+1); //the states we are interested in +1 for for saving the x-axis elements
  for (int i = 0; i < number; i++) //loops over the number of states we are interested in
  {
    states.col(i+1) = eigvecs.col(indices(i));
  }


  for (int i = 0; i < n; i++)//adds the x-axis, we need this for plotting
  {
    states(i, 0) = rho(i);
  }

  results2file("wavefunc_one_e", states); //writes the eigenvectors to file
  return 0;
}
