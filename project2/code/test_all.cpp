#include <iostream>
#include <cmath>
#include <armadillo>
//#include <cstdio>
#include "functions.hpp"

// This function runs all our test, and pints wether or not they are Correct

using namespace std;
using namespace arma;


//the following are UBUNTU/LINUX ONLY terminal color codes.
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */


int main()
{
  // do not use dynamic memory allocation so we can run the test without thinking of input
  int n = 50;       //size of matrix
  double epsilon = pow(10, -8); //tolerance
  double d = 2;  //use scalers as diagonals
  double a = -1; //and secondary diagonal
  int counter = 0; //define counter

  mat B = create_tridiagonal(n, a, d);   //create matrix
  vec A_eig = sort(eig_sym(B));          //calc the aigenvalues using the armadillo package
  vec analytic_eigenvalues = sort(eigenvalues_finder(n, a, d));   //calc the analytic eigenvalues using our analytic method
  vec jac_meth = sort(jacobys_method(B, epsilon, counter));       //the eigenvalues using our jacobi method
  cout << endl;

  cout << "ARMADILLO EIG_SYM METHOD VS ANALYTICAL:" << endl;
  bool sucess = true;

  //test if they are the same
  for (int t = 0; t < n; t++)
  {
    if (!(abs(analytic_eigenvalues(t)-A_eig(t)) < epsilon))
    {
      sucess = false;
    }
  }
  if (!sucess)
    cout << RED << "ERROR: Eigenvalues not correct." << RESET << endl;
  else
    cout << GREEN << "Correct eigenvalues." << RESET << endl;

  //prints wether or not it were a success


  //do exactly the same for our method
  cout << endl;
  cout << "ARMADILLO EIG_SYM METHOD VS JACOBYS METHOD:" << endl;
  sucess = true;

  for (int t = 0; t < n; t++)
  {
    if (!(abs(jac_meth(t)-A_eig(t)) < epsilon))
    {
      sucess = false;
    }
  }
  if (!sucess)
    cout << RED << "ERROR: Eigenvalues not correct." << RESET <<endl;
  else
    cout << GREEN << "Correct eigenvalues." << RESET << endl;
  cout << endl;

  /*
  TEST 2 test our eigenvectors
  */

  mat A = mat(3, 3, fill::eye)*2;
  A(0,1) = 1;
  A(1,0) = 1;
  A(1,2) = 1;
  A(2,1) = 1;

  mat eig = eigenvector(A);
  double tol = pow(10, -8);

  //this is the analytical eigenvalues

  bool check00 = (abs(eig(0,0) - 0.5) < tol);
  bool check10 = (abs(eig(1,0) - (-1.0/sqrt(2))) < tol);
  bool check20 = (abs(eig(2,0) - 0.5) < tol);

  bool check01 = (abs(eig(0,1) - 0.5) < tol);
  bool check11 = (abs(eig(1,1) - (1.0/sqrt(2))) < tol);
  bool check21 = (abs(eig(2,1) - 0.5) < tol);

  bool check02 = (abs(eig(0,2) - (-1.0/sqrt(2))) < tol);
  bool check12 = (abs(eig(1,2) - 0) < tol);
  bool check22 = (abs(eig(2,2) - (1.0/sqrt(2))) < tol);

  cout << "TEST OUR EIGENVECTORS METHOD:" << endl;
  if  (check00 && check10 && check20\
    && check01 && check11 && check21\
    && check02 && check12 && check22)
    cout << GREEN << "Correct eigenvectors." << RESET << endl;
  else
    cout << RED << "Eigenvectors are wrong." << RESET << endl;

    //Check that eigenvectors correspond to correct eigenvalue:
    double val1 = 2.0 - sqrt(2.0);
    double val2 = 2.0 + sqrt(2.0);
    double val3 = 2.0;


    n = 0;
    vec jac = jacobys_method(A, tol, n);
    tol = pow(10, -2);

    cout << endl;
    cout << "TEST OUR JACOBIE'S METHOD" << endl; //test the jacobies method
    if (abs(jac(0) - val1) < tol && abs(jac(1) - val2) < tol\
     && abs(jac(2) - val3) < tol)
      cout << GREEN << "Eigenvalues correspond." << RESET << endl;
    else
      cout << RED << "Eigenvalues are wrong." << RESET << endl;

      /*
      test if our off function works
      which collects the largest absolute value
      */

      n = 5;
      epsilon = 0;
      double largest_value = -72.9;

      mat C = mat(n, n, fill::eye); //create matrix

      //fill matrix with elements
      C(2, 3) = 5;
      C(3, 2) = 5;
      C(1, 4) = 72.8;
      C(4, 1) = 72.8;
      C(0, 3) = largest_value;
      C(3, 0) = largest_value;

      //define index
      int k = 0;
      int l = 0;

      //test if the largest value is correct
      sucess = false;
      off(C, k, l, epsilon, sucess);
      cout << endl;
      cout << "TEST IF FIND LARGEST NON-DIAG ELEMENT:" << endl;
      if  (C(k, l) == largest_value)
      {
        cout << GREEN << "Correct value found." << RESET << endl;
      }
      else
      {
        cout << RED << "Error:" << endl;
        cout << "Largest absolute value given in program: " << largest_value << endl;
        cout << "                    Largest value found: " << A(k, l) << RESET << endl;
      }
    cout << endl;

    //has now tested all the test functions
  return 0;
}
