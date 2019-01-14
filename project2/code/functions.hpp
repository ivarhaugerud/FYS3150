#ifndef FUNCTIONS
#define FUNCTIONS

//These are the functions we will use

#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

//to write results to file
void results2file(string filename, mat matrix);

//to calculate the largest non diag element
void off(mat matrix, int &k, int &l, double epsilon, bool &achieved);

//to rotate in jacobis method
void rotation(mat& A, double c, double s, int k, int l);

//find the analytical eigenvalues for the test
vec eigenvalues_finder(int n, double a, double d);

//create tridiagonal matrix with scalar input
mat create_tridiagonal(int n, double a, double d);

//create tridiagonal matrix with vector input
mat create_tridiagonal_vector(int n, vec a, vec d);

//to compute the whole of jacobys method
vec jacobys_method(mat matrix, double epsilon, int &counter);

//to find the eigenvectors of the input matrix using jaobys method
mat eigenvector(mat A);

#endif
