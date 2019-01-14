#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

/*
*/
double arctan(double x){
  return atan(x);
}
double derivate_1_double(double x, double h){
  return (arctan(x+h)-arctan(x))/h;
}

double derivate_2_double(double x, double h){
  return (arctan(x+h)-arctan(x-h))/(2*h);
}


double derivate_1_single(float x, float h){
  return (arctan(x+h)-arctan(x))/h;
}

double derivate_2_single(float x, float h){
  return (arctan(x+h)-arctan(x-h))/(2*h);
}

/*
*/

int main(int argc, char *argv[]){
  double x = sqrt(2);
  double h = 0;

  float single_method_one = 0;
  float single_method_two = 0;
  double double_method_one = 0;
  double double_method_two = 0;

  cout << " " << endl;
  cout << "arctan(x): " << arctan(x) << endl;
  cout << " " << endl;
  ofstream out_data("data.txt");


  for( int i = 0; i < argc-1; i ++ ) {
      /* COMPUTING DERIVATIVES */
      h = pow(10, atof(argv[i+1]));
      single_method_one = derivate_1_single(x, h);
      single_method_two = derivate_2_single(x, h);
      double_method_one = derivate_1_double(x, h);
      double_method_two = derivate_2_double(x, h);

      out_data << h << " " << single_method_one -1.0/3 << " " <<  single_method_two-1.0/3 << " " << double_method_one -1.0/3 << " " << double_method_two -1.0/3 << endl;
   }
  return 0;
}
