#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

/*
*/
double arctan2(double x, double h){
  return atan(x + h);
}
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

double richardsson(double x, double h){
  return (-arctan2(x, h) + 8*arctan2(x, h/2.0)-8*arctan2(x, -h/2.0) + arctan2(x, -h))/(6.0*h);
}

/*
*/

int main(int argc, char *argv[]){
  double x = sqrt(2);
  double h = 0;

  double richardsson_method = 0;
  float single_method_one = 0;
  float single_method_two = 0;
  double double_method_one = 0;
  double double_method_two = 0;

  ofstream out_data("data.txt");

  for( int i = 0; i < argc-1; i ++ ) {
      /* COMPUTING DERIVATIVES */
      h = pow(10, atof(argv[i+1]));
      single_method_one = derivate_1_single(x, h);
      single_method_two = derivate_2_single(x, h);
      double_method_one = derivate_1_double(x, h);
      double_method_two = derivate_2_double(x, h);
      richardsson_method = richardsson(x, h);
      cout << richardsson_method-1.0/3 << endl;
      out_data << h << " " << single_method_one -1.0/3 << " " <<  single_method_two -1.0/3 << " " << double_method_one -1.0/3 << " " << double_method_two -1.0/3 << " " << richardsson_method -1.0/3  << endl;
   }
  return 0;
}
