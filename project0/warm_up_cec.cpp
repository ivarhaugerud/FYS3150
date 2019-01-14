#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

struct Precision
{
  float flt;
  double dbl;
};

class Diff
{
private:
  double (*f)(double x);

public:
  Diff(double function(double x))
  {
    f = function;
  }

  double ddiff2c(double h, double x = sqrt(2))
  {
    double result = (f(x + h) - f(x))/h;
    return result;
  }

  double ddiff3c(double h, double x = sqrt(2))
  {
    double result = (f(x + h) - f(x - h))/(2.0*h);
    return result;
  }
};

double func(double x)
{
  return atan(x);
}

void results2file(string filename, int argc, double h[], Precision vec2c[],\
                  Precision vec3c[])
{
  ofstream outfile(filename);
  if (!outfile.is_open())
    cout<<"Could not open file" << endl;
  else
  {
    for (int i = 1; i < argc; i++)
      outfile << setw(15) << h[i] << setw(15) << abs(vec2c[i].flt - 1.0/3) \
      << setw(15) << abs(vec2c[i].dbl - 1.0/3) << setw(15) << abs(vec3c[i].flt - 1.0/3) \
      << setw(15) << abs(vec3c[i].dbl - 1.0/3) << endl;
  }
}

int main(int argc, char *argv[])
{
  double h[argc-1];
  Precision vec2c[argc-1];
  Precision vec3c[argc-1];

  for (int i = 1; i < argc; i++)
  {
    float hexp = atof(argv[i]);
    h[i-1] = pow(10, hexp);
  }

  Diff arcustangens(func);

  for (int i = 1; i < argc; i++)
  {
    vec2c[i].dbl = arcustangens.ddiff2c(h[i]);
    vec2c[i].flt = arcustangens.ddiff2c(h[i]);
    vec3c[i].dbl = arcustangens.ddiff3c(h[i]);
    vec3c[i].flt = arcustangens.ddiff3c(h[i]);
  }

  results2file("data_cec.txt", argc, h, vec2c, vec3c);

  return 0;
}
