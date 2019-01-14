//file containing the functions for the probability density and the energies
double psi_T1(mat r, double omega, double alpha, double beta)
{
  double r0 = r(0, 0)*r(0, 0) + r(0, 1)*r(0, 1) + r(0, 2)*r(0, 2);
  double r1 = r(1, 0)*r(1, 0) + r(1, 1)*r(1, 1) + r(1, 2)*r(1, 2);

  return exp(-alpha*omega*(r0 + r1));
  //removed 1/2 as we are going to square it either way
  //we do not need the constant C, as this will get canceld in the divition
}

double psi_T2(mat r, double omega, double alpha, double beta)
{
  double r0  = r(0, 0)*r(0, 0) + r(0, 1)*r(0, 1) + r(0, 2)*r(0, 2);
  double r1  = r(1, 0)*r(1, 0) + r(1, 1)*r(1, 1) + r(1, 2)*r(1, 2);
  double r12 = sqrt((r(0, 0)-r(1,0))*(r(0, 0)-r(1,0))
                  + (r(0, 1)-r(1,1))*(r(0, 1)-r(1,1))
                  + (r(0, 2)-r(1,2))*(r(0, 2)-r(1,2)));

  return exp(-alpha*omega*(r0+ r1))*exp(r12/(1+beta*r12));
  //removed 1/2 as we are going to square it either way
  //we do not need the constant C, as this will get canceld in the divition
}

double E_T0(mat r, double omega, double alpha, double beta)
{
  double r0 = r(0, 0)*r(0, 0) + r(0, 1)*r(0, 1) + r(0, 2)*r(0, 2);
  double r1 = r(1, 0)*r(1, 0) + r(1, 1)*r(1, 1) + r(1, 2)*r(1, 2);

  return 0.5*omega*omega*(r0 + r1)*(1 - alpha*alpha) + 3*alpha*omega;
}

double E_T1(mat r, double omega, double alpha, double beta)
{
  double r0 = r(0, 0)*r(0, 0) + r(0, 1)*r(0, 1) + r(0, 2)*r(0, 2);
  double r1 = r(1, 0)*r(1, 0) + r(1, 1)*r(1, 1) + r(1, 2)*r(1, 2);
  double r12 = sqrt((r(0, 0)-r(1,0))*(r(0, 0)-r(1,0))
                  + (r(0, 1)-r(1,1))*(r(0, 1)-r(1,1))
                  + (r(0, 2)-r(1,2))*(r(0, 2)-r(1,2)));

  return 0.5*omega*omega*(r0 + r1)*(1 - alpha*alpha) + 3*alpha*omega + 1.0/r12;
}

double E_T2(mat r, double omega, double alpha, double beta)
{
  double r0 = r(0, 0)*r(0, 0) + r(0, 1)*r(0, 1) + r(0, 2)*r(0, 2);
  double r1 = r(1, 0)*r(1, 0) + r(1, 1)*r(1, 1) + r(1, 2)*r(1, 2);
  double r12 = sqrt((r(0, 0)-r(1,0))*(r(0, 0)-r(1,0))
                  + (r(0, 1)-r(1,1))*(r(0, 1)-r(1,1))
                  +  (r(0, 2)-r(1,2))*(r(0, 2)-r(1,2)));
  double common_factor = (1+beta*r12);
  return 0.5*omega*omega*(r0 + r1)*(1 - alpha*alpha) + 3*alpha*omega + 1/r12
         + 1.0/(2*common_factor*common_factor)*(alpha*omega*r12
         - 1.0/(2*common_factor*common_factor) - 2.0/r12 + 2*beta/common_factor);
}
