#include "MainClass.hpp"

MainClass::MainClass()
{}

//all variabels we neeed for the class, in addition to the two functions for the energy and wave function
MainClass::MainClass(string save_name, int amount_of_data, double OMEGA,
                    int nr_particles, double ALPHA, double BETA, double STEP,
                    double wavefunc(mat, double, double, double),
                    double energyfunc(mat, double, double, double))
{
  //the file name for saving data, and number of lines when saving data
  filename   = save_name;
  data_lines = amount_of_data;

  //save the variabels for the class
  step = STEP;
  alpha = ALPHA;
  beta = BETA;
  omega = OMEGA;

  //the two functions we will use
  psi = wavefunc;
  calcE = energyfunc;

  // fixes the global variables
  E = 0;
  nr_accepted = 0;

  //matrix and vectors we need
  R_trial = mat(2, 3, fill::zeros);           //the  trial  position, 2 electrons and 3 dimensions
  R = mat(2, 3, fill::zeros);                 //the current position, 2 electrons and 3 dimensions
  expectationvalues = vec(5, fill::zeros);    //all expectationvalues
  average_expectations = vec(5, fill::zeros); //the average of the expectation values

  //distribute the position equally for each particle
  R(0, 0) =  0.1;
  R(1, 0) = -0.1;

  //the random number genrator
  std::mt19937 generator (std::clock());
  std::uniform_real_distribution<double> zero_to_one_distribution(0.0, 1.0);
}

//the metropolis algorithm
void MainClass::Metropolis()
{
    //update trial wave function of particle 0
    R_trial(0, 0) = R(0, 0) + (zero_to_one_distribution(generator)-0.5)*step;
    R_trial(0, 1) = R(0, 1) + (zero_to_one_distribution(generator)-0.5)*step;
    R_trial(0, 2) = R(0, 2) + (zero_to_one_distribution(generator)-0.5)*step;

    //update trial wave function of particle 1
    R_trial(1, 0) = R(1, 0) + (zero_to_one_distribution(generator)-0.5)*step;
    R_trial(1, 1) = R(1, 1) + (zero_to_one_distribution(generator)-0.5)*step;
    R_trial(1, 2) = R(1, 2) + (zero_to_one_distribution(generator)-0.5)*step;

    //calculate the probability of accepting new position
    double w = psi(R_trial, omega, alpha, beta)/psi(R, omega, alpha, beta);

    //the metropolis test if we should accept the new position
    if (w >= 1 || zero_to_one_distribution(generator) <= w)
      {
        //if accepted we need to update the poisition, the energy, and the number of accepts
        R = R_trial;
        E = calcE(R, omega, alpha, beta);
        nr_accepted += 1;
      }
}

//function for equilibration, only input is number of cycles for equilibrating
void MainClass::Equilibrate(int equilibration_cycles)
  {
    for (int equil_cycle = 0; equil_cycle < equilibration_cycles; equil_cycle++)
      {
        Metropolis();
      }
    reset(); //need to reset the expectation values
    E = 0;   //and set the energy to 0
  }

//calculate the potential energy, used for testing the virial theorem
double MainClass::calc_potential(mat r)
{
  double r12 = sqrt((r(0, 0)-r(1,0))*(r(0, 0)-r(1,0))
                  + (r(0, 1)-r(1,1))*(r(0, 1)-r(1,1))
                  + (r(0, 2)-r(1,2))*(r(0, 2)-r(1,2)));

  return 0.5*omega*omega*( r(0, 0)*r(0, 0) + r(0, 1)*r(0, 1) + r(0, 2)*r(0, 2)
                         + r(1, 0)*r(1, 0) + r(1, 1)*r(1, 1) + r(1, 2)*r(1, 2)) + 1/r12;
  //returns the potential energy for the given position
}

void MainClass::Run(int nr_cycles)
{
  //fix the number of digits used in writing the file
  stringstream alpha_string;
  alpha_string << fixed << setprecision(3) << alpha;

  stringstream beta_string;
  beta_string << fixed << setprecision(3) << beta;


  //open all files
  //the expectationvalues
  ofstream outfile_expect("../data/" +  filename  +"_expectationvalues_" + alpha_string.str() + "_" + beta_string.str() + ".txt", std::ios_base::app);
  if (!outfile_expect.is_open())
     cout<<"Could not open file" << endl;

   //the number of accepts
   ofstream outfile_accept("../data/" +  filename + "_accepts_" + alpha_string.str() + "_" + beta_string.str() + ".txt", std::ios_base::app);
   if (!outfile_accept.is_open())
      cout<<"Could not open file" << endl;

  //the average of the expectation values
  ofstream outfile("../data/" +  filename + "_" + alpha_string.str() + ".txt", std::ios_base::app);
  if (!outfile.is_open())
     cout<<"Could not open file" << endl;

  //how many lines of data we want
  int denominator = nr_cycles/data_lines;

  //initial calculation of energy
  E  = calcE(R, omega, alpha, beta);
  int counter = 0;

  //use energy and position to calculate the expectationvalues: energy ,energy2, distance, potential energy, 1/distance
  expectationvalues(0) += E;
  expectationvalues(1) += E*E;
  expectationvalues(2) += sqrt((R(0, 0)-R(1,0))*(R(0, 0)-R(1,0)) + (R(0, 1)-R(1,1))*(R(0, 1)-R(1,1)) +  (R(0, 2)-R(1,2))*(R(0, 2)-R(1,2)));
  expectationvalues(3) += calc_potential(R);
  expectationvalues(4) += 1/(sqrt((R(0, 0)-R(1,0))*(R(0, 0)-R(1,0)) + (R(0, 1)-R(1,1))*(R(0, 1)-R(1,1)) +  (R(0, 2)-R(1,2))*(R(0, 2)-R(1,2))));

  //run for the number of monte carlo cycles
  for (int cycle = 1; cycle < nr_cycles; cycle++)
  {
    //run metropolis
    Metropolis();

    //update expectation values
    expectationvalues(0) += E;
    expectationvalues(1) += E*E;
    expectationvalues(2) += sqrt((R(0, 0)-R(1,0))*(R(0, 0)-R(1,0)) + (R(0, 1)-R(1,1))*(R(0, 1)-R(1,1)) +  (R(0, 2)-R(1,2))*(R(0, 2)-R(1,2)));
    expectationvalues(3) += calc_potential(R);
    expectationvalues(4) += 1/(sqrt((R(0, 0)-R(1,0))*(R(0, 0)-R(1,0)) + (R(0, 1)-R(1,1))*(R(0, 1)-R(1,1)) +  (R(0, 2)-R(1,2))*(R(0, 2)-R(1,2))));
    int denominator = pow(1.1, counter);

    //calculate the averages the number amount of times
    if ( cycle%denominator == 0)
    {
      counter += 1;
      //update the averages
      average_expectations(0) += expectationvalues(0)/(cycle+1);
      average_expectations(1) += expectationvalues(1)/(cycle+1);
      average_expectations(2) += expectationvalues(2)/(cycle+1);
      average_expectations(3) += expectationvalues(3)/(cycle+1);
      average_expectations(4) += expectationvalues(4)/(cycle+1);

      //write the expectation values
      write_expectation_values(cycle+1, outfile_expect);
    }
  }
   //write the expectationvalues for the final state
  write_expectation_values(nr_cycles, outfile_expect);

  //update the expectation values
  average_expectations(0) += expectationvalues(0)/(nr_cycles);
  average_expectations(1) += expectationvalues(1)/(nr_cycles);
  average_expectations(2) += expectationvalues(2)/(nr_cycles);
  average_expectations(3) += expectationvalues(3)/(nr_cycles);
  average_expectations(4) += expectationvalues(4)/(nr_cycles);

  //write accepts and the the average of the expectation values
  write_accepts(nr_cycles, outfile_accept);
  write(outfile);
}

//reset all expectation values
void MainClass::reset()
{
  expectationvalues(0) = 0;
  expectationvalues(1) = 0;
  expectationvalues(2) = 0;
  expectationvalues(3) = 0;
  expectationvalues(4) = 0;

  average_expectations(0) = 0;
  average_expectations(1) = 0;
  average_expectations(2) = 0;
  average_expectations(3) = 0;
  average_expectations(4) = 0;

  nr_accepted = 0;
}

//write expectationvalues to file
void MainClass::write_expectation_values(int cycle_number, ofstream& OutputFile)
{
  OutputFile << cycle_number << " " << expectationvalues(0)/cycle_number << " " << expectationvalues(1)/cycle_number << " " << expectationvalues(2)/cycle_number << "\n";
}

//write number of accepts to file
void MainClass::write_accepts(int cycle_nr, ofstream& OutputFile)
{
  OutputFile << step << " " << " " << float(nr_accepted)/cycle_nr << "\n";
}

//write the averages to file
void MainClass::write(ofstream& OutputFile)
{
  OutputFile << beta << " " << omega << " " << average_expectations(0)/data_lines << " " << average_expectations(1)/data_lines << " " << average_expectations(2)/data_lines << " " << average_expectations(3)/data_lines << " " << average_expectations(4)/data_lines << "\n";
}

//change the value of alpha
void MainClass::change_alpha(double new_alpha)
{
  alpha = new_alpha;
}

//change the value of beta
void MainClass::change_beta(double new_beta)
{
  beta = new_beta;
}

//change the value of the step size
void MainClass::change_step(double new_step)
{
  step = new_step;
}

//change the value of omega
void MainClass::change_omega(double new_omega)
{
  omega = new_omega;
}

//overwrite an existing data file and delete it
void MainClass::DeleteDatafile()
{
  //fix the number of digits used in writing the file
  stringstream alpha_string;
  alpha_string << fixed << setprecision(3) << alpha;

  stringstream beta_string;
  beta_string << fixed << setprecision(3) << beta;

  //the expectationvalues
  ofstream outfile_expect("../data/" +  filename  +"_expectationvalues_" + alpha_string.str() + "_" + beta_string.str() + ".txt");
  if (!outfile_expect.is_open())
     cout<<"Could not open file" << endl;

   //the number of accepts
   ofstream outfile_accept("../data/" +  filename + "_accepts_" + alpha_string.str() + "_" + beta_string.str() + ".txt");
   if (!outfile_accept.is_open())
      cout<<"Could not open file" << endl;

  //open file
  ofstream outfile("../data/" +  filename + "_" + alpha_string.str() + ".txt");
  if (!outfile.is_open())
     cout<<"Could not open file" << endl;
}
