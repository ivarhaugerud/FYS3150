#include "MainClass.hpp"

MainClass::MainClass()
{}

MainClass::MainClass(int size, string save_name, int amount_of_data)
{
  //initalize everything which is decleared in .hpp
  L = size;
  M = 0;
  E = 0;
  count_accepts = 0;

  spins = Mat<int>(L, L);
  index_vector = Col<int>(L+2);
  expectationvalues = vec(5, fill::zeros);
  boltzmann_factor  = vec(17, fill::zeros);
  filename = save_name;
  data_lines = amount_of_data;

  std::mt19937 generator (std::clock());
  std::uniform_real_distribution<double> zero_to_one_distribution(0.0, 1.0);
  std::uniform_int_distribution<int> zero_to_L_distribution(0, L-1);

  //craete the index vector
  index_vector = Col<int>(L+2);
  for (int i = 1; i < L+2; i++)
  {
    index_vector(i) = i-1;
  }
  index_vector(0) = L-1;
  index_vector(L+1) = 0;
}

//gives a random initialization
void MainClass::initialize_random()
  {
    for (int i = 0; i < L; i++)
    {
      for (int j = 0; j < L; j++)
      {
        spins(i, j) = (zero_to_one_distribution(generator) > 0.5) ? 1 : -1;
      }
    }
  }

//initialize everything in the same direction
void MainClass::initialize_up()
  {
    for (int i = 0; i < L; i++)
    {
      for (int j = 0; j < L; j++)
      {
        spins(i, j) = 1;
      }
    }
  }

//calculates total energy of system
double MainClass::calc_energy()
  {
    int energy = 0;
    //loops over every interaction using the index vector
    for (int i = 0; i < L; i++)
    {
      for (int j= 0; j<L; j++)
      {
        energy -= spins(i, j)*( spins( index_vector(i+1), index_vector(j) ) + spins(index_vector(i), index_vector(j+1)));
      }
    }
    return energy;
  }

//calculates total magnetization of system
double MainClass::calc_magnetization()
  {
    int mag = 0;

    for (int i = 0; i < L; i++)
    {
      for (int j = 0; j < L; j++)
      {
        mag += spins(i, j);
      }
    }
    return mag;
  }

//the selection of a random spin, and using the metropolis requierment for checking if we flip or not
void MainClass::Metropolis()
    {
      int delta_E;
      int x_i;
      int y_i;
      int neighbour_sum;

      for (int i = 0; i < L*L; i++)
      {
        //choose random spin
        x_i = zero_to_L_distribution(generator)%L+1; //we plus with one since we will use them on "index_vector", which is shifted by +1
        y_i = zero_to_L_distribution(generator)%L+1; //we plus with one since we will use them on "index_vector", which is shifted by +1

        //calculates the sum of the neihbours
        neighbour_sum = spins( index_vector(y_i), index_vector(x_i+1) ) + spins( index_vector(y_i), index_vector(x_i-1) ) + spins( index_vector(y_i+1), index_vector(x_i) ) + spins( index_vector(y_i-1), index_vector(x_i) );

        //uses neighbour_sum to calculate change in energy
        delta_E = 2*neighbour_sum * spins( index_vector(y_i), index_vector(x_i) ); //J=1

        //metropolis check if we want to flip or not
        if (delta_E < 0 || zero_to_one_distribution(generator) <= boltzmann_factor(delta_E + 8))
        {
          E += delta_E; //update energy
          M -= 2*spins( index_vector(y_i), index_vector(x_i));  //update magnetization
          spins( index_vector(y_i) , index_vector(x_i) ) *= -1; //flip spin
        }
      }
    }

void MainClass::Run(double T, int nr_cycles)
{
  //start by calculating the energy and magnetization of configuration
  E = calc_energy();
  M = calc_magnetization();

  //used to write the correct amount of data to file
  int denominator = nr_cycles/data_lines;

  //used to have enough digits in file name
  stringstream stream;
  stream << fixed << setprecision(3) << T;

  //open data files
  ofstream outfile_spins("../final_data/" +  filename + "_" + stream.str() + "_spins.txt");
  if (!outfile_spins.is_open())
     cout<<"Could not open file" << endl;

  ofstream outfile_expect("../final_data/" + filename + "_" + stream.str() + "_expect.txt");
  if (!outfile_expect.is_open())
     cout<<"Could not open file" << endl;

  ofstream outfile_energy("../data/" + filename + "_" + stream.str() + "_energy.txt");
  if (!outfile_energy.is_open())
     cout<<"Could not open file" << endl;

  ofstream outfile_accepts("../final_data/" + filename + "_" + stream.str() + "_accepts.txt");
  if (!outfile_accepts.is_open())
     cout<<"Could not open file" << endl;

  //reduce number of flops
  int LL = L*L;
  int LLLL = LL*LL;
  count_accepts = 0;

  double beta = 1.0/(T); //k=J=1

  //so we only need to calcualte them once
  boltzmann_factor(-8+8) = exp(8.0*beta);
  boltzmann_factor(-4+8) = exp(4.0*beta);
  boltzmann_factor(0+8)  = exp(0);
  boltzmann_factor(4+8)  = exp(-4.0*beta);
  boltzmann_factor(8+8)  = exp(-8.0*beta);

  //write spins to file before we start, so we know initial configuration
  write_spinns(outfile_spins);

  for (int cycle = 0; cycle <= nr_cycles; cycle++)
  {
    //call metropolis
    Metropolis();
    //Metropolis_Count_Accept();

    //update expectation values after metropolis
    expectationvalues(0) += E/LL;
    expectationvalues(1) += E*E/LLLL;
    expectationvalues(2) += M/LL;
    expectationvalues(3) += M*M/LLLL;
    expectationvalues(4) += fabs(M)/LL;

    //write data to file
    if ( (cycle) % denominator == 0)
    {
      cout << cycle << endl;
      outfile_energy << E << "\n";
      write_spinns(outfile_spins);
      write_expectation_values(cycle+1, outfile_expect);
      outfile_accepts << count_accepts << "  " << (cycle+1)*LL<< "\n";
    }
  }
}

//equilibrates using magnetization
  void MainClass::equilibrate_magnetization(double T)
  {
    M = calc_magnetization();

    double beta = 1.0/(T); //k=J=1
    int counter = 1;

    //pre calcualtes boltzmann factors
    boltzmann_factor(-8+8) = exp(8.0*beta);
    boltzmann_factor(-4+8) = exp(4.0*beta);
    boltzmann_factor(0+8)  = exp(0);
    boltzmann_factor(4+8)  = exp(-4.0*beta);
    boltzmann_factor(8+8)  = exp(-8.0*beta);

    //variabels used to check if we have reached equilibrium
    int old_M = M;
    int delta_M = 0;
    double tolerance = L/10.0;
    double delta_M_average;

    //continue until false
    bool not_achievd = true;
    while (not_achievd)
      {
        //run metropolis and find change in magnetization
        Metropolis();
        delta_M += abs(old_M - M);
        old_M = M;

        //every 10 runs check what the change in magnetization is
        if (counter%(10)==0)
        {
          //find average change in magnetization
          delta_M_average = delta_M/10.0;

          //check this value vs the tolerance
          if (delta_M_average < tolerance)
          {
            not_achievd = false;
          }

          //if not, set delta_M to zero
          delta_M = 0;
        }
        //increase counter by one
        counter +=1;
      }
  }


//equilibrate by simply running a given number of MC cycles
void MainClass::equilibrate(double T, int nr_cycles)
{
  double beta = 1.0/(T); //k=J=1

  //Calculate the possible Boltzmann factors for the given temperature
  boltzmann_factor(-8+8) = exp(8.0*beta);
  boltzmann_factor(-4+8) = exp(4.0*beta);
  boltzmann_factor(0+8)  = exp(0);
  boltzmann_factor(4+8)  = exp(-4.0*beta);
  boltzmann_factor(8+8)  = exp(-8.0*beta);

  for (int i = 0; i < nr_cycles/10; i++) //divide by 10
    {
      Metropolis();
    }
  }

void MainClass::reset()
{
  //Function used to reset the expectationvalues
  expectationvalues(0) = 0;
  expectationvalues(1) = 0;
  expectationvalues(2) = 0;
  expectationvalues(3) = 0;
  expectationvalues(4) = 0;
}


void MainClass::write_spinns(ofstream& OutputFile)
  {
    //Function used to write the current spin configuration to file.
    for (int i = 0; i < L; i++)
    {
      for (int j = 0; j < L; j++)
      {
        OutputFile << spins(i, j) << " ";
      }
    }
    OutputFile << " \n";
  }


void MainClass::write_expectation_values(int cycle_nr, ofstream& OutputFile)
  {
    //Function used to write expectation values to file.
    vec dummy_expectationvalues = expectationvalues/(cycle_nr);
    OutputFile << cycle_nr << " " << dummy_expectationvalues(0) << " " << dummy_expectationvalues(1) << " " << dummy_expectationvalues(2) << " " << dummy_expectationvalues(3) << " " << dummy_expectationvalues(4) << "\n";
  }


//this is identical to the metropolis function, but this one counts how many accepts we have
//if unsure about a step, see metropolis function
void MainClass::Metropolis_Count_Accept()
    {
      //Function used to count the number of accepted flips.

      int delta_E;
      int x_i;
      int y_i;
      int neighbour_sum;

      for (int i = 0; i < L*L; i++)
      {
        //choose random spin
        //use this since the other one does not seem to work, there is an issue at github
        x_i = zero_to_L_distribution(generator)%L+1; //we plus with one since we will use them on "index_vector", which is shifted by +1
        y_i = zero_to_L_distribution(generator)%L+1; //we plus with one since we will use them on "index_vector", which is shifted by +1

        neighbour_sum = spins( index_vector(y_i), index_vector(x_i+1) ) + spins( index_vector(y_i), index_vector(x_i-1) ) + spins( index_vector(y_i+1), index_vector(x_i) ) + spins( index_vector(y_i-1), index_vector(x_i) );

        delta_E = 2*neighbour_sum * spins( index_vector(y_i), index_vector(x_i) ); //J=1

        if (delta_E < 0 || zero_to_one_distribution(generator) <= boltzmann_factor(delta_E + 8))
        {
          count_accepts += 1;
          E += delta_E;
          M -= 2*spins( index_vector(y_i), index_vector(x_i));
          spins( index_vector(y_i) , index_vector(x_i) ) *= -1;
        }
      }

    }

//times the function
void MainClass::Timer(string start_or_stop, string filename)
  {
    //if you wrote "start" we start the clock
    if (start_or_stop == "start")
      start = clock(); //start clock

      //if you wrote "stop" we stop the clock, and write to file
    if (start_or_stop == "stop")
    {
      stop = clock();//finish timer
      time_used = (stop-start)/CLOCKS_PER_SEC; //calculate time used

      //writes the time to file, with the lattice size
      std::ofstream outfile;
      outfile.open("../final_data/"+ filename +"_timer.txt", std::ios_base::app);
      if (!outfile.is_open())
        cout << "Could not open file" << endl;
      else
        outfile << setw(15) << L << "      " << time_used << endl;
      }
  }
