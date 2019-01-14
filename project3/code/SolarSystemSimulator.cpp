#include "SolarSystemSimulator.hpp"

SolarSystemSimulator::SolarSystemSimulator()
{}

SolarSystemSimulator::SolarSystemSimulator(vec mass, mat pos, mat vel,\
                    long long n, double T_max, string make_filename, bool E)
{
  G = 4*pi*pi;
  m = mass;
  //t = vec(n);
  lines_of_data = n/1000;

  // create time array
  filename = make_filename;
  energy = E;

  integrationpoints = n;
  dt =  T_max/(n-1);

  //for current time step
  x1 = mat(m.n_elem, 3); //Planets, time, space-coordinate
  v1 = mat(m.n_elem, 3); //Planets, time, space-coordinate

  //for next time step
  x2 = mat(m.n_elem, 3); //Planets, time, space-coordinate
  v2 = mat(m.n_elem, 3); //Planets, time, space-coordinate

  x1 = pos;     //Initial position
  v1 = vel;     //Initial velocity
}

mat SolarSystemSimulator::gravity(mat x) //j = planet numner, i = time index
{
  mat acc = mat(m.n_elem, 3, fill::zeros);
  mat forcevec = mat(1, 3);
  double force;
  for (int k = 0; k < m.n_elem; k++)
  {
    for (int j = 0; j < k; j++)
    {
      //sums over forces acting on each planet and adds them to the matrix
      force = m(k)*m(j)*pow(length((x1.row(j)-x1.row(k)).t()), -3.0);

      forcevec = force*(x.row(j) - x.row(k));
      acc.row(k) += forcevec/m(k);
      acc.row(j) -= forcevec/m(j);
    }
  }

  return acc*G; //return force matrix, now multiplied by G
}

double SolarSystemSimulator::AngularMomentum()
{
  //calculates total angular momentum of system through definition of angular momentum
  double L = 0;
  double l_x = 0;
  double l_y = 0;
  double l_z = 0;
  for (int j = 0; j < m.n_elem; j++)
  {
    l_x += m(j)*(x1(j,1)*v1(j,2) - x1(j,2)*v1(j,1));
    l_y += m(j)*(x1(j,2)*v1(j,0) - x1(j,0)*v1(j,2));
    l_z += m(j)*(x1(j,0)*v1(j,1) - x1(j,1)*v1(j,0));
  }
  L += sqrt(l_x*l_x + l_y*l_y + l_z*l_z);
  return L;
}

double SolarSystemSimulator::KineticEnergy()
{
  //calculates total kinetic energy of system through the definition of kinetic energy
  double K = 0;
  for (int i = 0; i < m.n_elem; i++)
  {
    K += 0.5*m(i)*(v1(i, 0)*v1(i, 0) + v1(i, 1)*v1(i, 1) + v1(i, 2)*v1(i, 2));
  }
  return K;
}

double SolarSystemSimulator::PotentialEnergy()
{
  //calculates total potential energy of system through definition of gravitational potential energy
  double V = 0;
  double distance;

  for (int i = 0; i < m.n_elem; i++)
  {
    for (int j = 0; j < i; j++)
    {
      distance = length((x1.row(i)-x1.row(j)).t());
      V += -G*m(i)*m(j)/distance;
    }
  }
  return V;
}

void SolarSystemSimulator::TotalEnergy(string filename)
{
  //calcualtes total energy and total angular momentum and writes the data to file, equal to variable filename
  mat energy_and_mom = mat(4,1);

  energy_and_mom(0,0) = KineticEnergy() + PotentialEnergy();
  energy_and_mom(1,0) = KineticEnergy();
  energy_and_mom(2,0) = PotentialEnergy();
  energy_and_mom(3,0) = AngularMomentum();

  WriteData(filename, energy_and_mom, "energy");
}


void SolarSystemSimulator::ChangeFileName(string new_filename)
{
  //if we want to change the name of data file
  filename = new_filename;
}

void SolarSystemSimulator::ResetInitialConditions(vec masses, mat pos, mat vel)
{
  //in case we want to change initial conditions
  m = masses;
  x1 = pos;
  v1 = vel;
}

//run forward Euler
void SolarSystemSimulator::ForwardEuler()
{
  mat a;

  //gets ready to write data to file, by altering filename and deleting existing datafiles
  filename += "_FE_";
  DeleteDatafile(filename, "coordinates");
  DeleteDatafile(filename + "energy",   "energy");

  for (long long i = 0; i < integrationpoints; i++)
  {
    a = gravity(x1);
    v2 = v1 + dt*a;
    x2 = x1 + dt*v1;

    //if we want to save positions -> save positions
    if (i%lines_of_data == 0)
    {
      WriteData(filename, x1, "coordinates");
      if (energy) //do we want energy as well?
        TotalEnergy(filename + "energy");
    }
    // update new and old positions and velocities
    x1 = x2;
    v1 = v2;
  }
  //always write last position to file, for analyzing
  WriteData(filename, x1, "coordinates");
  if (energy)
  TotalEnergy(filename + "energy");
}

void SolarSystemSimulator::VelocityVerlet()
{
  mat ai;
  mat aii;

  //gets ready to save data
  filename += "_VV_";
  DeleteDatafile(filename, "coordinates");
  if (energy) DeleteDatafile(filename + "energy", "energy");

  for (long long i = 0; i < integrationpoints; i++)
  {
    ai = gravity(x1);
    x2 = x1 + dt*v1 + dt*dt/2*ai;

    aii = gravity(x2);
    v2 = v1 + dt/2*(aii + ai);

    //if we want to save positions -> save positions
    if (i%lines_of_data == 0)
    {
      WriteData(filename, x1, "coordinates");
      if (energy) //save energy as well?
        TotalEnergy(filename + "energy");
    }
    //update positions
    x1 = x2;
    v1 = v2;
  }
  //always write last position and energy to file for analyzing
  WriteData(filename, x1, "coordinates");
  if (energy)
  TotalEnergy(filename + "energy");
}

//calcualte length of vector
double SolarSystemSimulator::length(vec r)
{
  double ans = 0;
  for (int i = 0; i < r.n_elem; i++)
  {
    ans += r(i)*r(i);
  }
  return sqrt(ans);
}

//write a matrix of data to file
void SolarSystemSimulator::WriteData(string filename, mat matrix, string type)
{
  if (type == "coordinates")
  {
    std::ofstream outfilex;
    std::ofstream outfiley;
    std::ofstream outfilez;

    outfilex.open("data/"+ filename +"_x.txt", std::ios_base::app);
    outfiley.open("data/"+ filename +"_y.txt", std::ios_base::app);
    outfilez.open("data/"+ filename +"_z.txt", std::ios_base::app);

    if (!outfilex.is_open() || !outfiley.is_open() || !outfilez.is_open())
      cout << "Could not open file" << endl;
    else
      for (int l = 0; l < matrix.n_rows; l++)
      {
        outfilex << setw(15) << matrix(l, 0);
        outfiley << setw(15) << matrix(l, 1);
        outfilez << setw(15) << matrix(l, 2);
      }
      outfilex << endl;
      outfiley << endl;
      outfilez << endl;
  }

  else if (type == "energy")
  {
    std::ofstream outfile;
    outfile.open("data/"+ filename +".txt", std::ios_base::app);

    if (!outfile.is_open())
      cout << "Could not open file" << endl;
    else
    {
      outfile << setw(30) << setprecision(20) << matrix(0,0) << setw(30) << setprecision(20) << matrix(1,0) << setw(30) << setprecision(20) << matrix(2,0) << setw(30) << setprecision(20) << matrix(3,0);
      outfile << endl;
    }
  }

}

//overwrite an existing data file
void SolarSystemSimulator::DeleteDatafile(string filename, string type)
{
  if (type == "coordinates")
  {
    ofstream outfilex("data/"+ filename +"_x.txt");
    ofstream outfiley("data/"+ filename +"_y.txt");
    ofstream outfilez("data/"+ filename +"_z.txt");
    if (!outfilex.is_open() || !outfiley.is_open() || !outfilez.is_open())
      cout<<"Could not open file" << endl;
  }
  else if (type == "energy")
  {
    ofstream outfile("data/"+ filename +".txt");
    if (!outfile.is_open())
      cout<<"Could not open file" << endl;
  }
}

//times algorithm
void SolarSystemSimulator::Timer(string start_or_stop, string filename)
{
  if (start_or_stop == "start")
    start = clock(); //start clock

  if (start_or_stop == "stop")
  {
    stop = clock();//finish timer
    time_used = (stop-start)/CLOCKS_PER_SEC; //calculate time used

    std::ofstream outfile;
    outfile.open("data/"+ filename +"_timer.txt", std::ios_base::app);
    if (!outfile.is_open())
      cout << "Could not open file" << endl;
    else
      outfile << setw(15) << integrationpoints << "      " << time_used << endl;
    }
}
