#include "ising1D.h"

int Ising1D :: Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= m_nspin) i = i - m_nspin;
    else if(i < 0) i = i + m_nspin;
    return i;
}

//Returns the boltzmann weight for a given energy
double Ising1D :: BoltzmannWeight(double energy){
    return exp(-m_beta*energy);
}

//Returns the energy difference when the spin in position k flips
double Ising1D :: EnergyDifference(int k){
  double delta_energy = 2*m_J*m_s[k]*( m_s[Pbc(k-1)] + m_s[Pbc(k+1)] ) + 2*m_h*m_s[k];
  return delta_energy;
}
//Returns the energy difference when the spin in position k flips if the k-th spin is sk
double Ising1D :: EnergyDifference(int k, int sk){
  double delta_energy = 2*m_J*sk*( m_s[Pbc(k-1)] + m_s[Pbc(k+1)] ) + 2*m_h*sk;
  return delta_energy;
}

double Ising1D ::  TotalEnergy(){
  double U=0;
  for(int i=0; i<m_nspin; i++){
    U+=(-m_J * m_s[i] * m_s[Pbc(i+1)] - 0.5 * m_h * (m_s[i] + m_s[Pbc(i+1)]));
  }
  return U;
}

double Ising1D :: TotalSpin(){
  double M=0.;
  for(int i=0; i<m_nspin; i++){
    M+=m_s[i];
  }
  return M;
}

Ising1D :: Ising1D(Random *rnd){
  m_rnd=rnd;

  m_block_U = 0.;
  m_block_M = 0.;
  m_block_U2 = 0.;
  m_block_M2 = 0.;

  m_ave_U = 0.;
  m_ave_M = 0.;
  m_ave_chi = 0.;
  m_ave_C = 0.;
  m_ave2_U = 0.;
  m_ave2_M = 0.;
  m_ave2_chi = 0.;
  m_ave2_C = 0.;

  m_block_index = 0;
  m_block_number = 0;

  m_attempted=0;
  m_accepted=0;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl;
  cout << "Nearest neighbour interaction      " << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl<<endl;

  ifstream ReadInput;
  //Read input data for the simulation
  ReadInput.open("input.dat");
  if(ReadInput.is_open()){
    ReadInput >> m_temp; //Temperature
    ReadInput >> m_nspin; //Number of spins
    ReadInput >> m_J; //Exchange interaction parameter
    ReadInput >> m_h; //External field
    ReadInput >> m_metro; //1 to use metropolis sampling, 0 to use gibbs sampling
    ReadInput >> m_nblk; //Number of blocks
    ReadInput >> m_nstep; //Number of steps in each block
    ReadInput >> m_restart; //1 to restart from old config, 0 otherwise
    ReadInput >> m_printall; //1 to print every value of U and M calculated, 0 otherwise
    ReadInput >> m_eq_nstep; //Number of step for equilibration
  }else cerr<<"Unable to open input.dat"<<endl;

  cout<<"Input data from input.dat"<<endl;
  cout<<"Temperature: "<<m_temp<<endl;
  cout<<"Number of spins: "<<m_nspin<<endl;
  cout<<"Exchange interaction: "<<m_J<<endl;
  cout<<"External field: "<<m_h<<endl;
  cout<<"Number of blocks: "<<m_nblk<<endl;
  cout<<"Number of steps in each block: "<<m_nstep<<endl;
  if(m_metro==1) cout << "Metropolis moves" << endl;
  else cout << "Gibbs moves" << endl;
  cout<<"Starting from old configuration: "<<m_restart<<endl;

  m_beta=1./m_temp;

  ReadInput.close();

  m_s = vector<int>(m_nspin);
  if(m_restart){
    RestartInitialization();
  }else{
    FirstInitialization();
  }
}

void Ising1D :: RestartInitialization(){
  //Read configuration r(t)
  ifstream ReadConf;
  cout<<"Reading initial configuration from config.final"<<endl;
  ReadConf.open("config.final");
  if(ReadConf.is_open()){
    for(int i=0; i<m_nspin; i++){
      ReadConf >> m_s[i];
    }
  }else cerr<<"Unable to open config.final"<<endl;
  ReadConf.close();
  cout<<"Energy per spin from config.final: " << TotalEnergy()/double(m_nspin) << endl;
}

void Ising1D :: FirstInitialization(){
  //Initialize spins randomly (i.e. T=\infty)
  for (int i=0; i<m_nspin; ++i)
  {
    if(m_rnd->Rannyu() >= 0.5) m_s[i] = 1;
    else m_s[i] = -1;
  }
  cout << "Energy per spin before equilibration: " << TotalEnergy()/double(m_nspin) << endl;
  for(int i=0; i<m_eq_nstep;i++){
    Move();
  }
  m_attempted=0;
  m_accepted=0;
  cout << "Energy per spin after equilibration: " << TotalEnergy()/double(m_nspin) << endl;
}


void Ising1D :: Move(){
  if(m_metro == 1){
    for(int i=0; i<m_nspin; i++){
      //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
      int k = (int)(m_rnd->Rannyu()*m_nspin);
      MoveMetropolis(k);
    }
  }else{
    for(int i=0; i<m_nspin; i++){
      MoveGibbs(i);
    }
  }
}


double Ising1D :: AcceptanceProbability(int k){
  double alpha = fmin(1, BoltzmannWeight(EnergyDifference(k)));
  m_alpha=alpha;
  return alpha;
}

void Ising1D :: MoveMetropolis(int k){
  m_attempted++;
  double alpha = AcceptanceProbability(k);
  double r = m_rnd -> Rannyu(); //Uniform in [0,1) to test acceptance
  if(r < alpha){
    m_s[k]=-m_s[k]; //Flip the spin k
    m_accepted++;
  } // else keep the same configuration
}

double Ising1D :: ConditionalProbability(int k, int sk){
  return 1./(1.+BoltzmannWeight(EnergyDifference(k,sk)));
}

void Ising1D :: MoveGibbs(int k){
  m_attempted++;
  m_accepted++;
  double p1 = ConditionalProbability(k,1);
  double r = m_rnd -> Rannyu(); //Uniform in [0,1) to sample conditional probability
  if(r < p1){
    m_s[k]=1;
  }else{
    m_s[k]=-1;
  }
}


void Ising1D :: Measure(){

  double U=0., M=0.;

  U=TotalEnergy();
  M=TotalSpin();

  m_block_U=m_block_index/double(m_block_index+1)*m_block_U+1./double(m_block_index+1)*U;
  m_block_U2=m_block_index/double(m_block_index+1)*m_block_U2+1./double(m_block_index+1)*U*U;
  m_block_M=m_block_index/double(m_block_index+1)*m_block_M+1./double(m_block_index+1)*M;
  m_block_M2=m_block_index/double(m_block_index+1)*m_block_M2+1./double(m_block_index+1)*M*M;

  m_block_index++;
  if(m_printall){
    ofstream Ustream, Mstream;

    Ustream.open("output_U.dat",ios::app); //ios::app appends at the end of the file
    Mstream.open("output_M.dat",ios::app);

    if(Ustream.is_open()){
      Ustream<<U/double(m_nspin)<<endl;
    }else cerr<<"Unable to open output_U.dat"<<endl;
    if(Mstream.is_open()){
      Mstream<<M/double(m_nspin)<<endl;
    }else cerr<<"Unable to open Mstream.dat"<<endl;


    Ustream.close();
    Mstream.close();
  }
}

void Ising1D :: ConfFinal(){
  ofstream WriteConf;
  cout<<"Printing final configuration to file config.final "<<endl;
  WriteConf.open("config.final");

  if(WriteConf.is_open()){
    for(int i=0; i<m_nspin; i++){
      WriteConf<<m_s[i]<<endl;
    }
  }else cerr<<"Unable to open config.final"<<endl;

  WriteConf.close();
}

void Ising1D :: MeasureBlock(){
  double error_U, error_M, error_C, error_chi;

  m_ave_U = m_block_number/double(m_block_number+1)*m_ave_U+1./double(m_block_number+1)*m_block_U;
  m_ave2_U = m_block_number/double(m_block_number+1)*m_ave2_U+1./double(m_block_number+1)*m_block_U*m_block_U;
  m_ave_M = m_block_number/double(m_block_number+1)*m_ave_M+1./double(m_block_number+1)*m_block_M;
  m_ave2_M = m_block_number/double(m_block_number+1)*m_ave2_M+1./double(m_block_number+1)*m_block_M*m_block_M;

  double C=m_beta*m_beta*(m_block_U2-m_block_U*m_block_U); //Value for the current block
  double chi=0;
  if(m_h == 0){
    //Assume <M>=0
    chi=m_beta*m_block_M2; //Value for the current block
  }else{
    chi=m_beta*(m_block_M2-m_block_M*m_block_M); //Value for the current block
  }

  m_ave_C = m_block_number/double(m_block_number+1)*m_ave_C+1./double(m_block_number+1)*C;
  m_ave2_C = m_block_number/double(m_block_number+1)*m_ave2_C+1./double(m_block_number+1)*C*C;
  m_ave_chi = m_block_number/double(m_block_number+1)*m_ave_chi+1./double(m_block_number+1)*chi;
  m_ave2_chi = m_block_number/double(m_block_number+1)*m_ave2_chi+1./double(m_block_number+1)*chi*chi;


  if(m_block_number == 0){
    error_U=0;
    error_M=0;
    error_C=0;
    error_chi=0;
  }else{
    error_U=sqrt((m_ave2_U-m_ave_U*m_ave_U)/m_block_number);
    error_M=sqrt((m_ave2_M-m_ave_M*m_ave_M)/m_block_number);
    error_C=sqrt((m_ave2_C-m_ave_C*m_ave_C)/m_block_number);
    error_chi=sqrt((m_ave2_chi-m_ave_chi*m_ave_chi)/m_block_number);
  }

  m_block_number++;


  //Print the averages calculated using m_block_number blocks
  ofstream Ustream, Mstream, Cstream, chistream;

  Ustream.open("ave_U.out",ios::app); //ios::app appends at the end of the file
  Mstream.open("ave_M.out",ios::app);
  Cstream.open("ave_C.out",ios::app);
  chistream.open("ave_chi.out",ios::app);

  if(Ustream.is_open()){
    Ustream<<m_block_number<<",\t"<<m_ave_U/double(m_nspin)<<",\t"<<error_U/double(m_nspin)<<endl;
  }else cerr<<"Unable to open ave_U.out"<<endl;
  if(Mstream.is_open()){
    Mstream<<m_block_number<<",\t"<<m_ave_M/double(m_nspin)<<",\t"<<error_M/double(m_nspin)<<endl;
  }else cerr<<"Unable to open ave_M.out"<<endl;
  if(Cstream.is_open()){
    Cstream<<m_block_number<<",\t"<<m_ave_C/double(m_nspin)<<",\t"<<error_C/double(m_nspin)<<endl;
  }else cerr<<"Unable to open ave_C.out"<<endl;
  if(chistream.is_open()){
    chistream<<m_block_number<<",\t"<<m_ave_chi/double(m_nspin)<<",\t"<<error_chi/double(m_nspin)<<endl;
  }else cerr<<"Unable to open ave_chi.out"<<endl;

  Ustream.close();
  Mstream.close();
  Cstream.close();
  chistream.close();

  //Reset for another block
  m_block_U = 0.;
  m_block_U2 = 0.;
  m_block_M = 0.;
  m_block_M2 = 0.;
  m_block_index = 0;

  //Print the temperature, the final block number and the final values with their errors
  if(m_block_number == m_nblk){
    const int wd=12;
    ofstream final_results;
    final_results.open("final.out",ios::app);
    if(final_results.is_open()){
      final_results<<m_temp<<","<<setw(wd)<<m_h<<","<<setw(wd)<<double(m_accepted)/double(m_attempted)<<",";
      final_results<<setw(wd)<<m_ave_U/double(m_nspin)<<","<<setw(wd)<<error_U/double(m_nspin)<<",";
      final_results<<setw(wd)<<m_ave_M/double(m_nspin)<<","<<setw(wd)<<error_M/double(m_nspin)<<",";
      final_results<<setw(wd)<<m_ave_C/double(m_nspin)<<","<<setw(wd)<<error_C/double(m_nspin)<<",";
      final_results<<setw(wd)<<m_ave_chi/double(m_nspin)<<","<<setw(wd)<<error_chi/double(m_nspin)<<","<<endl;
    }else cerr<<"Unable to open final.out"<<endl;
    final_results.close();
  }
}
