#include "MonteCarloNVT.h"
//Returns the boltzmann weight for a given energy
double MonteCarloNVT :: BoltzmannWeight(double energy){
    return exp(-m_beta*energy);
}

//Returns the interaction energy between particle ip and all the others if particle ip is in position (x,y,z)
double MonteCarloNVT :: InteractionEnergy(double x, double y, double z, int ip){
  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<m_npart; i++){
    if(i != ip){
      // distance ip-i in pbc
      dx = Pbc(x - m_x[i]);
      dy = Pbc(y - m_y[i]);
      dz = Pbc(z - m_z[i]);

      dr = sqrt(dx*dx + dy*dy + dz*dz);
      if(dr < m_rcut)
      {
        double dr6=pow(dr,6);
        double dr12=dr6*dr6;
        ene += 1.0/dr12 - 1.0/dr6;
      }
    }
  }

  return 4.0*ene;
}

MonteCarloNVT :: MonteCarloNVT(Random *rnd){
  m_rnd=rnd;

  cout << "Classic Lennard-Jones fluid             " << endl;
  cout << "Monte Carlo simulation             " << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl;
  cout << "The program uses Lennard-Jones units " << endl<<endl;

  ifstream ReadInput;
  //Read input data for the simulation
  ReadInput.open("input.dat");
  if(ReadInput.is_open()){
    ReadInput >> m_temp; //Temperature
    ReadInput >> m_npart; //Number of particles
    ReadInput >> m_rho; //Density
    ReadInput >> m_rcut; //Cutoff radius for potential evaluation
    ReadInput >> m_delta; //Width for uniform sampling [x-delta/2, x+delta/2)
    ReadInput >> m_nblk; //Number of blocks
    ReadInput >> m_nstep; //Number of steps in each block
    ReadInput >> m_restart; //1 to restart from old config, 0 otherwise
    ReadInput >> m_printall; //1 to print every value of U and M calculated, 0 otherwise
    ReadInput >> m_eq_nstep; //Number of step for equilibration
    ReadInput >> m_nbins; //Number of bins for g(r) histogram
  }else cerr<<"ERROR: Unable to open input.dat"<<endl;

  m_beta=1./m_temp;
  m_vol = (double)m_npart/m_rho; // Volume
  m_box = pow(m_vol,1./3.); // Cube box, linear dimension
  ReadInput.close();

  cout<<"Input data from input.dat (LJ units)"<<endl;
  cout<<"Temperature: "<<m_temp<<endl;
  cout<<"Number of particles: "<<m_npart<<endl;
  cout<<"Density: "<<m_rho<<endl;
  cout<<"Volume of the simulation box: "<<m_vol<<endl;
  cout<<"Edge of the simulation box: "<<m_box<<endl;
  cout<<"Cutoff radius: "<<m_rcut<<endl;
  cout<<"Number of blocks: "<<m_nblk<<endl;
  cout<<"Number of steps in each block: "<<m_nstep<<endl;
  cout<<"Starting from old configuration: "<<m_restart<<endl;
  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter: " << m_delta << endl;
  cout<<"Starting from old configuration: "<<m_restart<<endl;

  ReadInput.close();

  m_x = vector<double>(m_npart);
  m_y = vector<double>(m_npart);
  m_z = vector<double>(m_npart);

  m_block_v = 0.;
  m_block_w = 0.;
  m_ave_v = 0.;
  m_ave_w = 0.;
  m_ave2_v = 0.;
  m_ave2_w = 0.;
  m_block_g = vector<double>(m_nbins);
  m_ave_g = vector<double>(m_nbins);
  m_ave2_g = vector<double>(m_nbins);
  m_g = vector<double>(m_nbins);
  m_error_g = vector<double>(m_nbins);
  m_histo = new Histogram(m_nbins,0,m_box/2.);
  cout<<"Number of bins: "<<m_nbins<<endl;
  cout<<"Bin length: "<<m_histo->getBinLength()<<endl;
  m_block_index = 0;
  m_block_number = 0;

  m_attempted=0;
  m_accepted=0;

  //Tail corrections (per particle) for potential energy and virial
  m_vtail = (8.0*pi*m_rho)/(9.0*pow(m_rcut,9)) - (8.0*pi*m_rho)/(3.0*pow(m_rcut,3));
  m_wtail = (32.0*pi*m_rho)/(9.0*pow(m_rcut,9)) - (16.0*pi*m_rho)/(3.0*pow(m_rcut,3));
  cout <<endl<< "Tail correction for the potential energy: " << m_vtail << endl;
  cout << "Tail correction for the virial: " << m_wtail << endl;
  cout <<"Tail correction for pressure: " <<m_wtail*m_npart/m_vol << endl;;

  if(m_restart){
    RestartInitialization();
  }else{
    FirstInitialization();
  }

}

void MonteCarloNVT :: RestartInitialization(){
  //Read configuration r(t)
  ifstream ReadConf;
  cout<<"Reading initial configuration from config.final"<<endl;
  ReadConf.open("config.final");
  if(ReadConf.is_open()){
    for(int i=0; i<m_npart; i++){
      ReadConf >> m_x[i] >> m_y[i] >> m_z[i];
      m_x[i] = m_x[i] * m_box;
      m_y[i] = m_y[i] * m_box;
      m_z[i] = m_z[i] * m_box;
    }
  }else cerr<<"Unable to open config.final"<<endl;
  ReadConf.close();
  cout<<"No equilibration done"<<endl;
  Measure();
  cout<<"Initial potential energy (with tail corrections): " <<m_block_v/double(m_npart)+m_vtail<<endl;
  cout <<"Initial virial (with tail corrections): " <<m_block_w/double(m_npart)+m_wtail << endl;
  cout <<"Initial pressure (with tail corrections): " <<m_rho*m_temp+(m_block_w+m_wtail*m_npart)/m_vol << endl<<endl;
}

void MonteCarloNVT :: FirstInitialization(){
  ifstream ReadConf;
  cout<<"Reading initial configuration from config.0"<<endl;
  ReadConf.open("config.0");
  if(ReadConf.is_open()){
    for(int i=0; i<m_npart; i++){
      ReadConf >> m_x[i] >> m_y[i] >> m_z[i];
      m_x[i] = m_x[i] * m_box;
      m_y[i] = m_y[i] * m_box;
      m_z[i] = m_z[i] * m_box;
    }
  }else cerr<<"Unable to open config.0"<<endl;
  ReadConf.close();

  cout<<"Equilibration with " <<m_eq_nstep<<" steps"<<endl;
  for(int i=0; i<m_eq_nstep;i++){
    Move();
  }
  //Reset counter after equilibration
  m_accepted=0;
  m_attempted=0;

  Measure();
  cout<<"Initial potential energy (with tail corrections): " <<m_block_v/double(m_npart)+m_vtail<<endl;
  cout <<"Initial virial (with tail corrections): " <<m_block_w/double(m_npart)+m_wtail << endl;
  cout <<"Initial pressure (with tail corrections): " <<m_rho*m_temp+(m_block_w+m_wtail*m_npart)/m_vol << endl << endl;;

}


void MonteCarloNVT :: Move(){
  int o;
  double energy_old, energy_new;
  double xold, yold, zold, xnew, ynew, znew;

  for(int i=0; i<m_npart; ++i){
    //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
    o = (int)(m_rnd->Rannyu()*m_npart);

    //Old
    xold = m_x[o];
    yold = m_y[o];
    zold = m_z[o];

    energy_old = InteractionEnergy(xold,yold,zold,o);

    //New
    xnew = Pbc( m_x[o] + m_delta*(m_rnd->Rannyu() - 0.5) );
    ynew = Pbc( m_y[o] + m_delta*(m_rnd->Rannyu() - 0.5) );
    znew = Pbc( m_z[o] + m_delta*(m_rnd->Rannyu() - 0.5) );

    energy_new = InteractionEnergy(xnew,ynew,znew,o);

    //Metropolis test
    double p = BoltzmannWeight(energy_new-energy_old); //no need to calculate acceptance rate alpha

    if(p >= m_rnd->Rannyu()){
      //Update
      m_x[o] = xnew;
      m_y[o] = ynew;
      m_z[o] = znew;

      m_accepted++;
    }
    m_attempted++;
  }
}

void MonteCarloNVT :: Measure(){
  double v=0, w=0;
  double vij=0., wij=0.;
  double dx,dy,dz,dr;
  m_histo->reset(); //set all the bins to zero

  for(int i=0; i<m_npart-1; i++){
    for(int j=i+1; j<m_npart; j++){
      dx=Pbc(m_x[i]-m_x[j]);
      dy=Pbc(m_y[i]-m_y[j]);
      dz=Pbc(m_z[i]-m_z[j]);
      dr=sqrt(dx*dx+dy*dy+dz*dz);
      //update histogram g(r)
      m_histo->fill(dr);
      m_histo->fill(dr); //increase bin by 2

      if(dr < m_rcut){
        double dr6=pow(dr,6);
        double dr12=dr6*dr6;
        vij=1./dr12-1./dr6;
        wij=1./dr12-0.5/dr6;
        v+=vij;
        w+=wij;
      }
    }
  }
  v=4.*v;
  w=48.*w/3.;

  double bin_length = m_histo->getBinLength();
  for(int i=0; i<m_histo->getNumberOfBins(); i++){
    double norm = m_rho*m_npart*4./3.*pi*(pow((i+1)*bin_length,3)-pow(i*bin_length,3));
    m_g[i] = m_histo->getBin(i)/norm;
  }

  //add new value to average of block
  m_block_v=m_block_index/double(m_block_index+1)*m_block_v+1./double(m_block_index+1)*v;
  m_block_w=m_block_index/double(m_block_index+1)*m_block_w+1./double(m_block_index+1)*w;
  for(int i=0; i<m_nbins; i++){
    m_block_g[i]=m_block_index/double(m_block_index+1)*m_block_g[i]+1./double(m_block_index+1)*m_g[i];
  }

  m_block_index++;
  if(m_printall){
    ofstream Ustream, Pstream;

    Ustream.open("output_Uinst.dat",ios::app); //ios::app appends at the end of the file
    Pstream.open("output_Pinst.dat",ios::app);
    //Print with tail corrections
    if(Ustream.is_open()){
      Ustream<<v/double(m_npart)+m_vtail<<endl;
    }else cerr<<"Unable to open output_Uinst.dat"<<endl;
    if(Pstream.is_open()){
      double P=m_rho*m_temp+(w+m_wtail*m_npart)/m_vol; //Print pressure instead of virial
      Pstream<<P<<endl;
    }else cerr<<"Unable to open output_Pinst.dat"<<endl;

    Ustream.close();
    Pstream.close();
  }
}

void MonteCarloNVT :: ConfXYZ(int nconf){
  ofstream WriteXYZ;
  WriteXYZ.open("frames/config_"+to_string(nconf)+".xyz");

  if(WriteXYZ.is_open()){
    WriteXYZ<<m_npart<<endl;
    WriteXYZ<<"XYZ line for comment"<<endl;
    //<Atom_species> <x> <y> <z> is the XYZ format
    //All our atoms are the same, so we write a generic "LJ" instead of a specific element
    for(int i=0; i<m_npart; i++){
      WriteXYZ<<"LJ "<<Pbc(m_x[i])<<" "
              <<" "<<Pbc(m_y[i])<<" "
              <<" "<<Pbc(m_z[i])<<endl;
    }
  }else cerr<<"Unable to open frames/config_"+to_string(nconf)+".xyz"<<endl;
  WriteXYZ.close();
}

void MonteCarloNVT :: ConfFinal(){
  ofstream WriteConf;
  cout<<"Printing final configuration to file config.final "<<endl;
  WriteConf.open("config.final");

  if(WriteConf.is_open()){
    for(int i=0; i<m_npart; i++){
      WriteConf<<m_x[i]/m_box<<" "<<m_y[i]/m_box<<" "<<m_z[i]/m_box<<endl;
    }
  }else cerr<<"Unable to open config.final"<<endl;
  WriteConf.close();
}

void MonteCarloNVT :: MeasureBlock(){
  double error_v, error_w;

  //Aggiorno le medie progressive aggiungendo un nuovo blocco
  m_ave_v = m_block_number/double(m_block_number+1)*m_ave_v+1./double(m_block_number+1)*m_block_v;
  m_ave2_v = m_block_number/double(m_block_number+1)*m_ave2_v+1./double(m_block_number+1)*m_block_v*m_block_v;
  m_ave_w = m_block_number/double(m_block_number+1)*m_ave_w+1./double(m_block_number+1)*m_block_w;
  m_ave2_w = m_block_number/double(m_block_number+1)*m_ave2_w+1./double(m_block_number+1)*m_block_w*m_block_w;
  for(int i=0; i<m_nbins; i++){
    m_ave_g[i]=m_block_number/double(m_block_number+1)*m_ave_g[i]+1./double(m_block_number+1)*m_block_g[i];
    m_ave2_g[i]=m_block_number/double(m_block_number+1)*m_ave2_g[i]+1./double(m_block_number+1)*m_block_g[i]*m_block_g[i];
  }

  //Calcolo errori considerando tutti i blocchi fin qui misurati
  if(m_block_number == 0){
    error_v=0;
    error_w=0;
    for(auto& error_gi : m_error_g) error_gi=0;
  }else{
    error_v=sqrt((m_ave2_v-m_ave_v*m_ave_v)/m_block_number);
    error_w=sqrt((m_ave2_w-m_ave_w*m_ave_w)/m_block_number);
    for(int i=0; i<m_nbins; i++){
      m_error_g[i]=sqrt((m_ave2_g[i]-m_ave_g[i]*m_ave_g[i])/m_block_number);
    }
  }
  //Aggiorno il numero di blocchi dopo aver calcolato l'errore (quindi al denominatore utilizzo N-1)
  m_block_number++;
  double p=m_rho*m_temp+(m_ave_w+m_wtail*m_npart)/m_vol;
  double error_p=error_w/m_vol;

  ofstream vstream, pstream, gstream;
  vstream.open("output.uave.0",ios::app); //ios::app appends at the end of the file
  pstream.open("output.pave.0",ios::app);
  gstream.open("output.gofr.0",ios::app);

  //Stampo le medie di v e p considerando m_block_number blocchi
  if(vstream.is_open()){
    vstream<<m_block_number<<",\t"<<m_ave_v/double(m_npart)+m_vtail<<",\t"<<error_v/double(m_npart)<<endl;
  }else cerr<<"Unable to open output.uave.0"<<endl;
  if(pstream.is_open()){
    pstream<<m_block_number<<",\t"<<p<<",\t"<<error_p<<endl;
  }else cerr<<"Unable to open output.pave.0"<<endl;

  //Stampo il valore di g(r) blocco per blocco (non facendo la media)
  if(gstream.is_open()){
    double bin_length=m_histo->getBinLength();
    gstream<<"===== BLOCK NUMBER "<<m_block_number<<" ====="<<endl;
    for(int i=0; i<m_nbins; i++){
      gstream<<m_block_number<<",\t"<<i*bin_length<<",\t"<<(i+1)*bin_length<<",\t"<<m_block_g[i]<<",\t"<<endl;
    }
  }else cerr<<"Unable to open output.gofr.0"<<endl;

  vstream.close();
  pstream.close();
  gstream.close();

  cout << "Block number " << m_block_number << endl;
  cout << "Acceptance rate " << getAcceptanceRate() << endl;
  cout << "-------------------------"<<endl<<endl;


  //Reset for another block
  m_block_v = 0.;
  m_block_w = 0.;
  for(auto& block_gi : m_block_g) block_gi=0;
  m_block_index = 0;
  m_attempted = 0;
  m_accepted = 0;


  //Print the final values
  if(m_block_number == m_nblk){
    ofstream gfinal;
    gfinal.open("output.gave.0");
    if(gfinal.is_open()){
      gfinal<<"r_min,\t"<<"r_max,\t"<<"g(r),\t"<<"error"<<endl;
      for(int i=0; i<m_nbins; i++){
        gfinal<<i*m_histo->getBinLength()<<",\t"<<(i+1)*m_histo->getBinLength()<<",\t";
        gfinal<<m_ave_g[i]<<",\t"<<m_error_g[i]<<endl;
      }

    }else cerr<<"Unable to open output.gave.0"<<endl;
    gfinal.close();
  }

}
