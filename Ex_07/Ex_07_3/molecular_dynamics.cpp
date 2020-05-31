#include "molecular_dynamics.h"

MolecularDynamics :: MolecularDynamics(){
  m_block_epot = 0.;
  m_block_ekin = 0.;
  m_block_etot = 0.;
  m_block_temp = 0.;
  m_block_w = 0.;

  m_ave_epot = 0.;
  m_ave_ekin = 0.;
  m_ave_etot = 0.;
  m_ave_temp = 0.;
  m_ave_w = 0.;
  m_ave2_epot = 0.;
  m_ave2_ekin = 0.;
  m_ave2_etot = 0.;
  m_ave2_temp = 0.;
  m_ave2_w = 0.;

  m_block_index = 0;
  m_block_number = 0;
 
  ifstream ReadInput;
  //Read input data for the simulation
  ReadInput.open("input.dat");
  if(ReadInput.is_open()){
    ReadInput >> m_temp; //Temperature
    ReadInput >> m_npart; //Number of particles
    ReadInput >> m_rho; //Density
    ReadInput >> m_rcut; //Cutoff radius for potential evaluation
    ReadInput >> m_dt; //Time step
    ReadInput >> m_nsteps; //Number of steps for simulation
    ReadInput >> m_restart; //Restart from old config
    ReadInput >> m_printall; //Restart from old config
    ReadInput >> m_nbins; //Number of bins for g(r) histogram
  }else cerr<<"Unable to open input.dat"<<endl;

  cout<<"Input data from input.dat (LJ units)"<<endl;
  cout<<"Temperature: "<<m_temp<<endl;
  cout<<"Number of particles: "<<m_npart<<endl;
  cout<<"Density: "<<m_rho<<endl;

  m_vol = (double)m_npart/m_rho; // Volume
  m_box = pow(m_vol,1./3.); // Cube box, linear dimension
  cout<<"Volume of the simulation box: "<<m_vol<<endl;
  cout<<"Edge of the simulation box: "<<m_box<<endl;

  cout<<"Cutoff radius: "<<m_rcut<<endl;
  cout<<"Time step: "<<m_dt<<endl;
  cout<<"Number of steps: "<<m_nsteps<<endl;
  cout<<"Starting from old configuration: "<<m_restart<<endl;


  ReadInput.close();


  m_x = vector<double>(m_npart);
  m_y = vector<double>(m_npart);
  m_z = vector<double>(m_npart);
  m_xold = vector<double>(m_npart);
  m_yold = vector<double>(m_npart);
  m_zold = vector<double>(m_npart);
  m_vx = vector<double>(m_npart);
  m_vy = vector<double>(m_npart);
  m_vz = vector<double>(m_npart);


  m_block_g = vector<double>(m_nbins);
  m_ave_g = vector<double>(m_nbins);
  m_ave2_g = vector<double>(m_nbins);
  m_g = vector<double>(m_nbins);
  m_error_g = vector<double>(m_nbins);
  m_histo = new Histogram(m_nbins,0,m_box/2.);
  cout<<"Number of bins: "<<m_nbins<<endl;
  cout<<"Bin length: "<<m_histo->getBinLength()<<endl;


  if(m_restart){
    RestartInitialization();
  }else{
    FirstInitialization();
  }

  m_vtail = (8.0*pi*m_rho)/(9.0*pow(m_rcut,9)) - (8.0*pi*m_rho)/(3.0*pow(m_rcut,3));
  m_wtail = (32.0*pi*m_rho)/(9.0*pow(m_rcut,9)) - (16.0*pi*m_rho)/(3.0*pow(m_rcut,3));
  cout <<endl<< "Tail correction for the potential energy: " << m_vtail << endl;
  cout << "Tail correction for the virial: " << m_wtail << endl << endl;
}

void MolecularDynamics :: RestartInitialization(){
  //Read configuration r(t)
  ifstream ReadConf;
  cout<<"Reading initial configuration r(t) from old.0"<<endl;
  ReadConf.open("old.0");
  if(ReadConf.is_open()){
    for(int i=0; i<m_npart; i++){
      ReadConf >> m_x[i] >> m_y[i] >> m_z[i];
      m_x[i] = m_x[i] * m_box;
      m_y[i] = m_y[i] * m_box;
      m_z[i] = m_z[i] * m_box;
    }
  }else cerr<<"Unable to open old.0"<<endl;
  ReadConf.close();

  //Read configuration r(t-dt)
  cout<<"Reading initial configuration r(t-dt) from old.final"<<endl;
  ReadConf.open("old.final");
  if(ReadConf.is_open()){
    for(int i=0; i<m_npart; i++){
      ReadConf >> m_xold[i] >> m_yold[i] >> m_zold[i];
      m_xold[i] = m_xold[i] * m_box;
      m_yold[i] = m_yold[i] * m_box;
      m_zold[i] = m_zold[i] * m_box;
    }
  }else cerr<<"Unable to open old.final"<<endl;
  ReadConf.close();

  Move(); //Compute r(t+dt)

  //Velocities at t+dt/2
  for(int i=0; i<m_npart; i++){
    m_vx[i]=Pbc(m_x[i]-m_xold[i])/(m_dt);
    m_vy[i]=Pbc(m_y[i]-m_yold[i])/(m_dt);
    m_vz[i]=Pbc(m_z[i]-m_zold[i])/(m_dt);
  }


  //Compute mean of square velocities
  double avg_v2=0;
  for(int i=0; i<m_npart; i++){
    //Square velocity of particle i
    double v2 = pow(m_vx[i],2) + pow(m_vy[i],2) + pow(m_vz[i],2);
    //Compute mean square velocity on the fly
    avg_v2 = i/double(i+1)*avg_v2+1./double(i+1)*v2;
  }

  //Rescale the velocities in order to obtain the input temperature
  // fs = velocity scale factor
  double fs = sqrt(3 * m_temp / avg_v2);
  for(int i=0; i<m_npart; i++){
    //Rescale velocities
    m_vx[i] = fs * m_vx[i];
    m_vy[i] = fs * m_vy[i];
    m_vz[i] = fs * m_vz[i];

    //Correct the r(t-dt) configuration
    m_xold[i] = Pbc(m_x[i]-m_vx[i]*m_dt);
    m_yold[i] = Pbc(m_y[i]-m_vy[i]*m_dt);
    m_zold[i] = Pbc(m_z[i]-m_vz[i]*m_dt);
  }

}

void MolecularDynamics :: FirstInitialization(){
  //Read the initial configuration of the particles
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

  srand(1);
  vector<double> avg_v(3); //Already initialized with zeros
  for(int i=0; i<m_npart; i++){
    //Extract velocity components uniformly between -0.5 and 0.5
    m_vx[i] = rand()/double(RAND_MAX) - 0.5;
    m_vy[i] = rand()/double(RAND_MAX) - 0.5;
    m_vz[i] = rand()/double(RAND_MAX) - 0.5;

    avg_v[0] = i/double(i+1)*avg_v[0]+1./double(i+1)*m_vx[i];
    avg_v[1] = i/double(i+1)*avg_v[1]+1./double(i+1)*m_vy[i];
    avg_v[2] = i/double(i+1)*avg_v[2]+1./double(i+1)*m_vz[i];
  }

  double avg_v2=0;
  for(int i=0; i<m_npart; i++){
    //Translate the velocities in order to have zero centre of mass velocity
    m_vx[i] = m_vx[i]-avg_v[0];
    m_vy[i] = m_vy[i]-avg_v[1];
    m_vz[i] = m_vz[i]-avg_v[2];

    //Square velocity of particle i
    double v2 = pow(m_vx[i],2) + pow(m_vy[i],2) + pow(m_vz[i],2);
    //Compute mean square velocity on the fly
    avg_v2 = i/double(i+1)*avg_v2+1./double(i+1)*v2;
  }

  //Rescale the velocities in order to obtain the input temperature
  // fs = velocity scale factor
  double fs = sqrt(3 * m_temp / avg_v2);
  for(int i=0; i<m_npart; i++){
    m_vx[i] = fs * m_vx[i]; 
    m_vy[i] = fs * m_vy[i];
    m_vz[i] = fs * m_vz[i];

    m_xold[i] = Pbc(m_x[i]-m_vx[i]*m_dt);
    m_yold[i] = Pbc(m_y[i]-m_vy[i]*m_dt);
    m_zold[i] = Pbc(m_z[i]-m_vz[i]*m_dt);

  }
}


double MolecularDynamics :: Force(int ipart, int idir){
  double force_component=0.; //Force component in direction idir
  double dr; //Distance between two particles
  //vector<double> dvec(3); //Vector between two particles
  double dvec[3];
  for(int i=0; i<m_npart; i++){
    if(i != ipart){
      //Computing distance between two particles using pbc
      dvec[0] = Pbc(m_x[ipart]-m_x[i]);
      dvec[1] = Pbc(m_y[ipart]-m_y[i]);
      dvec[2] = Pbc(m_z[ipart]-m_z[i]);
      dr = dvec[0]*dvec[0]+dvec[1]*dvec[1]+dvec[2]*dvec[2];
      dr=sqrt(dr);

      if(dr < m_rcut){ //Check if particle i is inside the cutoff radius of particle ipart
        //Computing the contribute of particle i to the force acting on particle ipart in direction idir
        force_component+= dvec[idir]*(48./pow(dr,14) - 24./pow(dr,8)); //- gradient of potential for particle ipart
      }
    }
  }
  return force_component;
}


void MolecularDynamics :: Move(){
  //Move particles with Verlet algorithm
  double xnew, ynew, znew;
  vector<double> fx(m_npart);
  vector<double> fy(m_npart);
  vector<double> fz(m_npart);
  
  for(int i=0;i<m_npart; i++){
    fx[i]=Force(i,0);
    fy[i]=Force(i,1);
    fz[i]=Force(i,2);
  }

  for(int i=0; i<m_npart; i++){    //Verlet algorithm

    //Positions at t+dt
    xnew = Pbc(2.*m_x[i]-m_xold[i]+fx[i]*m_dt*m_dt);
    ynew = Pbc(2.*m_y[i]-m_yold[i]+fy[i]*m_dt*m_dt);
    znew = Pbc(2.*m_z[i]-m_zold[i]+fz[i]*m_dt*m_dt);

    //Velocities at t
    m_vx[i]=Pbc(xnew-m_xold[i])/(2.*m_dt);
    m_vy[i]=Pbc(ynew-m_yold[i])/(2.*m_dt);
    m_vz[i]=Pbc(znew-m_zold[i])/(2.*m_dt);

    //Update positions
    m_xold[i]=m_x[i];
    m_yold[i]=m_y[i];
    m_zold[i]=m_z[i];
    m_x[i]=xnew;
    m_y[i]=ynew;
    m_z[i]=znew;
  }
}

void MolecularDynamics :: Measure(){
  double dx,dy,dz,dr;
  double v = 0.; //Total potential energy
  double w = 0.; //Virial
  m_histo->reset(); //set all the bins to zero

  //Compute potential energy, cycle over pairs of particles
  for(int i=0; i<m_npart-1; i++){
    for(int j=i+1; j<m_npart; j++){
      dx = Pbc(m_x[i] - m_x[j]);
      dy = Pbc(m_y[i] - m_y[j]);
      dz = Pbc(m_z[i] - m_z[j]);
      dr = sqrt(dx*dx + dy*dy +dz*dz);
      
      m_histo->fill(dr);
      m_histo->fill(dr); //increase bin by 2
      
      if(dr<m_rcut){
        double dr6=pow(dr,6);
        double dr12=dr6*dr6;
        double vij=1./dr12-1./dr6;
        double wij=1./dr12-0.5/dr6;
        v+=vij;
        w+=wij;
      }
    }
    
  }
  v=4.*v+m_vtail;
  w=48.*w/3.+m_wtail*m_npart;

  double t=0.; //Total kinetic energy
  for(int i=0; i<m_npart; i++){
    t+= 0.5 * (m_vx[i]*m_vx[i] + m_vy[i]*m_vy[i] + m_vz[i]*m_vz[i]);
  }

  double epot_per_part = v/double(m_npart);
  double ekin_per_part = t/double(m_npart);
  double temp = 2./3. * ekin_per_part; //Equipartition theorem
  double etot_per_part = epot_per_part+ekin_per_part;

  m_block_epot=m_block_index/double(m_block_index+1)*m_block_epot+1./double(m_block_index+1)*epot_per_part;
  m_block_ekin=m_block_index/double(m_block_index+1)*m_block_ekin+1./double(m_block_index+1)*ekin_per_part;
  m_block_etot=m_block_index/double(m_block_index+1)*m_block_etot+1./double(m_block_index+1)*etot_per_part;
  m_block_temp=m_block_index/double(m_block_index+1)*m_block_temp+1./double(m_block_index+1)*temp;
  m_block_w=m_block_index/double(m_block_index+1)*m_block_w+1./double(m_block_index+1)*w;
  

  double bin_length = m_histo->getBinLength();
  for(int i=0; i<m_histo->getNumberOfBins(); i++){
    double norm = m_rho*m_npart*4./3.*pi*(pow((i+1)*bin_length,3)-pow(i*bin_length,3));
    m_g[i] = m_histo->getBin(i)/norm;
  }
  
  for(int i=0; i<m_nbins; i++){
    m_block_g[i]=m_block_index/double(m_block_index+1)*m_block_g[i]+1./double(m_block_index+1)*m_g[i];
  }

  m_block_index++;
  if(m_printall){
    ofstream Epot, Ekin, Etot, Temp, Pres;

    Epot.open("output_epot.dat",ios::app); //ios::app appends at the end of the file
    Ekin.open("output_ekin.dat",ios::app);
    Temp.open("output_temp.dat",ios::app);
    Etot.open("output_etot.dat",ios::app);
    Pres.open("output_pres.dat",ios::app);

    if(Epot.is_open()){
      Epot<<epot_per_part<<endl;
    }else cerr<<"Unable to open output_epot.dat"<<endl;
    if(Ekin.is_open()){
      Ekin<<ekin_per_part<<endl;
    }else cerr<<"Unable to open output_ekin.dat"<<endl;
    if(Temp.is_open()){
      Temp<<temp<<endl;
    }else cerr<<"Unable to open output_temp.dat"<<endl;
    if(Etot.is_open()){
      Etot<<etot_per_part<<endl;
    }else cerr<<"Unable to open output_etot.dat"<<endl;
    if(Pres.is_open()){
      double P=m_rho*m_temp+w/m_vol;
      Pres<<P<<endl;
    }else cerr<<"Unable to open output_pres.dat"<<endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pres.close();
  }
}


void MolecularDynamics :: ConfXYZ(int nconf){
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


void MolecularDynamics :: ConfFinal(){
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


void MolecularDynamics :: SaveOldConfig(){
  ofstream WriteConf;

  cout<<"Printing final configuration r(t-dt) to file old.final"<<endl;
  WriteConf.open("old.final");
  if(WriteConf.is_open()){
    for(int i=0; i<m_npart; i++){
      WriteConf<<m_xold[i]/m_box<<" "<<m_yold[i]/m_box<<" "<<m_zold[i]/m_box<<endl;
    }
  }else cerr<<"Unable to open old.final"<<endl;
  WriteConf.close();

  cout<<"Printing final configuration r(t) to file old.0"<<endl;
  WriteConf.open("old.0");
  if(WriteConf.is_open()){
    for(int i=0; i<m_npart; i++){
      WriteConf<<m_x[i]/m_box<<" "<<m_y[i]/m_box<<" "<<m_z[i]/m_box<<endl;
    }
  }else cerr<<"Unable to open old.0"<<endl;
  WriteConf.close();

}

void MolecularDynamics :: MeasureBlock(){
  double error_epot, error_ekin, error_temp, error_etot, error_w;
  m_ave_epot = m_block_number/double(m_block_number+1)*m_ave_epot+1./double(m_block_number+1)*m_block_epot;
  m_ave2_epot = m_block_number/double(m_block_number+1)*m_ave2_epot+1./double(m_block_number+1)*m_block_epot*m_block_epot;
  m_ave_ekin = m_block_number/double(m_block_number+1)*m_ave_ekin+1./double(m_block_number+1)*m_block_ekin;
  m_ave2_ekin = m_block_number/double(m_block_number+1)*m_ave2_ekin+1./double(m_block_number+1)*m_block_ekin*m_block_ekin;
  m_ave_temp = m_block_number/double(m_block_number+1)*m_ave_temp+1./double(m_block_number+1)*m_block_temp;
  m_ave2_temp = m_block_number/double(m_block_number+1)*m_ave2_temp+1./double(m_block_number+1)*m_block_temp*m_block_temp;
  m_ave_etot = m_block_number/double(m_block_number+1)*m_ave_etot+1./double(m_block_number+1)*m_block_etot;
  m_ave2_etot = m_block_number/double(m_block_number+1)*m_ave2_etot+1./double(m_block_number+1)*m_block_etot*m_block_etot;
  m_ave_w = m_block_number/double(m_block_number+1)*m_ave_w+1./double(m_block_number+1)*m_block_w;
  m_ave2_w = m_block_number/double(m_block_number+1)*m_ave2_w+1./double(m_block_number+1)*m_block_w*m_block_w;
  for(int i=0; i<m_nbins; i++){
    m_ave_g[i]=m_block_number/double(m_block_number+1)*m_ave_g[i]+1./double(m_block_number+1)*m_block_g[i];
    m_ave2_g[i]=m_block_number/double(m_block_number+1)*m_ave2_g[i]+1./double(m_block_number+1)*m_block_g[i]*m_block_g[i];
  }

  if(m_block_number == 0){
    error_epot=0;
    error_ekin=0;
    error_temp=0;
    error_etot=0;
    error_w=0;
    for(auto& error_gi : m_error_g) error_gi=0;

  }else{
    error_epot=sqrt((m_ave2_epot-m_ave_epot*m_ave_epot)/m_block_number);
    error_ekin=sqrt((m_ave2_ekin-m_ave_ekin*m_ave_ekin)/m_block_number);
    error_temp=sqrt((m_ave2_temp-m_ave_temp*m_ave_temp)/m_block_number);
    error_etot=sqrt((m_ave2_etot-m_ave_etot*m_ave_etot)/m_block_number);
    error_w=sqrt((m_ave2_w-m_ave_w*m_ave_w)/m_block_number);
    for(int i=0; i<m_nbins; i++){
      m_error_g[i]=sqrt((m_ave2_g[i]-m_ave_g[i]*m_ave_g[i])/m_block_number);
    }
  }

  m_block_number++;
 
  double p=m_rho*m_temp+m_ave_w/m_vol;
  double error_p=error_w/m_vol;

  ofstream Epot, Ekin, Etot, Temp, Pres, gofr, gave;

  Epot.open("ave_epot.out",ios::app); //ios::app appends at the end of the file
  Ekin.open("ave_ekin.out",ios::app);
  Temp.open("ave_temp.out",ios::app);
  Etot.open("ave_etot.out",ios::app);
  Pres.open("ave_pres.out",ios::app);
  gofr.open("output.gofr.0",ios::app);
  gave.open("output.gave.0"); //overwrite the last one every time in order to have only the "best" values

  if(Epot.is_open()){
    Epot<<m_block_number<<",\t"<<m_ave_epot<<",\t"<<error_epot<<endl;
  }else cerr<<"Unable to open ave_epot.out"<<endl;
  if(Ekin.is_open()){
    Ekin<<m_block_number<<",\t"<<m_ave_ekin<<",\t"<<error_ekin<<endl;
  }else cerr<<"Unable to open ave_ekin.out"<<endl;
  if(Temp.is_open()){
    Temp<<m_block_number<<",\t"<<m_ave_temp<<",\t"<<error_temp<<endl;
  }else cerr<<"Unable to open ave_temp.out"<<endl;
  if(Etot.is_open()){
    Etot.precision(10); //Double value is truncated too early otherwise
    Etot<<m_block_number<<",\t"<<m_ave_etot<<",\t"<<error_etot<<endl;
  }else cerr<<"Unable to open ave_etot.out"<<endl;
  if(Pres.is_open()){
    //Pres.precision(10); //Double value is truncated too early otherwise
    Pres<<m_block_number<<",\t"<<p<<",\t"<<error_p<<endl;
  }else cerr<<"Unable to open ave_pres.out"<<endl;

  //Stampo il valore di g(r) blocco per blocco (non facendo la media)
  if(gofr.is_open()){
    double bin_length=m_histo->getBinLength();
    gofr<<"===== BLOCK NUMBER "<<m_block_number<<" ====="<<endl;
    for(int i=0; i<m_nbins; i++){
      gofr<<m_block_number<<",\t"<<i*bin_length<<",\t"<<(i+1)*bin_length<<",\t"<<m_block_g[i]<<",\t"<<endl;
    }
  }else cerr<<"Unable to open output.gofr.0"<<endl;

  if(gave.is_open()){
      gave<<"r_min,\t"<<"r_max,\t"<<"g(r),\t"<<"error"<<endl;
      for(int i=0; i<m_nbins; i++){
        gave<<i*m_histo->getBinLength()<<",\t"<<(i+1)*m_histo->getBinLength()<<",\t";
        gave<<m_ave_g[i]<<",\t"<<m_error_g[i]<<endl;
      }

    }else cerr<<"Unable to open output.gave.0"<<endl;

  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  Pres.close();
  gofr.close();
  gave.close();

  m_block_epot = 0.;
  m_block_ekin = 0.;
  m_block_etot = 0.;
  m_block_temp = 0.;
  m_block_w=0.;
  m_block_index = 0;
  for(auto& block_gi : m_block_g) block_gi=0;

}
