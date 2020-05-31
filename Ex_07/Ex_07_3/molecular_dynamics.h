#ifndef __MolecularDynamics__
#define __MolecularDynamics__

#include<cmath>
#include<string>
#include<cstdlib>
#include<vector>
#include<iostream>
#include<fstream>
#include "histogram.h"

using namespace std;

const double pi = 3.14159265358979323846;

class MolecularDynamics {

private:
  double m_temp, m_rho, m_rcut, m_dt, m_vol, m_box;
  int m_npart, m_nsteps;
  int m_block_index; //Index for each block
  int m_block_number; //Index for the total number of blocks
 
  double m_vtail, m_wtail; //tail corrections fo rpotential energy and virial


  bool m_restart, m_printall;

  vector<double> m_x,m_y,m_z;
  vector<double> m_xold,m_yold,m_zold;
  vector<double> m_vx,m_vy,m_vz;
  
  double m_block_epot, m_block_ekin, m_block_etot, m_block_temp, m_block_w;
  double m_ave_epot, m_ave_ekin, m_ave_etot, m_ave_temp, m_ave_w;
  double m_ave2_epot, m_ave2_ekin, m_ave2_etot, m_ave2_temp, m_ave2_w;

  Histogram *m_histo;
  int m_nbins;

  vector<double> m_block_g;
  vector<double> m_ave_g;
  vector<double> m_ave2_g;
  vector<double> m_g;
  vector<double> m_error_g;

  void FirstInitialization();
  void RestartInitialization();

  double Pbc(double r){ return r-m_box*rint(r/m_box);}
  double Force(int ipart, int idir);


public:
	MolecularDynamics();
	~MolecularDynamics(){};

  void Input();
  void Measure();
  void Move();
  void ConfXYZ(int nconf); //Print configuration in .xyz format
  void ConfFinal();
  void SaveOldConfig(); //Print configuration in order to be used again by the class to restart
  int getSteps(){return m_nsteps;}
  void MeasureBlock();

};

#endif // __MolecularDynamics__
