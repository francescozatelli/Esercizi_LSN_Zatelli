#ifndef __MolecularDynamics__
#define __MolecularDynamics__

#include<cmath>
#include<string>
#include<cstdlib>
#include<vector>
#include<iostream>
#include<fstream>

using namespace std;

class MolecularDynamics {

private:
  double m_temp, m_rho, m_rcut, m_dt, m_vol, m_box;
  int m_npart, m_nsteps;
  bool m_restart, m_printall;

  vector<double> m_x,m_y,m_z;
  vector<double> m_xold,m_yold,m_zold;
  vector<double> m_vx,m_vy,m_vz;

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
  void ConfFinal(); //Print configuration in order to be used again by the class to restart
  void SaveOldConfig();

  int getSteps(){return m_nsteps;}
};

#endif // __MolecularDynamics__
