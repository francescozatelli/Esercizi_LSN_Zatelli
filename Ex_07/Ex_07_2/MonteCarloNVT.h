#ifndef __MonteCarloNVT__
#define __MonteCarloNVT__

#include<cmath>
#include<string>
#include<cstdlib>
#include<vector>
#include<iostream>
#include<fstream>
#include "random.h"
#include "histogram.h"

using namespace std;

const double pi = 3.14159265358979323846;
 
class MonteCarloNVT {

private:
  Random *m_rnd;

  vector<double> m_g;
  vector<double> m_error_g;

  Histogram *m_histo;
  int m_nbins;

  vector<double> m_x,m_y,m_z;

  double m_temp, m_beta, m_delta, m_rho, m_rcut, m_vol, m_box;
  int m_npart, m_nblk, m_nstep, m_eq_nstep;
  bool m_restart, m_printall;
  int m_attempted, m_accepted;

  double m_vtail, m_wtail; //tail corrections fo rpotential energy and virial

  int m_block_index; //Index for each block
  int m_block_number; //Index for the total number of blocks


  //v->potential energy, w->virial (P=rho*kb*T+<w>/volume)
  double m_block_v, m_block_w;
  double m_ave_v, m_ave_w;
  double m_ave2_v, m_ave2_w;

  vector<double> m_block_g;
  vector<double> m_ave_g;
  vector<double> m_ave2_g;

  void FirstInitialization();
  void RestartInitialization();

  double Pbc(double r){ return r-m_box*rint(r/m_box);}

  //Returns the boltzmann weight for a given energy
  double BoltzmannWeight(double energy);
  double InteractionEnergy(double x, double y, double z, int ip);


public:
	MonteCarloNVT(Random *rnd);
	~MonteCarloNVT(){delete m_histo;};

  void Move();
  void Measure();
  void ConfFinal();
  void ConfXYZ(int nconf);
  void MeasureBlock();
  double getAcceptanceRate(){return m_accepted/double(m_attempted);};
  int getNumberOfBlocks(){return m_nblk;};
  int getNumberOfStepsPerBlock(){return m_nstep;};


};

#endif // __MonteCarloNVT__
