#ifndef __Ising1D__
#define __Ising1D__

#include<cmath>
#include<string>
#include<cstdlib>
#include<vector>
#include<iostream>
#include<fstream>
#include "random.h"

using namespace std;

class Ising1D {

private:
  Random *m_rnd;
  double m_temp, m_beta, m_J, m_h;
  int m_nspin, m_nblk, m_nstep, m_eq_nstep;
  bool m_metro, m_restart, m_printall;
  int m_attempted, m_accepted;

  double m_alpha;


  int m_block_index; //Index for each block
  int m_block_number; //Index for the total number of blocks


  vector<int> m_s; //spins array


  double m_block_U, m_block_M;
  double m_block_U2, m_block_M2; //necessary to calculate C and chi
  double m_ave_U, m_ave_C, m_ave_chi, m_ave_M;
  double m_ave2_U, m_ave2_C, m_ave2_chi, m_ave2_M;

  void FirstInitialization();
  void RestartInitialization();

  int Pbc(int i);
  //Returns the boltzmann weight for a given energy
  double BoltzmannWeight(double energy);
  //Returns the energy difference when the spin in position k flips
  double EnergyDifference(int k);
  //Returns the energy difference when the spin in position k flips if the k-th spin is sk
  double EnergyDifference(int k, int sk);
  double TotalEnergy();
  double TotalSpin();

  void MoveMetropolis(int k);
  double AcceptanceProbability(int k);

  void MoveGibbs(int k);
  double ConditionalProbability(int k, int sk);

public:
	Ising1D(Random *rnd);
	~Ising1D(){};

  void Move();
  void Measure();
  void ConfFinal();
  void MeasureBlock();
  double getAlpha(){return m_accepted/double(m_attempted);};
  int getNumberOfBlocks(){return m_nblk;};
  int getNumberOfStepsPerBlock(){return m_nstep;};


};

#endif // __Ising1D__
