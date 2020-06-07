#include "random.h"
#include "genetics.h"
#include <vector>
#include <tuple>
#include <cmath>
#include<iostream>
#include<fstream>

#ifndef __Annealer__
#define __Annealer__

using namespace std;

class Annealer{

private:
  Random *m_rnd;

  int m_accepted, m_attempted;
  Chromosome<int> m_chromosome;
  Chromosome<int> m_trial_chromosome;

  vector<pair<double,double>> m_coordinates;

  double BoltzmannWeight(double beta, double energy);
  double calculateLength(Chromosome<int> &chromosome);
  double distance(int town1, int town2) const;
  int randomGeneIndex(){ return m_rnd->UniformInteger(0,m_chromosome.size()-1);};

public:
  Annealer(int size, Random *rnd, vector<pair<double,double>> coordinates);
  Annealer(Chromosome<int> chromosome, Random *rnd, vector<pair<double,double>> coordinates);
  ~Annealer(){};

  void annealingStep(double beta, int nsteps);

  double getAcceptance() const{return m_accepted/double(m_attempted);};
  double getLength()const{ return 1./m_chromosome.getFitness();};
  double getFitness()const{ return m_chromosome.getFitness();};

  Chromosome<int> getChromosome()const{ return m_chromosome;};
};

#endif //__Annealer__

//Functions to generate random towns, not strictly related to annealer, placed here for convenience.

vector<pair<double,double>> townsOnCycle(int N, double R, Random *rnd);
vector<pair<double,double>> townsInSquare(int N, double edge, Random *rnd);
