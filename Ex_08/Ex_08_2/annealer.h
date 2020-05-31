#include "random.h"
#include "metropolis.h"
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

  double m_mu, m_sigma, m_energy; //values for the last accepted metropolis step
  int m_accepted, m_attempted;

  int m_metrosteps; //number of steps for each energy calculation
  double m_metrowidth; //width of the metropolis steps for energy calculation

  vector<double> m_energies; //energies of the last accepted metropolis step
  vector<double> m_energies_new; //energies for each metropolis step (it replaces m_energies if the step is accepted)

  vector<double> m_weights; //weights of the last accepted metropolis step
  vector<double> m_weights_new;  //weights for each metropolis step (it replaces m_weights  if the step is accepted)

  vector<double> m_points; //contains the sampled path
  vector<double> m_pevals; //contains the p(x) for each point in the sampled path
  ProbabilityDistribution m_psample; //p(x) used to sample the path

  double BoltzmannWeight(double beta, double energy);
  
  void samplePath(double mu, double sigma); //sample p_{mu,sigma}(x)
  double CalculateEnergy(double mu, double sigma); //compute the energy with the sampled path and reweighting
  double sampleAndCalculateEnergy(double mu, double sigma); //samples and computes the energy (no reweighting, slow)

public:
  Annealer(double mu0, double sigma0, double metrowidth, int metrosteps, Random *rnd);
  ~Annealer(){};

  void annealingStep(double beta, int nsteps, double mu_width, double sigma_width, double metrowidth);
  tuple<double, double> dataBlockingEnergy(int block_size); //returns average and error of last accepted metropolis step using m_energies and m_weights
  double getMu() const{return m_mu;};
  double getSigma() const{return m_sigma;};
  double getEnergyTemp() const{return m_energy;};
  double getAcceptance() const{return m_accepted/double(m_attempted);};
};

#endif //__Annealer__
