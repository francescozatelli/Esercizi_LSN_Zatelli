#include "random.h"
#include <cmath>
#include<iostream>

#ifndef __Metropolis__
#define __Metropolis__

using namespace std;


class Function{
public:
    virtual double Eval(double x) const=0;
};

class LocalEnergy : public Function{ //Function to evaluate (local energy), it is not needed for the sampling (here just for convenience)
private:
  double m_sigma;
  double m_mu;
public:
  LocalEnergy(){m_sigma=1; m_mu=1;};
  LocalEnergy(double mu, double sigma){m_mu=mu; m_sigma=sigma;};
  ~LocalEnergy(){};
  virtual double Eval(double x) const;
};

class ProbabilityDistribution : public Function{ //Probability distribution to sample
private:
  double m_sigma;
  double m_mu;
public:
  ProbabilityDistribution(){m_sigma=1; m_mu=1;};
  ProbabilityDistribution(double mu, double sigma){m_mu=mu; m_sigma=sigma;};
  ~ProbabilityDistribution(){};
  virtual double Eval(double x) const;
};

class LogarithmicDerivativeP : public Function{ //Logarithmic derivative of the distribution to sample
private:
  double m_sigma;
  double m_mu;
public:  
  LogarithmicDerivativeP(){m_sigma=1; m_mu=1;};
  LogarithmicDerivativeP(double mu, double sigma){m_mu=mu; m_sigma=sigma;};
  ~LogarithmicDerivativeP(){};
  virtual double Eval(double x) const;
};

//NB: The sigma above is unrelated to the sigma below
//The first one refers to the trial function, the other to the size of a metropolis step

class Metropolis{

private:
  Random *m_rnd;
  Function *m_p; //probability distribution to sample
  Function *m_dlogp; //logarithmic derivative of the probability distribution to sample

  double m_x;
  int m_attempted;
  int m_accepted;
  
  double GaussFunction(double x, double mu, double sigma) const; // Used for T(X|Y) and T(Y|X), not normalized

public:
  Metropolis(double x, Function *p, Function *dlogp, Random *rnd);
  ~Metropolis(){};

  void StepUnif(double delta);
  void StepGauss(double sigma);
  void StepSmart(double sigma);

  double getX() const{return m_x;}
  double getAlpha() const{return m_accepted/double(m_attempted);}
  void reset(){m_attempted=0; m_accepted=0;};

};

#endif //__Metropolis__
