#include "random.h"
#include <cmath>
#include<iostream>

#ifndef __Metropolis__
#define __Metropolis__

using namespace std;

class Function{
public:
    virtual double Eval(double x, double y, double z) const=0;
};

class F100 : public Function{
public:
  F100(){};
  ~F100(){};
  virtual double Eval(double x, double y, double z) const;
};

class F210 : public Function{
public:
  F210(){};
  ~F210(){};
  virtual double Eval(double x, double y, double z) const;
};

class Metropolis{

private:
  Random *m_rnd;
  double m_x, m_y, m_z;
  double m_alpha;
  Function *m_f;

  double TransitionUnif(double yi, double delta); //yi is one of the coordinates of the last point
  double TransitionGauss(double yi, double sigma);
  double AcceptanceProbability(double xnew, double ynew, double znew);

public:
  Metropolis(double x, double y, double z, Function *f, Random *rnd);
  ~Metropolis(){};

  void StepUnif(double delta);
  void StepGauss(double sigma);

  double getX() const{return m_x;}
  double getY() const{return m_y;}
  double getZ() const{return m_z;}

  double getAlpha() const{return m_alpha;}
};

#endif //__Metropolis__
