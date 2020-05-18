#include "metropolis.h"


double F100 :: Eval(double x, double y, double z) const{
  double r=sqrt(x*x+y*y+z*z);
  return exp(-2*r);
}

double F210 :: Eval(double x, double y, double z) const{
  double r=sqrt(x*x+y*y+z*z);
  return z*z*exp(-r);
}

double Metropolis :: TransitionUnif(double yi, double delta){
  return m_rnd->Rannyu(yi-delta,yi+delta);
}

double Metropolis :: TransitionGauss(double yi, double sigma){
  return m_rnd->Gauss(yi,sigma);
}

double Metropolis :: AcceptanceProbability(double xnew, double ynew, double znew){
  double alpha = fmin(1, m_f->Eval(xnew,ynew,znew)/m_f->Eval(m_x, m_y, m_z));
  m_alpha=alpha;
  return alpha;
}


Metropolis :: Metropolis(double x, double y, double z, Function *f, Random *rnd){
  m_x=x;
  m_y=y;
  m_z=z;
  m_f=f;
  m_rnd=rnd;
}


void Metropolis :: StepUnif(double delta){
  double xnew=TransitionUnif(m_x, delta);
  double ynew=TransitionUnif(m_y, delta);
  double znew=TransitionUnif(m_z, delta);

  double alpha=AcceptanceProbability(xnew,ynew,znew);
  double r = m_rnd -> Rannyu(); //Uniform in [0,1) to test acceptance
  if(r < alpha){
    m_x=xnew;
    m_y=ynew;
    m_z=znew;
  } // else use the same point
}

void Metropolis :: StepGauss(double sigma){
  double xnew=TransitionGauss(m_x, sigma);
  double ynew=TransitionGauss(m_y, sigma);
  double znew=TransitionGauss(m_z, sigma);

  double alpha=AcceptanceProbability(xnew,ynew,znew);
  double r = m_rnd -> Rannyu(); //Uniform in [0,1) to test acceptance
  if(r < alpha){
    m_x=xnew;
    m_y=ynew;
    m_z=znew;
  } // else use the same point
}
