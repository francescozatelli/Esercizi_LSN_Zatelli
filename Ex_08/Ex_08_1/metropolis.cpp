#include "metropolis.h"

double LocalEnergy :: Eval(double x) const{
  double x2=x*x;
  double x4=x2*x2;
  double sigma2=m_sigma*m_sigma;
  double sigma4=sigma2*sigma2;
  double value=-(x2+m_mu*m_mu-sigma2-2*m_mu*x*tanh(m_mu*x/sigma2))/(2*sigma4)+x4-5./2.*x2;
  return value;
}

double ProbabilityDistribution :: Eval(double x) const{
  double value=exp(-(x-m_mu)*(x-m_mu)/(2*m_sigma*m_sigma))+exp(-(x+m_mu)*(x+m_mu)/(2*m_sigma*m_sigma));
  value=value*value;
  return value;
}

double LogarithmicDerivativeP :: Eval(double x) const{
  double value=-2*(x-m_mu*tanh(m_mu*x/(m_sigma*m_sigma)))/(m_sigma*m_sigma);
  return value;
}


double Metropolis :: GaussFunction(double x, double mu, double sigma) const{
  return exp(-pow(x-mu,2)/(2*pow(sigma,2)));
}
 

Metropolis :: Metropolis(double x, Function *p, Function *dlogp, Random *rnd){
  m_x=x;
  m_p=p;
  m_dlogp=dlogp;
  m_accepted = 0;
  m_attempted = 0;
  m_rnd=rnd;
}


void Metropolis :: StepUnif(double delta){
  double xnew=m_rnd->Rannyu(m_x-delta,m_x+delta);
  double q=m_p->Eval(xnew)/m_p->Eval(m_x);
  if(q >= m_rnd->Rannyu()){
    m_x=xnew;
    m_accepted++;
  } // else use the same point
  m_attempted++;
}

void Metropolis :: StepGauss(double sigma){
  double xnew=m_rnd->Gauss(m_x,sigma);
  double q=m_p->Eval(xnew)/m_p->Eval(m_x);
  if(q >= m_rnd->Rannyu()){
    m_x=xnew;
    m_accepted++;
  } // else use the same point
  m_attempted++;
}

void Metropolis :: StepSmart(double sigma){
  double dlogp_old=m_dlogp->Eval(m_x); //D ln p(xold)
  double xnew=m_rnd->Gauss(m_x+sigma*sigma/2.*dlogp_old, sigma); //Sample with drift
  double dlogp_new=m_dlogp->Eval(xnew); //D ln p(xnew)

  double txy=GaussFunction(xnew, m_x+sigma*sigma/2.*dlogp_old, sigma); //T(x|y)
  double tyx=GaussFunction(m_x, xnew+sigma*sigma/2.*dlogp_new, sigma); //T(y|x)
  
  double q = tyx/txy*m_p->Eval(xnew)/m_p->Eval(m_x); //acceptance=min(1,q)

  if(q >= m_rnd->Rannyu()){ //Test acceptance
    m_x=xnew;
    m_accepted++;
  } // else use the same point
  m_attempted++;
}