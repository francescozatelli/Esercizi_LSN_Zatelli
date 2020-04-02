#ifndef __BlackScholes__
#define __BlackScholes__

#include "random.h"

class BlackScholes {

private:
  double m_S0, m_K, m_r, m_sigma, m_T;
  Random *m_rnd;
  double N(double x){ return 0.5*(1+erf(x/sqrt(2)));}

public:
  //Asset price at t=0, strike price, risk-free interest rate, volatility, delivery time
	BlackScholes(double S0,double T, double K, double r, double sigma, Random *rnd);
	~BlackScholes(){};

  double callDirect(int n_repetitions);
  double putDirect(int n_repetitions);
  double callDiscretized(int n_steps, int n_repetitions);
  double putDiscretized(int n_steps, int n_repetitions);
  double callAnalytic();
  double putAnalytic();
};

#endif // __BlackScholes__
