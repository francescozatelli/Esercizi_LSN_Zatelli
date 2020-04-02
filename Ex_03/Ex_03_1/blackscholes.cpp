#include <iostream>
#include <cmath>
#include "random.h"
#include "blackscholes.h"

using namespace std;

BlackScholes :: BlackScholes(double S0, double T, double K, double r, double sigma, Random *rnd){
	m_S0=S0;
	m_T=T;
	m_K=K;
	m_r=r;
	m_sigma=sigma;
	m_rnd=rnd;
}

double BlackScholes :: callDirect(int n_repetitions){
	double avg_call=0;
	for(int i=0;i<n_repetitions;i++){
		double W=m_rnd->Gauss(0,m_T);
		double S=m_S0*exp((m_r-0.5*m_sigma*m_sigma)*m_T+m_sigma*W);
		double n_call=exp(-m_r*m_T)*fmax(0,S-m_K);
		avg_call=i/double(i+1)*avg_call+1/double(i+1)*n_call;
	}
	return avg_call;
}

double BlackScholes :: putDirect(int n_repetitions){
	double avg_call=0;
	for(int i=0;i<n_repetitions;i++){
		double W=m_rnd->Gauss(0,m_T);
		double S=m_S0*exp((m_r-0.5*m_sigma*m_sigma)*m_T+m_sigma*W);
		double n_call=exp(-m_r*m_T)*fmax(0,m_K-S);
		avg_call=i/double(i+1)*avg_call+1/double(i+1)*n_call;
	}
	return avg_call;
}



double BlackScholes :: callDiscretized(int n_steps, int n_repetitions){
	double avg_call=0;
	double time_step=m_T/double(n_steps);

	for(int i=0;i<n_repetitions;i++){
		double S=m_S0;

		for(int j=0;j<n_steps;j++){
			double Z=m_rnd->Gauss(0,1);
			S=S*exp((m_r-0.5*m_sigma*m_sigma)*time_step+m_sigma*Z*sqrt(time_step));
		}

		double n_call=exp(-m_r*m_T)*fmax(0,S-m_K);
		avg_call=i/double(i+1)*avg_call+1/double(i+1)*n_call;
	}
	return avg_call;
}


double BlackScholes :: putDiscretized(int n_steps, int n_repetitions){
	double avg_call=0;
	double time_step=m_T/double(n_steps);

	for(int i=0;i<n_repetitions;i++){
		double S=m_S0;

		for(int j=0;j<n_steps;j++){
			double Z=m_rnd->Gauss(0,1);
			S=S*exp((m_r-0.5*m_sigma*m_sigma)*time_step+m_sigma*Z*sqrt(time_step));
		}

		double n_call=exp(-m_r*m_T)*fmax(0,m_K-S);
		avg_call=i/double(i+1)*avg_call+1/double(i+1)*n_call;
	}
	return avg_call;
}

double BlackScholes :: callAnalytic(){
	double d1=1./(m_sigma*sqrt(m_T))*(log(m_S0/m_K)+m_r+0.5*m_sigma*m_sigma*m_T);
	double d2=d1-m_sigma*sqrt(m_T);
	double call=m_S0*N(d1)-m_K*exp(-m_r*m_T)*N(d2);
	return call;
}

double BlackScholes :: putAnalytic(){
	double d1=1./(m_sigma*sqrt(m_T))*(log(m_S0/m_K)+m_r+0.5*m_sigma*m_sigma*m_T);
	double d2=d1-m_sigma*sqrt(m_T);
	double put=m_S0*(N(d1)-1)-m_K*exp(-m_r*m_T)*(N(d2)-1);
	return put;
}
