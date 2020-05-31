#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "metropolis.h"
#include "annealer.h"
#include "functions.h"
 
using namespace std;

int main (int argc, char *argv[]){

	Random *rnd=new Random();
	initRandom(*rnd); //Sposto dal main l'inizializzazione del generatore (seed ecc.)

	double mu0, sigma0;
	double samplingmove_width, sigmamove_width, mumove_width;
	double beta0, beta_step;
	int nsteps_beta, nsteps_musigma, nsteps_sampling, block_size;


	ifstream ReadInput;
	//Read input data for the simulation
	ReadInput.open("input.dat");
	if(ReadInput.is_open()){
		ReadInput>>mu0;
		ReadInput>>sigma0;
		ReadInput>>samplingmove_width;
		ReadInput>>sigmamove_width;
		ReadInput>>mumove_width;
		ReadInput>>beta0;
		ReadInput>>beta_step;
		ReadInput>>nsteps_beta;
		ReadInput>>nsteps_musigma;
		ReadInput>>nsteps_sampling;
		ReadInput>>block_size;
	}else cerr<<"ERROR: Unable to open input.dat"<<endl;
	ReadInput.close();

	cout<<"Simulated annealing for trial wave function"<<endl;
	cout<<"Initial mu: "<<mu0<<endl;
	cout<<"Initial sigma: "<<sigma0<<endl;
	cout<<"Initial beta: "<<beta0<<" - Initial T: "<<1./beta0<<endl;
	cout<<"Beta step: "<<beta_step<<endl;
	cout<<"Number of beta steps: "<<nsteps_beta<<endl;
	cout<<"Final beta: "<<nsteps_beta*beta_step+beta0<<endl;
	cout<<"Mu move width: "<<mumove_width<<endl;
	cout<<"Sigma move width: "<<sigmamove_width<<endl;
	cout<<"Number of moves (mu,sigma) for each beta step: "<<nsteps_musigma<<endl;
	cout<<"Smart sampling width: "<<samplingmove_width<<endl;
	cout<<"Number of throws for each energy calculation: "<<nsteps_sampling<<endl;
	cout<<"Number of blocks for energy estimation: "<<block_size<<endl;

	ofstream fileout;

	fileout.open("annealing.out");
	fileout<<"Beta,\t T,\t Mu,\t Sigma,\t Energy,\t EnergyError,\tSigmaMuAcceptance"<<endl;


	Annealer annealer(mu0, sigma0, samplingmove_width, nsteps_sampling, rnd);
	for(int i=0; i<=nsteps_beta; i++){
		cout<<"======================"<<endl;
		double beta = beta_step*i+beta0; 
		cout<<"Beta = "<<beta<<"; T = "<<1./beta<<endl;
		annealer.annealingStep(beta, nsteps_musigma, mumove_width, sigmamove_width, samplingmove_width);
		double energy=0;
		double energy_err=0;
		cout<<"Acceptance for a (sigma, mu) step: "<<annealer.getAcceptance()<<endl;
		tie(energy, energy_err)=annealer.dataBlockingEnergy(block_size);
		cout<<"Energy: "<<energy<<" +/- "<<energy_err<<endl;
		cout<<"Mu: "<<annealer.getMu()<<" - Sigma: "<<annealer.getSigma()<<endl; 
		fileout<<beta<<",\t"<<1./beta<<",\t"<<annealer.getMu()<<",\t"<<annealer.getSigma()<<",\t"<<energy<<",\t"<<energy_err<<",\t"<<annealer.getAcceptance()<<endl;
	}
	cout<<endl<<"========== Finished sampling =============="<<endl;
	rnd->SaveSeed();
	fileout.close();

	ofstream optout;
	optout.open("optimalmusigma.out");
	optout<<annealer.getMu()<<endl;
	optout<<annealer.getSigma()<<endl;
	optout.close();
	cout<<"Optimal mu and sigma printed in optimalmusigma.out"<<endl;
	
	
	return 0;
}
