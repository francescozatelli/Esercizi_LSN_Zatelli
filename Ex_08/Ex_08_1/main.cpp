#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "metropolis.h"
#include "functions.h"

using namespace std;

int main (int argc, char *argv[]){

	int M = 1E6; //Number of throws
	int L = 10000; //Block size
	int N=M/L; //Number of blocks

	vector<double> eunif_throws(M); //single throw
	vector<double> esmart_throws(M);
	vector<double> eunif_blocks(N); //block values
	vector<double> esmart_blocks(N);
	vector<double> eunif_avg(N); //progressive average of blocks
	vector<double> esmart_avg(N);
	vector<double> eunif_avg2(N); //progressive average of blocks^2
	vector<double> esmart_avg2(N);
	vector<double> eunif_err(N); //progressive standard deviation
	vector<double> esmart_err(N);

	Random *rnd=new Random();
	initRandom(*rnd); //Sposto dal main l'inizializzazione del generatore (seed ecc.)
 
	//Parameters for trial wavefunction
	//double mu=0.792873; 
	//double sigma=0.615049;
	double mu=1; 
	double sigma=0.75;
	LocalEnergy *eloc= new LocalEnergy(mu, sigma); //local energy for trial wave function
	ProbabilityDistribution *p= new ProbabilityDistribution(mu, sigma); //probability distribution for trial wave function
	LogarithmicDerivativeP *dlogp= new LogarithmicDerivativeP(mu, sigma); //logarithmic derivative for trial wave function

	//Parameters for metropolis algorithm
	double delta=3.; //Uniform width
	double tau=.75; //Smart width
	double xstart=mu; //Starting point 

	Metropolis DumbMonteCarlo(xstart, p, dlogp, rnd);
	Metropolis SmartMonteCarlo(xstart, p, dlogp, rnd);

	cout<<"Energy expectation value for trial wavefunction with - Unifom transition probability vs smart Monte Carlo "<<endl;
	cout<<"Number of throws: "<<M<<endl;
	cout<<"Number of throws in each block: "<<L<<endl;
	cout<<"Number of blocks: "<<N<<endl;
	cout<<"Uniform distribution width: "<<delta<<endl;
	cout<<"Smart gaussian distribution width: "<<tau<<endl;
	cout<<"Starting point: "<<xstart<<endl;

	for(int i=0; i<M;i++){
		DumbMonteCarlo.StepUnif(delta);
		SmartMonteCarlo.StepSmart(tau);

		double xdumb=DumbMonteCarlo.getX();
		double xsmart=SmartMonteCarlo.getX();

		double eunif=eloc->Eval(xdumb);
		double esmart=eloc->Eval(xsmart);

		eunif_throws[i]=eunif;
		esmart_throws[i]=esmart;
	}
	cout<<endl<<"========== Finished sampling =============="<<endl;
	cout<<"Average acceptance rate for uniform Monte Carlo: "<<DumbMonteCarlo.getAlpha()<<endl;
	cout<<"Average accaptance rate for smart Monte Carlo: "<<SmartMonteCarlo.getAlpha()<<endl;
	rnd->SaveSeed();


	//Calcolo i valori per gli N blocchi (le N "misure")
	for(int i=0; i<N; i++){
		//Fa la media di una "fetta" di vettore che parte da i*L ed è lunga L
		//Per ogni blocco considero una sequenza di L=M/N numeri estratti casualmente
		eunif_blocks[i]=averageSlice(eunif_throws, i*L, L);
		esmart_blocks[i]=averageSlice(esmart_throws, i*L, L);
	}

	//printVector(eunif_blocks, "blocksunif.txt");
	//printVector(esmart_blocks, "blocks.txt");


	//Average and error
	for(int i=0; i<N; i++){
		eunif_avg[i]=averageSlice(eunif_blocks,0,i+1); //Calcola la media dei primi i+1 esperimenti
		eunif_avg2[i]=averageQuadSlice(eunif_blocks,0,i+1); //Calcola la media dei quadrati dei primi i+1 esperimenti
		eunif_err[i]=error(eunif_avg,eunif_avg2,i); //Calcola la varianza della media per i primi i+1 esperimenti

		esmart_avg[i]=averageSlice(esmart_blocks,0,i+1); //Calcola la media dei primi i+1 esperimenti
		esmart_avg2[i]=averageQuadSlice(esmart_blocks,0,i+1); //Calcola la media dei quadrati dei primi i+1 esperimenti
		esmart_err[i]=error(esmart_avg,esmart_avg2,i); //Calcola la varianza della media per i primi i+1 esperimenti
	}



	//Output dei risultati
	ofstream fileout;
	fileout.open("out_eunif.txt");
	if(fileout.is_open()){
		fileout<<"N,\tAverage,\tError"<<endl;
		for(int i=0; i<N; i++){
			fileout<<i+1<<",\t"<<eunif_avg[i]<<",\t"<<eunif_err[i]<<endl;
		}
	}else cerr<< "Unable to open output file"<<endl;

	fileout.close();

	fileout.open("out_esmart.txt");
	if(fileout.is_open()){
		fileout<<"N,\tAverage,\tError"<<endl;
		for(int i=0; i<N; i++){
			fileout<<i+1<<",\t"<<esmart_avg[i]<<",\t"<<esmart_err[i]<<endl;
		}
	}else cerr<< "Unable to open output file"<<endl;
	fileout.close();

	cout<<"Output files for estimated energy: out_eunif.txt, out_esmart.txt"<<endl;
	return 0;
}
