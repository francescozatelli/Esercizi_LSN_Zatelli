#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "functions.h"
#include "blackscholes.h"

using namespace std;

int main (int argc, char *argv[]){

	Random *rnd=new Random();
	initRandom(*rnd); //Sposto dal main l'inizializzazione del generatore (seed ecc.)

	double S0=100; //Asset price at t=0
	double T=1; //Delivery time
	int n_steps=100; //Number of steps for discretized GBM
	double K=100; //Strike price
	double r=0.1; //Risk-free interest rate
	double sigma=0.25; //Volatiliy
	BlackScholes black_scholes(S0,T,K,r,sigma,rnd);


	int M=10000; //Number of asset prices
	int N=100; //Number of blocks
	int L=M/N; //Number of "samples" per block

	vector<double> call_direct_blocks(N); //Singole "misure" per ognuno degli N blocchi/esperimenti
	vector<double> put_direct_blocks(N); //Singole "misure" per ognuno degli N blocchi/esperimenti
	vector<double> call_discretized_blocks(N); //Singole "misure" per ognuno degli N blocchi/esperimenti
	vector<double> put_discretized_blocks(N); //Singole "misure" per ognuno degli N blocchi/esperimenti

	vector<double> call_direct_avg(N); //Medie delle misure per ogni numero di esperimenti considerati da 0 a N
	vector<double> call_direct_avg2(N); //Medie dei quadrati delle misure per ogni numero di esperimenti considerati da 0 a N
	vector<double> call_direct_err(N); //Deviazione standard della media per ogni numero di esperimenti considerati da 0 a N
	vector<double> put_direct_avg(N); //Medie delle misure per ogni numero di esperimenti considerati da 0 a N
	vector<double> put_direct_avg2(N); //Medie dei quadrati delle misure per ogni numero di esperimenti considerati da 0 a N
	vector<double> put_direct_err(N); //Deviazione standard della media per ogni numero di esperimenti considerati da 0 a N
	vector<double> call_discretized_avg(N); //Medie delle misure per ogni numero di esperimenti considerati da 0 a N
	vector<double> call_discretized_avg2(N); //Medie dei quadrati delle misure per ogni numero di esperimenti considerati da 0 a N
	vector<double> call_discretized_err(N); //Deviazione standard della media per ogni numero di esperimenti considerati da 0 a N
	vector<double> put_discretized_avg(N); //Medie delle misure per ogni numero di esperimenti considerati da 0 a N
	vector<double> put_discretized_avg2(N); //Medie dei quadrati delle misure per ogni numero di esperimenti considerati da 0 a N
	vector<double> put_discretized_err(N); //Deviazione standard della media per ogni numero di esperimenti considerati da 0 a N



	//Calcolo i valori di ciascuno degli N blocchi
	for(int i=0;i<N; i++){
		//Il valore di ogni blocco viene ottenuto facendo la media di L valori (media fatta all'interno delle funzioni)
		call_direct_blocks[i]=black_scholes.callDirect(L);
		put_direct_blocks[i]=black_scholes.putDirect(L);
		call_discretized_blocks[i]=black_scholes.callDiscretized(n_steps,L);
		put_discretized_blocks[i]=black_scholes.putDiscretized(n_steps,L);
	}

	rnd->SaveSeed();

	for(int i=0; i<N; i++){
		call_direct_avg[i]=averageSlice(call_direct_blocks,0,i+1); //Calcola la media dei primi i+1 esperimenti
		call_direct_avg2[i]=averageQuadSlice(call_direct_blocks,0,i+1); //Calcola la media dei quadrati dei primi i+1 esperimenti
		call_direct_err[i]=error(call_direct_avg,call_direct_avg2,i); //Calcola la varianza della media per i primi i+1 esperimenti

		put_direct_avg[i]=averageSlice(put_direct_blocks,0,i+1); //Calcola la media dei primi i+1 esperimenti
		put_direct_avg2[i]=averageQuadSlice(put_direct_blocks,0,i+1); //Calcola la media dei quadrati dei primi i+1 esperimenti
		put_direct_err[i]=error(put_direct_avg,put_direct_avg2,i); //Calcola la varianza della media per i primi i+1 esperimenti

		call_discretized_avg[i]=averageSlice(call_discretized_blocks,0,i+1); //Calcola la media dei primi i+1 esperimenti
		call_discretized_avg2[i]=averageQuadSlice(call_discretized_blocks,0,i+1); //Calcola la media dei quadrati dei primi i+1 esperimenti
		call_discretized_err[i]=error(call_discretized_avg,call_discretized_avg2,i); //Calcola la varianza della media per i primi i+1 esperimenti

		put_discretized_avg[i]=averageSlice(put_discretized_blocks,0,i+1); //Calcola la media dei primi i+1 esperimenti
		put_discretized_avg2[i]=averageQuadSlice(put_discretized_blocks,0,i+1); //Calcola la media dei quadrati dei primi i+1 esperimenti
		put_discretized_err[i]=error(put_discretized_avg,put_discretized_avg2,i); //Calcola la varianza della media per i primi i+1 esperimenti

	}

	double call_analytic=black_scholes.callAnalytic();
	double put_analytic=black_scholes.putAnalytic();

	//Output dei risultati
	ofstream fileout;

	fileout.open("out_call_direct.txt");
	if(fileout.is_open()){
		fileout<<"N,\tAverage,\tError"<<endl;
		for(int i=0; i<N; i++){
			fileout<<i+1<<",\t"<<call_direct_avg[i]<<",\t"<<call_direct_err[i]<<endl;
		}
	}else cerr<< "Unable to open output file"<<endl;
	fileout.close();

	fileout.open("out_put_direct.txt");
	if(fileout.is_open()){
		fileout<<"N,\tAverage,\tError"<<endl;
		for(int i=0; i<N; i++){
			fileout<<i+1<<",\t"<<put_direct_avg[i]<<",\t"<<put_direct_err[i]<<endl;
		}
	}else cerr<< "Unable to open output file"<<endl;
	fileout.close();

	fileout.open("out_call_discretized.txt");
	if(fileout.is_open()){
		fileout<<"N,\tAverage,\tError"<<endl;
		for(int i=0; i<N; i++){
			fileout<<i+1<<",\t"<<call_discretized_avg[i]<<",\t"<<call_discretized_err[i]<<endl;
		}
	}else cerr<< "Unable to open output file"<<endl;
	fileout.close();

	fileout.open("out_put_discretized.txt");
	if(fileout.is_open()){
		fileout<<"N,\tAverage,\tError"<<endl;
		for(int i=0; i<N; i++){
			fileout<<i+1<<",\t"<<put_discretized_avg[i]<<",\t"<<put_discretized_err[i]<<endl;
		}
	}else cerr<< "Unable to open output file"<<endl;
	fileout.close();

	fileout.open("out_call_put_analytic.txt");
	if(fileout.is_open()){
		fileout<<"Call - Analytic:"<<endl<<call_analytic<<endl;
		fileout<<"Put - Analytic:"<<endl<<put_analytic<<endl;
	}else cerr<< "Unable to open output file"<<endl;
	fileout.close();

	return 0;
}
