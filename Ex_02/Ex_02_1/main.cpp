#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "functions.h"

using namespace std;
 

double f_unif(double x){
	return M_PI/2*cos(M_PI/2*x);
}

double f_imp(double x){
	return M_PI*cos(M_PI/2*x)/(4.*(1-x));
}

int main (int argc, char *argv[]){

	Random rnd;
	initRandom(rnd); //Sposto dal main l'inizializzazione del generatore (seed ecc.)
	
	int M=10000; //Numero di estrazioni
	int N=100; //Numero di blocchi
	int L=M/N; //Numero di numeri casuali da considerare in ogni blocco

	vector<double> uniform_sampling(M);
	vector<double> importance_sampling(M);

	vector<double> unif_values_block(N); //Singole "misure" per ognuno degli N blocchi/esperimenti
	vector<double> imp_values_block(N); //Singole "misure" per ognuno degli N blocchi/esperimenti

	vector<double> unif_avg(N); //Medie delle misure per ogni numero di esperimenti considerati da 0 a N
	vector<double> unif_avg2(N); //Medie dei quadrati delle misure per ogni numero di esperimenti considerati da 0 a N
	vector<double> unif_err(N); //Deviazione standard della media per ogni numero di esperimenti considerati da 0 a N
	
	vector<double> imp_avg(N); //Medie delle misure per ogni numero di esperimenti considerati da 0 a N
	vector<double> imp_avg2(N); //Medie dei quadrati delle misure per ogni numero di esperimenti considerati da 0 a N
	vector<double> imp_err(N); //Deviazione standard della media per ogni numero di esperimenti considerati da 0 a N
	
	//Estraggo M numeri casuali estratti uniformemente in [0,1)
	for(int i=0; i<M; i++){
		double extracted_uniform=rnd.Rannyu();
		uniform_sampling[i]=f_unif(extracted_uniform);
		importance_sampling[i]=f_imp(1-sqrt(1-extracted_uniform));
	}
	rnd.SaveSeed();


	//Calcolo i valori per gli N blocchi (le N "misure")
	for(int i=0; i<N; i++){
		//Fa la media di una "fetta" di vettore che parte da i*L ed è lunga L
		//Per ogni blocco considero una sequenza di L=M/N numeri estratti casualmente
		unif_values_block[i]=averageSlice(uniform_sampling, i*L, L);
		imp_values_block[i]=averageSlice(importance_sampling, i*L, L);

	}
	
	for(int i=0; i<N; i++){
		unif_avg[i]=averageSlice(unif_values_block,0,i+1); //Calcola la media dei primi i+1 esperimenti
		unif_avg2[i]=averageQuadSlice(unif_values_block,0,i+1); //Calcola la media dei quadrati dei primi i+1 esperimenti
		unif_err[i]=error(unif_avg,unif_avg2,i); //Calcola la varianza della media per i primi i+1 esperimenti
		
		imp_avg[i]=averageSlice(imp_values_block,0,i+1); //Calcola la media dei primi i+1 esperimenti
		imp_avg2[i]=averageQuadSlice(imp_values_block,0,i+1); //Calcola la media dei quadrati dei primi i+1 esperimenti
		imp_err[i]=error(imp_avg,imp_avg2,i); //Calcola la varianza della media per i primi i+1 esperimenti
	}
	
	//Output dei risultati
	ofstream fileout;
	fileout.open("out_unif_02_1.txt");
	if(fileout.is_open()){
		fileout<<"N,\tAverage,\tError"<<endl;
		for(int i=0; i<N; i++){
			fileout<<i+1<<",\t"<<unif_avg[i]<<",\t"<<unif_err[i]<<endl;
		}
	}else cerr<< "Unable to open output file"<<endl;
	
	fileout.close();
	
	fileout.open("out_imp_02_1.txt");
	if(fileout.is_open()){
		fileout<<"N,\tAverage,\tError"<<endl;
		for(int i=0; i<N; i++){
			fileout<<i+1<<",\t"<<imp_avg[i]<<",\t"<<imp_err[i]<<endl;
		}
	}else cerr<< "Unable to open output file"<<endl;
	
	fileout.close();
	
	return 0;
}
