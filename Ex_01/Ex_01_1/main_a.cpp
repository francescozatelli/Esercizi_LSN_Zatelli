#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "functions.h"

using namespace std;

int main (int argc, char *argv[]){

	Random rnd;
	initRandom(rnd); //Sposto dal main l'inizializzazione del generatore (seed ecc.)

	int M=10000; //Numero di estrazioni
	int N=100; //Numero di blocchi
	int L=M/N; //Numero di numeri casuali da considerare in ogni blocco

	vector<double> extracted_numbers(M);
	vector<double> values_block(N); //Singole "misure" per ognuno degli N blocchi/esperimenti
	vector<double> avg(N); //Medie delle misure per ogni numero di esperimenti considerati da 0 a N
	vector<double> avg2(N); //Medie dei quadratidelle misure per ogni numero di esperimenti considerati da 0 a N
	vector<double> err(N); //Deviazione standard della media per ogni numero di esperimenti considerati da 0 a N


  cout<<"Test average of uniform distribution in [0,1)"<<endl;
  cout<<"Random numbers extracted: "<<M<<endl;
  cout<<"Random numbers per block: "<<L<<endl;

	//Estraggo M numeri casuali estratti uniformemente in [0,1)
	for(int i=0; i<M; i++){
		extracted_numbers[i]=rnd.Rannyu();
	}
	rnd.SaveSeed();



	//Calcolo i valori per gli N blocchi (le N "misure")
	for(int i=0; i<N; i++){
		//Fa la media di una "fetta" di vettore che parte da i*L ed è lunga L
		//Per ogni blocco considero una sequenza di L=M/N numeri estratti casualmente
		values_block[i]=averageSlice(extracted_numbers, i*L, L);
	}

	for(int i=0; i<N; i++){
		avg[i]=averageSlice(values_block,0,i+1); //Calcola la media dei primi i+1 esperimenti
		avg2[i]=averageQuadSlice(values_block,0,i+1); //Calcola la media dei quadrati dei primi i+1 esperimenti
		err[i]=error(avg,avg2,i); //Calcola la varianza della media per i primi i+1 esperimenti
	}

	//Output dei risultati
	ofstream fileout;
	fileout.open("out_01_1a.txt");
  cout<<"Output file: out_01_1a.txt"<<endl;
	if(fileout.is_open()){
		fileout<<"N,\tAverage,\tError"<<endl;
		for(int i=0; i<N; i++){
			fileout<<i+1<<",\t"<<avg[i]<<",\t"<<err[i]<<endl;
		}
	}else cerr<< "Unable to open output file"<<endl;

	fileout.close();
	return 0;
}
