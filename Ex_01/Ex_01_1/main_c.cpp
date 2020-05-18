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

	int M=100; //Numero di sottointervalli
	int n=10000; //Numero di lanci
	int n_repetitions=100; //Quante volte calcolo il chi-quadro

	vector<double> chi(n_repetitions);
	vector<double> extracted_numbers(n);

	cout<<"Chi-square test"<<endl;
	cout<<"Number of sub-intervals in [0,1): "<<M<<endl;
	cout<<"Extracted numbers for each chi-square calculation: "<<n<<endl;
	cout<<"Number of chi-square calculated: "<<n_repetitions<<endl;

	for(int j=0; j<n_repetitions; j++){
		for(int i=0; i<n; i++){
			//Estraggo gli i-esimi n numeri casuali
			extracted_numbers[i]=rnd.Rannyu();
		}
		//Con createHistogram divido i dati in extracted numbers in M canali, il secondo argomento è la lunghezza dell'intervallo a cui appartengono i dati in input
		//Il vettore così ottenuto si usa per calcolare il chi quadro (in questo caso considero un numero di occorrenze atteso n/M uguale per tutti i sottointervalli).
		chi[j]=calculateChiSquare(createHistogram(extracted_numbers,1,M),double(n)/M);
	}
	rnd.SaveSeed();



	//Output dei risultati
	ofstream fileout;
	fileout.open("out_01_1c.txt");
	cout<<"Output file: out_01_1c.txt"<<endl;
	if(fileout.is_open()){
		fileout<<"j,\tChiSquare"<<endl;
		for(int j=0; j<n_repetitions; j++){
			fileout<<j+1<<",\t"<<chi[j]<<endl;
		}
	}else cerr<< "Unable to open output file"<<endl;

	fileout.close();
	return 0;
}
