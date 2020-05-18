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

	int M=100000; //Numero di lanci
	int N=100; //Numero di blocchi
	int K=M/N; //Numero di numeri casuali da considerare in ogni blocco

	double d=1.; //Distanza tra le righe
	double L=0.8; //Lunghezza ago

	vector<double> extracted_centres(M);
	vector<double> extracted_thetas(M);


	vector<double> values_block(N); //Singole "misure" per ognuno degli N blocchi/esperimenti
	vector<double> avg(N); //Medie delle misure per ogni numero di esperimenti considerati da 1 a N
	vector<double> avg2(N); //Medie dei quadratidelle misure per ogni numero di esperimenti considerati da 1 a N
	vector<double> err(N); //Deviazione standard della media per ogni numero di esperimenti considerati da 1 a N

	cout<<"Buffon experiment"<<endl;
	cout<<"Needle length: "<<L<<endl;
	cout<<"Distance between two lines: "<<d<<endl;
	cout<<"Number of throws: "<<M<<endl;
	cout<<"Number of throws per block: "<<K<<endl;

	for(int i=0; i<M; i++){
		//Scelgo la cruna dell'ago in [0,d)
		extracted_centres[i]=d*rnd.Rannyu();

		double x=0;
		double y=0;
		do{
			//Estraggo due numeri casuali uniformi in [0,1)x[0,1)
			x=-1+2*rnd.Rannyu();
			y=-1+2*rnd.Rannyu();
			//Se uno dei due è zero o se il modulo è maggiore di 1, riestraggo la coppia di numeri
		}while((x==0 and y==0) or x*x+y*y>1);

		//Atan2 restituisce l'angolo delle coordinate polari
		extracted_thetas[i]=atan2(y,x);
	}
	rnd.SaveSeed();

	//Calcolo i valori per gli N blocchi (le N "misure")
	for(int i=0; i<N; i++){
		int n_hit=0;
		//Per ogni blocco considero una sequenza di K=M/N numeri estratti casualmente
		for(int j=i*K; j<(i+1)*K; j++){
			//Calcolo la proiezione della posizione della punta dell'ago
			double point=extracted_centres[j]+L*cos(extracted_thetas[j]);
			if(point<0 or point>d){
				n_hit++;
			}
		}
		//Stima di pi greco per un esperimento, K=N_throws
		values_block[i]=2*L*K/(n_hit*d);
	}

	for(int i=0; i<N; i++){
		avg[i]=averageSlice(values_block,0,i+1); //Calcola la media dei primi i+1 esperimenti
		avg2[i]=averageQuadSlice(values_block,0,i+1); //Calcola la media dei quadrati dei primi i+1 esperimenti
		err[i]=error(avg,avg2,i); //Calcola la deviazione standard della media per i primi i+1 esperimenti
	}

	//Output dei risultati
	ofstream fileout;
	fileout.open("out_01_3.txt");
	cout<<"Output file: out_01_3.txt"<<endl;
	if(fileout.is_open()){
		fileout<<"N,\tAverage,\tError"<<endl;
		for(int i=0; i<N; i++){
			fileout<<i+1<<",\t"<<avg[i]<<",\t"<<err[i]<<endl;
		}
	}else cerr<< "Unable to open output file"<<endl;

	fileout.close();
	return 0;
}
