#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "functions.h"
#include "rw.h"

using namespace std;

int main (int argc, char *argv[]){

	Random *rnd=new Random();
	initRandom(*rnd); //Sposto dal main l'inizializzazione del generatore (seed ecc.)
	int M=10000; //Numero di random walks
	int N=100; //Numero di blocchi
	int L=M/N; //Numero di random walks per blocco
	int nsteps=100; //Numero di passi per ogni random walk
	int a_discrete=1; //Distanza reticolare (caso discreto)
	double a_continuous=1.; //Distanza percorsa per step (caso continuo)

	//Valori dei blocchi di L random walks nel caso di spazio discreto e continuo in funzione del numero di step
	//(Per ogni blocco conservo l'andamento da i=0 a i=100)
	vector<vector <double>> discr_squared_distances_block(N,vector<double>(nsteps));
	vector<vector <double>> cont_squared_distances_block(N,vector<double>(nsteps));

	vector<double> discr_avg(nsteps); //Valor medio dei quadrati della distanza percorsa in funzione del numero di step
	vector<double> discr_avg2(nsteps); //Valor medio dei quadrati dei quadrati della distanza percorsa in funzione del numero di step
	vector<double> discr_err(nsteps); //Deviazione standard della media in funzione del numero di step

	vector<double> cont_avg(nsteps); //Come sopra
	vector<double> cont_avg2(nsteps);
	vector<double> cont_err(nsteps);

	RandomWalk<int> rw_discrete(rnd,a_discrete);
	RandomWalk<double> rw_continuous(rnd,a_continuous);

	cout<<"Random walks"<<endl;
	cout<<"Number of random walks: "<<M<<endl;
	cout<<"Number of blocks: "<<N<<endl;
	cout<<"Number of random walks per block: "<<L<<endl;
	cout<<"Number of steps for each random walk: "<<nsteps<<endl;
	cout<<"Space step for discrete random walk a = "<<a_discrete<<endl;
	cout<<"Space step for continuous random walk a = "<<a_continuous<<endl;


	/*Questa parte di codice è un po' contorta ma serve per evitare di creare M vettori contenenti nsteps elementi.
	*	Praticamente determino i valori dei singoli blocchi a mano a mano che simulo gli M random walks.
	*/
	for(int k=0; k<N;k++){ //Ciclo su tutti i blocchi di random walk (sugli "esperimenti")
		for(int j=0; j<L; j++){
			//Ogni ciclo su j corrisponde a un singolo random walk che andrò a mediare facendolo confluire in un singolo blocco
			for(int i=0; i<nsteps; i++){
				//Ogni ciclo su i corrisponde a uno step del random walk
				rw_discrete.makeDiscreteStep();
				rw_continuous.makeContinuousStep();
				//Aggiorno l'i-esimo step del k-esimo blocco di random walks aggiungendo il j-esimo random walk
				discr_squared_distances_block[k][i]=j/double(j+1)*discr_squared_distances_block[k][i]+1/double(j+1)*pow(rw_discrete.getDistance(),2);
				cont_squared_distances_block[k][i]=j/double(j+1)*cont_squared_distances_block[k][i]+1/double(j+1)*pow(rw_continuous.getDistance(),2);

			}
			rw_discrete.setPosition(0,0,0);
			rw_continuous.setPosition(0.,0.,0.);
		}
	}

	rnd->SaveSeed();

	//CASO DISCRETO
	for(int i=0; i<nsteps; i++){
		double sum=0;
		double sum2=0;
		for(int j=0; j<N; j++){
			sum+=discr_squared_distances_block[j][i];
			sum2+=pow(discr_squared_distances_block[j][i],2);
		}
		discr_avg[i]=sum/N; //<|r_N|^2>
		discr_avg2[i]=sum2/N; //<|r_N|^4>
		discr_err[i]=sqrt((discr_avg2[i]-discr_avg[i]*discr_avg[i])/(N-1)); //Errore su <|r_N|^2>
	}

	//CASO CONTINUO
	for(int i=0; i<nsteps; i++){
		double sum=0;
		double sum2=0;
		for(int j=0; j<N; j++){
			sum+=cont_squared_distances_block[j][i];
			sum2+=pow(cont_squared_distances_block[j][i],2);
		}
		cont_avg[i]=sum/N;  //<|r_N|^2>
		cont_avg2[i]=sum2/N; //<|r_N|^4>
		cont_err[i]=sqrt((cont_avg2[i]-cont_avg[i]*cont_avg[i])/(N-1)); //Errore su <|r_N|^2>
	}


	//Output dei risultati
	ofstream fileout;
	fileout.open("out_discr.txt");
	cout<<"Output file for discrete random walk: out_discr.txt"<<endl;
	if(fileout.is_open()){
		fileout<<"N,\tAverage,\tError"<<endl;
		fileout<<"0,\t0,\t0"<<endl; //Il primo step lo metto manualmente in quanto banale
		for(int i=0; i<N; i++){
			fileout<<i+1<<",\t"<<discr_avg[i]<<",\t"<<discr_err[i]<<endl;
		}
	}else cerr<< "Unable to open output file"<<endl;

	fileout.close();

	fileout.open("out_cont.txt");
	cout<<"Output file for continuous random walk: out_cont.txt"<<endl;
	if(fileout.is_open()){
		fileout<<"N,\tAverage,\tError"<<endl;
		fileout<<"0,\t0,\t0"<<endl; //Il primo step lo metto manualmente in quanto banale
		for(int i=0; i<N; i++){
			fileout<<i+1<<",\t"<<cont_avg[i]<<",\t"<<cont_err[i]<<endl;
		}
	}else cerr<< "Unable to open output file"<<endl;

	fileout.close();

	return 0;
}
