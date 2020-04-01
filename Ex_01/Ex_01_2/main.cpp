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
	

	vector<double> avg_std(M); //Vettore che conterrà i valori della variabile aleatoria S_N per il dado standard
	vector<double> avg_exp(M); //Vettore che conterrà i valori della variabile aleatoria S_N per il dado esponenziale
	vector<double> avg_cauchy(M); //Vettore che conterrà i valori della variabile aleatoria S_N per il dado Cauchy

	
	//STANDARD DICE
	for(int i=0;i<100;i++){
		for(int j=0; j<M; j++){
			double unif_rand=6*rnd.Rannyu(); //Estraggo un numero in [0,6)
			//Arrotondo per eccesso al primo intero maggiore o uguale al numero estratto
			//Media dato per dato, per ogni ciclo su i aggiungo un dato in più alla media.
			avg_std[j]=(i*avg_std[j]+ceil(unif_rand))/double(i+1);
		}
		rnd.SaveSeed();
		
		if(i==0){ //N=1
			printDice(avg_std,"out_std_S1.txt");
		}else if(i==1){ //N=2
			printDice(avg_std,"out_std_S2.txt");
		}else if(i==9){ //N=10
			printDice(avg_std,"out_std_S10.txt");
		}else if(i==99){ //N=100
			printDice(avg_std,"out_std_S100.txt");
		}
	}

	//DADO ESPONENZIALE
	for(int i=0;i<100;i++){
		
		for(int j=0; j<M; j++){
			//Media dato per dato, per ogni ciclo su i aggiungo un dato in più alla media.
			avg_exp[j]=(i*avg_exp[j]+rnd.Exp(1.))/double(i+1);
		}
		rnd.SaveSeed();
		
		if(i==0){ //N=1
			printDice(avg_exp,"out_exp_S1.txt");
		}else if(i==1){ //N=2
			printDice(avg_exp,"out_exp_S2.txt");
		}else if(i==9){ //N=10
			printDice(avg_exp,"out_exp_S10.txt");
		}else if(i==99){ //N=100
			printDice(avg_exp,"out_exp_S100.txt");
		}
	}

	//DADO CAUCHY
	for(int i=0;i<100;i++){
		
		for(int j=0; j<M; j++){
			//Media dato per dato, per ogni ciclo su i aggiungo un dato in più alla media.
			avg_cauchy[j]=(i*avg_cauchy[j]+rnd.Cauchy(0.,1.))/double(i+1);
		}
		rnd.SaveSeed();
		
		if(i==0){ //N=1
			printDice(avg_cauchy,"out_cauchy_S1.txt");
		}else if(i==1){ //N=2
			printDice(avg_cauchy,"out_cauchy_S2.txt");
		}else if(i==9){ //N=10
			printDice(avg_cauchy,"out_cauchy_S10.txt");
		}else if(i==99){ //N=100
			printDice(avg_cauchy,"out_cauchy_S100.txt");
		}
	}
	
	return 0;
}
