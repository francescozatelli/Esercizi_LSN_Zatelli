#include "genetics.h"
#include "functions.h"
#include "annealer.h"
#include<iostream>
#include<fstream>

using namespace std;

int main(){
    Random *rnd=new Random();
	initRandom(*rnd); //Sposto dal main l'inizializzazione del generatore (seed ecc.)

    bool is_cycle;
    double edge, beta0, beta_coeff;
	int ntowns, nsteps_beta, ntrials_per_beta;

    ifstream ReadInput;

    ReadInput.open("input.dat");
    if(!ReadInput.is_open()){
        cerr<<"Unable to open input.dat"<<endl;
        exit(-1);
    }
    ReadInput >> is_cycle; //1 for random towns on cycle, 0 for random towns in square
    ReadInput >> edge; //radius for the cycle, edge for the square
    ReadInput >> ntowns;
    ReadInput>>beta0;
    ReadInput>>beta_coeff;
    ReadInput>>nsteps_beta;
    ReadInput>>ntrials_per_beta;
    ReadInput.close();

    //Random coordinates
    cout<<"TSP with simulated annealing"<<endl<<endl;
    vector<pair<double,double>> coordinates; 
    if(is_cycle){
        coordinates=townsOnCycle(32,edge,rnd);
        cout<<"Towns on a cycle"<<endl;
        cout<<"Radius: "<<edge<<endl;
    }else{
        coordinates=townsInSquare(32,edge,rnd);
        cout<<"Towns in a square"<<endl;
        cout<<"Edge: "<<edge<<endl;
    }
    cout<<"Random coordinates printed in coordinates.out"<<endl<<endl;;
    cout<<"Number of towns: "<<ntowns<<endl;
    cout<<"Initial beta: "<<beta0<<" - Initial T: "<<1./beta0<<endl;
	cout<<"Beta increase coefficient: "<<beta_coeff<<endl;
	cout<<"Number of beta steps: "<<nsteps_beta<<endl;
	cout<<"Final beta: "<<beta0*pow(beta_coeff, nsteps_beta)<<endl;
    cout<<"Number of trial chromosomes for each beta: "<<ntrials_per_beta<<endl;

    //Print the random coordinates
    ofstream coordinates_stream("coordinates.out");
    for(auto& town : coordinates){
        coordinates_stream<<town.first<<", "<<town.second<<endl;
    }
    coordinates_stream.close();


    ofstream fileout, chromosomes_stream;

	fileout.open("annealing.out");
	fileout<<"Beta,\t T,\t Length,\tMutationAcceptance"<<endl;

    chromosomes_stream.open("fittest.out");


    //ntowns-1 genes are needed, the first town is fixed, initial random chromosome
	Annealer annealer(ntowns-1,rnd, coordinates);
    double beta=beta0;

    //Print initial data
	for(int i=0; i<=nsteps_beta; i++){
		cout<<"======================"<<endl;
		cout<<"Beta = "<<beta<<"; T = "<<1./beta<<endl;
		annealer.annealingStep(beta, ntrials_per_beta);
		cout<<"Acceptance for a trial chromosome: "<<annealer.getAcceptance()<<endl;
		double length = annealer.getLength();
		cout<<"Length: "<<length<<endl;
		fileout<<beta<<",\t"<<1./beta<<",\t"<<length<<",\t"<<annealer.getAcceptance()<<endl;

        chromosomes_stream<<0<<", ";
        for(auto& gene : annealer.getChromosome()){
            chromosomes_stream<<gene<<", ";
        }
        chromosomes_stream<<0<<endl;
        beta = beta*beta_coeff; 
	}
	cout<<endl<<"========== Finished sampling =============="<<endl;
	rnd->SaveSeed();
	fileout.close();
    chromosomes_stream.close();

    return 0;
}