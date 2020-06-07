#include "genetics.h"
#include "functions.h"
#include "tsp.h"
#include<iostream>
#include<fstream>

using namespace std;

int main(){
    Random *rnd=new Random();
	initRandom(*rnd); //Sposto dal main l'inizializzazione del generatore (seed ecc.)

    //Dummy initial values, they are replaced when reading input.dat
    bool is_cycle=false;
    int population_size=100;
    int ntowns=32;
    int ngenerations=1000;
    double edge=10;
    double selection_parameter=2;
    double crossover_probability=0.75;
    double mutation_probability=0.05;

    ifstream ReadInput;

    ReadInput.open("input.dat");
    if(ReadInput.is_open()){
        ReadInput >> is_cycle; //1 for random towns on cycle, 0 for random towns in square
        ReadInput >> edge; //radius for the cycle, edge for the square
        ReadInput >> ntowns;
        ReadInput >> population_size;
        ReadInput >> selection_parameter;
        ReadInput >> crossover_probability;
        ReadInput >> mutation_probability;
        ReadInput >> ngenerations;
    }else{
        cerr<<"Unable to open input.dat"<<endl;
        exit(-1);
    } 
    ReadInput.close();
    
    //Random coordinates
    cout<<"TSP with genetic algorithms"<<endl<<endl;
    vector<pair<double,double>> coordinates; 
    if(is_cycle){
        coordinates=townsOnCycle(ntowns,edge,rnd);
        cout<<"Towns on a cycle"<<endl;
        cout<<"Radius: "<<edge<<endl;
    }else{
        coordinates=townsInSquare(ntowns,edge,rnd);
        cout<<"Towns in a square"<<endl;
        cout<<"Edge: "<<edge<<endl;
    }
    cout<<"Random coordinates printed in coordinates.out"<<endl<<endl;;
    cout<<"Number of towns: "<<ntowns<<endl;
    cout<<"Number of chromosomes: "<<population_size<<endl;
    cout<<"Selection parameter: "<<selection_parameter<<endl;
    cout<<"Mutations probability: "<<mutation_probability<<endl;
    cout<<"Crossover probability: "<<crossover_probability<<endl;
    cout<<"Number of generations: "<<ngenerations<<endl<<endl;

    //Print the random coordinates
    ofstream coordinates_stream("coordinates.out");
    for(auto& town : coordinates){
        coordinates_stream<<town.first<<", "<<town.second<<endl;
    }
    coordinates_stream.close();

    //ntowns-1 genes are needed, the first town is fixed
    PopulationTSP pop(population_size,ntowns-1, rnd, coordinates);

    
    ofstream chromosomes_stream, length_stream, bestlength_stream;
    chromosomes_stream.open("fittest.out");
    length_stream.open("path_length.out");
    bestlength_stream.open("best_length.out");

    for(int i=0; i<ngenerations; i++){
        pop.newGeneration(selection_parameter, crossover_probability, mutation_probability);
        length_stream<<pop.averageDistanceBestHalf()<<endl;
        chromosomes_stream<<0<<", ";
        for(auto& gene : pop[0]){
            chromosomes_stream<<gene<<", ";
        }
        chromosomes_stream<<0<<endl;
        bestlength_stream<<1./pop[0].getFitness()<<endl;
        if(i%(ngenerations/10)==0){
            cout<<"Generation: "<<i<<"/"<<ngenerations<<endl;
        }
    }
    chromosomes_stream.close();
    length_stream.close();
    bestlength_stream.close();
    
    cout<<"================================="<<endl;
    cout<<endl<<"Final average of the path length for the best half of the population: "<<endl;
    cout<<pop.averageDistanceBestHalf()<<endl<<endl;
    cout<<endl<<"Length of the final best path: "<<endl;
    cout<<1./pop[0].getFitness()<<endl<<endl;
    cout<<"================================="<<endl;


    cout<<"Fittest chromosome in each generation: fittest.out"<<endl;
    cout<<"Path length of the fittest chromosome in each generation: best_length.out"<<endl;
    cout<<"Evolution of path length for the best half of the chromosomes: path_length.out"<<endl;

    return 0;
}