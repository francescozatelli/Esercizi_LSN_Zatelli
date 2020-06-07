#include "genetics.h"
#include "functions.h"
#include "tsp.h"
#include "mpi.h"
#include<iostream>
#include<fstream>
#include <string> 

using namespace std;

int main(){
    MPI_Init(NULL, NULL);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double tstart = MPI_Wtime();
    Random *rnd=new Random();
    initRandomParallel(*rnd, rank);
    
    int is_cycle;
    int population_size, ntowns, ngenerations, nmigration;
    double edge, selection_parameter, crossover_probability, mutation_probability;

    if(rank == 0){ //Only the first node reads the external output
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
            ReadInput >> nmigration;
        }else cerr<<"Unable to open input.dat"<<endl;
        ReadInput.close();
    }

    //Broadcast all the variables to the nodes
    // (the data types in the slides didn't work, lead to "MPI_ERR_TYPE: invalid datatype")
    MPI_Bcast(&is_cycle,1,MPI_INT,0, MPI_COMM_WORLD);
    MPI_Bcast(&population_size,1,MPI_INT,0, MPI_COMM_WORLD);
    MPI_Bcast(&ntowns,1,MPI_INT,0, MPI_COMM_WORLD);
    MPI_Bcast(&ngenerations,1,MPI_INT,0, MPI_COMM_WORLD);
    MPI_Bcast(&nmigration,1,MPI_INT,0, MPI_COMM_WORLD);
    MPI_Bcast(&edge,1,MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Bcast(&selection_parameter,1,MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Bcast(&crossover_probability,1,MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Bcast(&mutation_probability,1,MPI_DOUBLE,0, MPI_COMM_WORLD);

    vector<pair<double,double>> coordinates(ntowns); 
    vector<double> x(ntowns);
    vector<double> y(ntowns);
    if(rank == 0){ //Coordinates chosen by the first node
        cout<<"Parallelized TSP with genetic algorithms, rank "<<rank<<" out of "<<size<<" is speaking"<<endl<<endl; 
        if(is_cycle==1){
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
        cout<<"Number of generations: "<<ngenerations<<endl;
        cout<<"Migration every "<<nmigration<<" generations"<<endl<<endl;
        
        //Print the chosen coordinates
        ofstream coordinates_stream("coordinates.out");
        for(auto& town : coordinates){
            coordinates_stream<<town.first<<", "<<town.second<<endl;
        }
        coordinates_stream.close();

        for(int i=0; i<ntowns; i++){ //unpack the coordinates in order to broadcast array of double
            x[i]=coordinates[i].first;
            y[i]=coordinates[i].second;
        }
    }
    
    MPI_Bcast(&ngenerations,1,MPI_INT,0, MPI_COMM_WORLD);

    //Broadcast x and y array separately
    MPI_Bcast(&x[0],ntowns,MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Bcast(&y[0],ntowns,MPI_DOUBLE,0, MPI_COMM_WORLD);

    //Pack again the coordinates into pairs (needed for the implementation of the GA)
    for(int i=0; i<ntowns; i++){
        coordinates[i] = tie(x[i],y[i]);
    }

   
    vector<int> gathered_array; //not needed for rank!=0
    if(rank==0){
        gathered_array = vector<int>((ntowns-1)*size);
    }

    vector<int> scattered_array(ntowns-1);

    ofstream chromosomes_stream, length_stream, bestlength_stream;
    chromosomes_stream.open("fittest.out"+to_string(rank));
    length_stream.open("path_length.out"+to_string(rank));
    bestlength_stream.open("best_length.out"+to_string(rank));

    //Actual genetic evolution
    PopulationTSP pop(population_size,ntowns-1, rnd, coordinates);
    for(int i=0; i<ngenerations; i++){
        pop.newGeneration(selection_parameter, crossover_probability, mutation_probability);

        length_stream<<pop.averageDistanceBestHalf()<<endl;
        chromosomes_stream<<0<<", ";
        for(auto& gene : pop[0]){
            chromosomes_stream<<gene<<", ";
        }
        chromosomes_stream<<0<<endl;
        bestlength_stream<<1./pop[0].getFitness()<<endl;

        if(i%nmigration==0 and i!=0){
            MPI_Gather(&pop[0][0], ntowns-1, MPI_INT, &gathered_array[0], ntowns-1, MPI_INT, 0, MPI_COMM_WORLD);
            if(rank == 0){

                vector<vector<int>> separated_gathered_array(size, vector<int>(ntowns-1));

                //Separate the gathered array
                for(int j=0; j<(ntowns-1)*size; j++){
                    separated_gathered_array[j/((ntowns-1))][j%(ntowns-1)]=gathered_array[j];
                }

                //Random permutation of the best chromosomes gathered 
                for(int j=0; j<size; j++){
                    int index1 = rnd->UniformInteger(0,size-1);
                    int index2 = rnd->UniformInteger(0,size-1);
                    while(index2 == index1){
                        index2 = rnd->UniformInteger(0,size-1);
                    }
                    swap(separated_gathered_array[index1],separated_gathered_array[index2]);
                }
  
                //Reinsert the shuffled chromosomes into the monodimensional vector
                for(int j=0; j<(ntowns-1)*size; j++){
                    gathered_array[j]=separated_gathered_array[j/((ntowns-1))][j%(ntowns-1)];
                }
                /*
                for(int j=0;j <size; j++){
                    cout<<"Chromosome shuffled "<<j<<endl;
                    for(int k=0; k<ntowns-1; k++){
                        cout<<gathered_array[k+j*(ntowns-1)]<<" - ";
                        //cout<<separated_gathered_array[j][k]<<" - ";
                    }
                    cout<<endl;
                }*/
            }


            //Scatter the shuffled chromosomes into all the nodes
            MPI_Scatter(&gathered_array[0], ntowns-1, MPI_INT, &scattered_array[0], ntowns-1, MPI_INT, 0, MPI_COMM_WORLD);

            //Replace the best chromosome (index 0) in each population with the scattered one            
            pop.replaceChromosome(0, scattered_array);
        }
    }

    printf("Rank %d. Best final path: %f, Final average length of the best half: %f \n", rank, 1./pop[0].getFitness(), pop.averageDistanceBestHalf());
    double tend=MPI_Wtime();
    double dt=tend-tstart;
    printf("Rank %d. Time: %f \n", rank, dt);

    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}