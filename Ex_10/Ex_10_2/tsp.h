#include "genetics.h"
#include <cmath>

vector<pair<double,double>> townsOnCycle(int N, double R, Random *rnd){
    vector<pair<double,double>> coordinates(N);

    for(int i=0; i<N; i++){

        double x = rnd->Rannyu(-R,R);
        double y = rnd->Rannyu(-R,R);
        double norm=sqrt(x*x+y*y);

        while(norm > 1){
            x = rnd->Rannyu(-R,R);
            y = rnd->Rannyu(-R,R);
            norm=sqrt(x*x+y*y);
        }

        coordinates[i].first=R*x/norm;
        coordinates[i].second=R*y/norm;
    }
    return coordinates;
}

vector<pair<double,double>> townsInSquare(int N, double edge, Random *rnd){
    vector<pair<double,double>> coordinates(N);
    for(int i=0; i<N; i++){
        coordinates[i].first=rnd->Rannyu(-edge/2,edge/2);
        coordinates[i].second=rnd->Rannyu(-edge/2,edge/2);
    }
    return coordinates;
}

class PopulationTSP : public Population<int>{
private:
    vector<pair<double,double>> m_coordinates;
    double distance(int town1, int town2) const;
public:
    PopulationTSP(int population_size, int chromosomes_size, Random *rnd, vector<pair<double,double>> coordinates);
    void calculateFitness(Chromosome<int> &chromosome); //It depends on the implementation
    double averageDistanceBestHalf() const;
    pair<Chromosome<int>,Chromosome<int>> crossover(Chromosome<int>& chromosome1, Chromosome<int>& chromosome2, int site_cut) const;

};

PopulationTSP :: PopulationTSP(int population_size, int chromosomes_size, Random *rnd, vector<pair<double,double>> coordinates) :
 Population<int>(population_size, chromosomes_size, rnd){
    m_coordinates = coordinates;
    int counter=1;
    //fill every chromosome with every town 1, 2, 3, ...
    for(auto& chromosome : m_chromosomes){
        for(auto& gene : chromosome){
            gene=counter;
            counter++;
        }
        counter=1;
    } 
   
    //shuffle the genes of every chromosome randomly
    for(auto& chromosome : m_chromosomes){
        for(int i=0; i<chromosomes_size; i++){
            chromosome.swapGenes(randomGeneIndex(), randomGeneIndex());
        }
    }
    sortByFitness();  //it includes calculateFitness() for each chromosome   
};

double PopulationTSP :: distance(int town1, int town2) const{
    return sqrt(pow(m_coordinates[town2].first-m_coordinates[town1].first,2)+pow(m_coordinates[town2].second-m_coordinates[town1].second,2));
}

void PopulationTSP :: calculateFitness(Chromosome<int> &chromosome){
    double total_distance=0;
    int last_gene=0; //the first town is always 0
    for(auto& gene : chromosome){
        total_distance+=distance(last_gene, gene);
        last_gene=gene;
    }
    total_distance+=distance(last_gene, 0); //the last town is the first town
    chromosome.setFitness(1./total_distance);
}

double PopulationTSP :: averageDistanceBestHalf() const{
    double avg_distance=0;
    for(auto& chromosome : m_chromosomes){
        avg_distance+=1./chromosome.getFitness();
    }
    avg_distance/=m_population_size;
    return avg_distance;
}

pair<Chromosome<int>,Chromosome<int>> PopulationTSP :: crossover(Chromosome<int>& chromosome1, Chromosome<int>& chromosome2, int site_cut) const{
    pair<Chromosome<int>,Chromosome<int>> children (chromosome1, chromosome2);

    int index=0;
    for(auto& elem2 : chromosome2){
        //for every element in the second chromosome, check if it is contained in the cut part of the first chromosome
        if(find(chromosome1.begin()+site_cut, chromosome1.end(), elem2) != chromosome1.end()){
            //if it is, place it in the first child chromosome
            children.first[site_cut+index]=elem2;
            index++;
        }
    }

    //do the same for the second child
    index=0;
    for(auto& elem1 : chromosome1){
        //for every element in the first chromosome, check if it is contained in the cut part of the second chromosome
        if(find(chromosome2.begin()+site_cut, chromosome2.end(), elem1) != chromosome2.end()){
            //if it is, place it in the first child chromosome
            children.second[site_cut+index]=elem1;
            index++;
        }
    }

    return children;
};






bool checkChromosome(const Chromosome<int>& chromosome){
    Chromosome<int> test_vector(chromosome.size());
    for(int i=0; i<test_vector.size(); i++){
        test_vector[i]=i+1;
    }
    return is_permutation(chromosome.begin(), chromosome.end(), test_vector.begin());
}

bool checkPopulation(const PopulationTSP& population){
    for(auto& chromosome : population){
        if(!checkChromosome(chromosome)){
            return false;
        }
    }
    return true;
}


