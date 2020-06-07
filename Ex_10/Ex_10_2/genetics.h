#ifndef __Genetics__
#define __Genetics__

#include <iterator>
#include<iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include "random.h"

using namespace std;
 

template <typename T> class Chromosome{
protected:
    vector<T> m_genes;
    double m_fitness;
public:
    Chromosome(int size){ m_genes=vector<T>(size); m_fitness=0; };

    T& operator [](int i) { return m_genes[i]; }; //access the "raw" vector
    const T& operator [](int i) const{ return m_genes[i]; }; //access the "raw" vector

    void setFitness(double fitness){ m_fitness = fitness; };
    double getFitness() const{ return m_fitness; };
    
    bool operator < (const Chromosome<T>& chromosome) const{
        return (getFitness() < chromosome.getFitness());
    }

    void swapGenes(int gene1, int gene2);
    void swapBlockOfGenes(int firstGene1, int firstGene2, int block_size);
    void shiftGenes(int gene1, int genes_shifted, int shift);
    void reverseGenes(int gene1, int genes_reversed);

    int size() const{return m_genes.size();};
    typename vector<T>::iterator begin(){ return m_genes.begin(); };
    typename vector<T>::const_iterator begin() const{ return m_genes.end(); };
    typename vector<T>::iterator end(){ return m_genes.end(); };
    typename vector<T>::const_iterator end() const{ return m_genes.end(); };
};

template <typename T> class Population{
protected:
    vector<Chromosome<T>> m_chromosomes;
    int m_population_size;
    int m_chromosomes_size;
    int m_generation_number;
    Random *m_rnd;
public:
    Population(int population_size, int chromosomes_size, Random *rnd);
    Chromosome<T>& operator [](int i) { return m_chromosomes[i]; }; //access the "raw" vector
    const Chromosome<T>& operator [](int i) const { return m_chromosomes[i]; }; //access the "raw" vector

    virtual void calculateFitness(Chromosome<T> &chromosome)=0; //It depends on the implementation
    void sortByFitness();
    double averageFitness() const;

    Chromosome<T>& randomChromosome();
    T& randomGene(Chromosome<T>& chromosome);
    int randomChromosomeIndex() const;
    int randomGeneIndex() const;

    Chromosome<T> rankSelection(double p);
    Chromosome<T> tournamentSelection(int tournament_size);

    void replaceChromosome(int replaced_index, const vector<T>& new_chromosome);

    void newGeneration(double parameter_selection, double probability_crossover, double probability_mutation);
    int generationNumber() const{ return m_generation_number;};
    virtual pair<Chromosome<T>,Chromosome<T>> crossover(Chromosome<T>& chromosome1, Chromosome<T>& chromosome2, int site_cut) const;

    int size() const{return m_chromosomes.size();}
    typename vector<Chromosome<T>>::iterator begin(){ return m_chromosomes.begin(); };
    typename vector<Chromosome<T>>::const_iterator begin() const{ return m_chromosomes.end(); };
    typename vector<Chromosome<T>>::iterator end(){ return m_chromosomes.end(); };
    typename vector<Chromosome<T>>::const_iterator end() const{ return m_chromosomes.end(); };
};

// =====================================================
// CHROMOSOMES
// =====================================================

template <typename T>
void Chromosome<T> :: swapGenes(int gene1, int gene2){
    swap(m_genes[gene1], m_genes[gene2]);
}

template <typename T>
void Chromosome<T> :: shiftGenes(int gene1, int genes_shifted, int shift){
    vector<T> temp=vector<T>(genes_shifted);
    for(int i=0; i<genes_shifted; i++){
        int index_copy = (i+gene1)%size(); //periodic boundary conditions
        temp[i]=m_genes[index_copy];
    }

    for(int i=0; i<shift; i++){
        m_genes[(i+gene1)%size()]=m_genes[(i+gene1+genes_shifted)%size()];
    }
    
    for(int i=0; i<genes_shifted; i++){
        m_genes[(i+gene1+shift)%size()]=temp[i];
    }
}

template <typename T>
void Chromosome<T> :: swapBlockOfGenes(int firstGene1, int firstGene2, int block_size){
    for(int i=0; i<block_size; ++i){
        swap(m_genes[(firstGene1+i)%size()],m_genes[(firstGene2+i)%size()]);
    }
}
template <typename T>
void Chromosome<T> :: reverseGenes(int gene1, int genes_reversed){
    for(int i=0; i<genes_reversed/2; i++){
        swap(m_genes[(gene1+i)%size()], m_genes[(gene1+genes_reversed-i-1)%size()]);
    }
}
   
// =====================================================
// POPULATION
// =====================================================

template <typename T>
Population<T> :: Population(int population_size, int chromosomes_size, Random *rnd){
    m_chromosomes=vector<Chromosome<T>>(population_size, Chromosome<T>(chromosomes_size));
    m_rnd=rnd;
    m_population_size=population_size;
    m_chromosomes_size=chromosomes_size;
    m_generation_number=0;
};

template <typename T>
double Population<T> :: averageFitness() const{
    double avg_fitness=0;
    for(auto& chromosome : m_chromosomes){
        avg_fitness+=chromosome.getFitness();
    }
    avg_fitness/=m_population_size;
    return avg_fitness;
}


template <typename T>
Chromosome<T>& Population<T> :: randomChromosome(){
    int random_index=m_rnd->UniformInteger(0,size()-1);
    return m_chromosomes[random_index];
};

template <typename T>
T& Population<T> :: randomGene(Chromosome<T>& chromosome){
    int random_index=m_rnd->UniformInteger(0,chromosome.size()-1);
    return chromosome[random_index];
};

template <typename T>
int Population<T> :: randomGeneIndex() const{
    return m_rnd->UniformInteger(0,m_chromosomes_size-1);
};

template <typename T>
int Population<T> :: randomChromosomeIndex() const{
    return m_rnd->UniformInteger(0,m_population_size-1);
};

template <typename T>
Chromosome<T> Population<T> :: rankSelection(double p){

    return m_chromosomes[int(m_population_size*pow(m_rnd->Rannyu(),p))];

};

template <typename T>
Chromosome<T> Population<T> :: tournamentSelection(int tournament_size){
    vector<Chromosome<T>> fighters(tournament_size, Chromosome<T>(m_chromosomes_size));
    for(auto& chromosome : fighters){
        chromosome = randomChromosome();
    }

    sort(fighters.rbegin(), fighters.rend());

    return fighters[0];
};

template <typename T>
pair<Chromosome<T>,Chromosome<T>> Population<T> :: crossover(Chromosome<T>& chromosome1, Chromosome<T>& chromosome2, int site_cut) const{
    pair<Chromosome<T>,Chromosome<T>> children (chromosome1, chromosome2);
    swap_ranges(children.first.begin(), children.first.begin()+site_cut, children.second.begin());
    return children;
};


template <typename T>
void Population<T> :: newGeneration(double parameter_selection, double probability_crossover, double probability_mutation){
    //m_chromosomes are already sorted by fitness
    
    vector<Chromosome<T>> new_chromosomes(m_population_size, Chromosome<T>(m_chromosomes_size));;
    
    for(int i=0; i<m_population_size; i++){
        new_chromosomes[i]=rankSelection(parameter_selection);
        //new_chromosomes[i]=tournamentSelection(int(parameter_selection)); //using tournamentselection it is not strictly needed to sort the whole population by fitness
    }
    
    for(int i=0; i<m_population_size/2; i++){
        //Choose pair of selected chromosomes for crossover
        int index1=i;
        int index2=m_population_size-1-i;
        if(m_rnd->Rannyu()<probability_crossover){
            pair<Chromosome<T>,Chromosome<T>> children = crossover(new_chromosomes[index1],new_chromosomes[index2],randomGeneIndex());
            new_chromosomes[index1]=children.first;
            new_chromosomes[index2]=children.second;
        }
    }

    for(auto& chromosome : new_chromosomes){

        if(m_rnd->Rannyu()<probability_mutation){
            chromosome.swapGenes(randomGeneIndex(), randomGeneIndex());
        }
        if(m_rnd->Rannyu()<probability_mutation){
            chromosome.swapBlockOfGenes(randomGeneIndex(), randomGeneIndex(), m_rnd->UniformInteger(0,m_chromosomes_size/2));
        }
        if(m_rnd->Rannyu()<probability_mutation){
            chromosome.shiftGenes(randomGeneIndex(), m_rnd->UniformInteger(0,m_chromosomes_size), m_rnd->UniformInteger(0,m_chromosomes_size));
        }
        if(m_rnd->Rannyu()<probability_mutation){
            chromosome.reverseGenes(randomGeneIndex(), m_rnd->UniformInteger(0,m_chromosomes_size));
        }
    }
   
    m_chromosomes=new_chromosomes;
    sortByFitness();
    m_generation_number++;
};

template <typename T>
void Population<T> :: sortByFitness(){ 

     for(auto& chromosome : m_chromosomes){
        calculateFitness(chromosome);
    }

    sort(m_chromosomes.rbegin(), m_chromosomes.rend()); //sort in descending order by fitness
}; 

template <typename T>
void Population<T> :: replaceChromosome(int replaced_index, const vector<T>& new_chromosome){
    for(int i=0; i<m_chromosomes_size; i++){
        m_chromosomes[replaced_index][i]=new_chromosome[i];
    }
    sortByFitness();
}




#endif //__Genetics__