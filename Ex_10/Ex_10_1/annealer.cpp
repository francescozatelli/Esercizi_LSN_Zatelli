#include "annealer.h"


//Returns the boltzmann weight for a given energy
double Annealer :: BoltzmannWeight(double beta, double energy){
    return exp(-beta*energy);
}


double Annealer :: distance(int town1, int town2) const{
    return sqrt(pow(m_coordinates[town2].first-m_coordinates[town1].first,2)+pow(m_coordinates[town2].second-m_coordinates[town1].second,2));
}

double Annealer :: calculateLength(Chromosome<int> &chromosome){
    double total_distance=0;
    int last_gene=0; //the first town is always 0
    for(auto& gene : chromosome){
        total_distance+=distance(last_gene, gene);
        last_gene=gene;
    }
    total_distance+=distance(last_gene, 0); //the last town is the first town
    chromosome.setFitness(1./total_distance);
    return total_distance;
}


Annealer :: Annealer(int size, Random *rnd, vector<pair<double,double>> coordinates) : 
  m_chromosome(size),
  m_trial_chromosome(size)
  {
  m_rnd=rnd;
  m_coordinates = coordinates;
  int counter=1;
  //fill chromosome with every town 1, 2, 3, ...
  for(auto& gene : m_chromosome){
    gene=counter;
    counter++;
  }
  //shuffle the genes randomly
  for(int i=0; i<size; i++){
    m_chromosome.swapGenes(randomGeneIndex(), randomGeneIndex());
  }

  m_accepted=0;
  m_attempted=0;
  calculateLength(m_chromosome);
  calculateLength(m_trial_chromosome);
}

Annealer :: Annealer(Chromosome<int> chromosome, Random *rnd,vector<pair<double,double>> coordinates) : 
  m_chromosome(chromosome),
  m_trial_chromosome(chromosome)
  {
  m_coordinates = coordinates;
  m_rnd=rnd;
  m_accepted=0;
  m_attempted=0;
  calculateLength(m_chromosome);
  calculateLength(m_trial_chromosome);
}

void Annealer :: annealingStep(double beta, int nsteps){

  m_accepted=0;
  m_attempted=0;

  for(int i=0; i<nsteps; i++){
    m_trial_chromosome=m_chromosome;
    
    //Mutate the trial chromosome
    
    int choosemutation = m_rnd->UniformInteger(0,3);

    switch(choosemutation){
      case 0:
        m_trial_chromosome.swapGenes(randomGeneIndex(), randomGeneIndex());
      case 1:
        m_trial_chromosome.swapBlockOfGenes(randomGeneIndex(), randomGeneIndex(), m_rnd->UniformInteger(0,m_trial_chromosome.size()/2));
      case 2:
        m_trial_chromosome.shiftGenes(randomGeneIndex(), m_rnd->UniformInteger(0,m_trial_chromosome.size()), m_rnd->UniformInteger(0,m_trial_chromosome.size()));
      case 3: 
        m_trial_chromosome.reverseGenes(randomGeneIndex(), m_rnd->UniformInteger(2,m_trial_chromosome.size()-2)); //with these parameters it always mutates
    }

    double new_length = calculateLength(m_trial_chromosome);
    double old_length = getLength(); //length of m_chromosome already calculated
    double q = BoltzmannWeight(beta, new_length-old_length);
    if(q >= m_rnd->Rannyu()){
      m_chromosome=m_trial_chromosome;
      m_accepted++;
    } // else use the same chromosome
    m_attempted++;
  }
}



//Functions to generate random towns, not strictly related to annealer, placed here for convenience.

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