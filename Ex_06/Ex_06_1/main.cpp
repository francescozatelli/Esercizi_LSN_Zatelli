#include <iostream>
#include "ising1D.h"
#include "random.h"
#include "functions.h"

using namespace std;

int main(){

  Random *rnd=new Random();
	initRandom(*rnd); //Sposto dal main l'inizializzazione del generatore (seed ecc.)


  Ising1D simulation(rnd);

  for(int i=0; i<simulation.getNumberOfBlocks(); i++){
    for(int j = 0; j<simulation.getNumberOfStepsPerBlock(); j++){
      simulation.Move();
      simulation.Measure();
    }
    simulation.MeasureBlock();
  }
  cout<<"Average acceptance rate: "<<simulation.getAlpha()<<endl;
  simulation.ConfFinal();
  rnd->SaveSeed();
  
  return 0;
}
