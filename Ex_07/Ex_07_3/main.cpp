#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "molecular_dynamics.h"

using namespace std;

int main(){

  MolecularDynamics simulation;

  int nconf=1;
  for(int istep=1; istep<=simulation.getSteps(); istep++){
    if(istep%1000 == 0) cout <<"Number of time steps: "<<istep<<endl;
    if(istep%20 == 0){
      simulation.Measure();
      //simulation.ConfXYZ(nconf);
      nconf+=1;
    }
    if(istep%2000 == 0){
      simulation.MeasureBlock();
    }
    simulation.Move();
  }
  simulation.ConfFinal();
  return 0;
}
