#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"

using namespace std;

#ifndef __functions__
#define __functions__

void initRandom(Random& rnd);
//Calcola la media di una "fetta" di vettore (fetta parte da first ed è lunga interval)
double averageSlice(const vector<double>& v, int first, int interval);
//Calcola la media dei quadrati di una "fetta" di vettore (fetta parte da first ed è lunga interval)
double averageQuadSlice(const vector<double>& v, int first, int interval);
//Calcola la deviazione standard della media
double error(const vector<double>& avg, const vector<double>& avg2, int n);
//Divide i dati in input in nbins canali. Chiedo in ingresso anche la lunghezza dell'intervallo entro cui sono distribuiti i dati.
vector<int> createHistogram(const vector<double>& data, double interval, int nbins);
//Per il quadrato considero un valore unico per i valori attesi, si può facilmente estendere a considerare un vettore anche per i valori attesi
double calculateChiSquare(const vector<int>& observed, double expected);
//Stampa il vettore contenente gli esiti del lancio del dado
void printDice(const vector<double>& sn, string filename);

#endif // __functions__
