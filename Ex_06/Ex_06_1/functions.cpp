#include "functions.h"

void initRandom(Random& rnd){
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
	   Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
	   while ( !input.eof() ){
		  input >> property;
		  if( property == "RANDOMSEED" ){
			 input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
			 rnd.SetRandom(seed,p1,p2);
		  }
	   }
	   input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
}

double averageSlice(const vector<double>& v, int first, int interval){
	double sum=0;
	for(int i=first; i<first+interval; i++){
		sum+=v[i];
	}
	return sum/interval;
}

double averageQuadSlice(const vector<double>& v, int first, int interval){
	double sum=0;
	for(int i=first; i<first+interval; i++){
		sum+=v[i]*v[i];
	}
	return sum/interval;
}

double error(const vector<double>& avg, const vector<double>& avg2, int n){
	if(n==0){
		return 0;
	}else{
		return sqrt((avg2[n]-avg[n]*avg[n])/n);
	}
}

vector<int> createHistogram(const vector<double>& data, double interval, int nbins){
	double bin_length=interval/nbins;
	vector<int> histo(nbins);

	for(auto& data_i : data){
		int bin_index=floor(data_i/bin_length);
		histo[bin_index]++;
	}

	return histo;
}


double calculateChiSquare(const vector<int>& observed, double expected){
	double chi=0;
	for(auto& obs_i : observed){
		chi+=pow(obs_i-expected,2)/expected;
	}
	return chi;
}
void printVector(const vector<double>& v, string filename){
	ofstream fileout;
	fileout.open(filename);
	if(fileout.is_open()){
		for(auto & v_i : v){
			fileout<<v_i<<endl;
		}
	}else cerr<< "Unable to open output file"<<endl;

	fileout.close();
}
