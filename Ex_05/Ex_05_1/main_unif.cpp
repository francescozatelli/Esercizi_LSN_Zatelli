#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "metropolis.h"
#include "functions.h"

using namespace std;

int main (int argc, char *argv[]){

	double M_eq=1000; //Number of throws to test equilibrium
	double M = 1E6; //Number of throws
	double L = 10000; //Block size
	double N=M/L; //Number of blocks

	vector<double> r100_throws(M); //single throw
	vector<double> r210_throws(M);
	vector<double> r100_blocks(N); //block values
	vector<double> r210_blocks(N);
	vector<double> r100_avg(N); //progressive average of blocks
	vector<double> r210_avg(N);
	vector<double> r100_avg2(N); //progressive average of blocks^2
	vector<double> r210_avg2(N);
	vector<double> r100_err(N); //progressive standard deviation
	vector<double> r210_err(N);

	Random *rnd=new Random();
	initRandom(*rnd); //Sposto dal main l'inizializzazione del generatore (seed ecc.)
	F100 *f100= new F100(); //1s not normalized probability density
	F210 *f210= new F210(); //2p not normalized probability density

	//Starting point for 100
	double xstart100=0;
	double ystart100=0;
	double zstart100=0;

	//Starting point for 210
	double xstart210=0;
	double ystart210=0;
	double zstart210=2;


	double delta100=1.225;
	double delta210=2.98;

	double alpha100_avg=0;
	double alpha210_avg=0;


	Metropolis p100(xstart100, ystart100, zstart100, f100, rnd);
	Metropolis p210(xstart210, ystart210, zstart210, f210, rnd);

	cout<<"Sampling from hydrogen 100 and 210 probability distribution using Metropolis algorithm"<<endl;
	cout<<"Uniform transition probability"<<endl<<endl;
	cout<<"Number of throws to test equilibrium: "<<M_eq<<endl;
	cout<<"Number of throws: "<<M<<endl;
	cout<<"Number of throws in each block: "<<L<<endl;
	cout<<"Number of blocks: "<<N<<endl;
	cout<<"Uniform distribution step width for 100: "<<delta100<<endl;
	cout<<"Starting point for 100: ("<<xstart100<<", "<<ystart100<<", "<<zstart100<<")"<<endl;
	cout<<"Uniform distribution step width for 210: "<<delta210<<endl;
	cout<<"Starting point for 210: ("<<xstart210<<", "<<ystart210<<", "<<zstart210<<")"<<endl;

	//EQUILIBRATION
	ofstream fileout_eq_points_100;
	ofstream fileout_eq_points_210;
	ofstream fileout_eq_r_100;
	ofstream fileout_eq_r_210;

	fileout_eq_points_100.open("eq100.points");
	if(!fileout_eq_points_100.is_open()){
		cout<<"Unable to open output file eq100.points"<<endl;
		return -1;
	}

	fileout_eq_points_210.open("eq210.points");
	if(!fileout_eq_points_210.is_open()){
		cout<<"Unable to open output file eq100.points"<<endl;
		return -1;
	}

	fileout_eq_r_100.open("out_eq_r100.txt");
	if(!fileout_eq_points_100.is_open()){
		cout<<"Unable to open output file out_eq_r100.txt"<<endl;
		return -1;
	}

	fileout_eq_r_210.open("out_eq_r210.txt");
	if(!fileout_eq_points_210.is_open()){
		cout<<"Unable to open output file out_eq_r210.txt"<<endl;
		return -1;
	}


	fileout_eq_points_100<<"x,\ty,\tz"<<endl;
	fileout_eq_points_210<<"x,\ty,\tz"<<endl;

	for(int i=0; i<M_eq;i++){
		p100.StepUnif(delta100);
		p210.StepUnif(delta210);

		double x100=p100.getX();
		double y100=p100.getY();
		double z100=p100.getZ();
		double x210=p210.getX();
		double y210=p210.getY();
		double z210=p210.getZ();

		double r100=sqrt(x100*x100+y100*y100+z100*z100);
		double r210=sqrt(x210*x210+y210*y210+z210*z210);

		alpha100_avg=i/double(i+1)*alpha100_avg+1./double(i+1)*p100.getAlpha();
		alpha210_avg=i/double(i+1)*alpha210_avg+1./double(i+1)*p210.getAlpha();

		fileout_eq_points_100<<x100<<",\t"<<y100<<",\t"<<z100<<endl;
		fileout_eq_points_210<<x210<<",\t"<<y210<<",\t"<<z210<<endl;
		fileout_eq_r_100<<r100<<endl;
		fileout_eq_r_210<<r210<<endl;
	}
	rnd->SaveSeed();
	fileout_eq_points_100.close();
	fileout_eq_points_210.close();
	fileout_eq_r_100.close();
	fileout_eq_r_210.close();
	cout<<endl<<"========== Finished equilibrium test sampling =============="<<endl;
	cout<<"Output files for points coordinates: eq100.points, eq210.points"<<endl;
	cout<<"Output files for r: out_eq_r100.txt, out_eq_r210.txt"<<endl;

	//SAMPLING AFTER EQUILIBRATION
	ofstream fileout_points_100;
	ofstream fileout_points_210;

	fileout_points_100.open("100.points");
	if(!fileout_points_100.is_open()){
		cout<<"Unable to open output file 100.points"<<endl;
		return -1;
	}

	fileout_points_210.open("210.points");
	if(!fileout_points_210.is_open()){
		cout<<"Unable to open output file 100.points"<<endl;
		return -1;
	}

	fileout_points_100<<"x,\ty,\tz"<<endl;
	fileout_points_210<<"x,\ty,\tz"<<endl;

	for(int i=0; i<M;i++){
		p100.StepUnif(delta100);
		p210.StepUnif(delta210);

		double x100=p100.getX();
		double y100=p100.getY();
		double z100=p100.getZ();
		double x210=p210.getX();
		double y210=p210.getY();
		double z210=p210.getZ();

		double r100=sqrt(x100*x100+y100*y100+z100*z100);
		double r210=sqrt(x210*x210+y210*y210+z210*z210);

		r100_throws[i]=r100;
		r210_throws[i]=r210;

		alpha100_avg=i/double(i+1)*alpha100_avg+1./double(i+1)*p100.getAlpha();
		alpha210_avg=i/double(i+1)*alpha210_avg+1./double(i+1)*p210.getAlpha();

		if(i%1000 == 0){
			fileout_points_100<<x100<<",\t"<<y100<<",\t"<<z100<<endl;
			fileout_points_210<<x210<<",\t"<<y210<<",\t"<<z210<<endl;
		}
	}

	cout<<endl<<"========== Finished sampling =============="<<endl;
	cout<<"Average alpha for 100: "<<alpha100_avg<<endl;
	cout<<"Average alpha for 210: "<<alpha210_avg<<endl;
	rnd->SaveSeed();
	fileout_points_100.close();
	fileout_points_210.close();


	//Calcolo i valori per gli N blocchi (le N "misure")
	for(int i=0; i<N; i++){
		//Fa la media di una "fetta" di vettore che parte da i*L ed è lunga L
		//Per ogni blocco considero una sequenza di L=M/N numeri estratti casualmente
		r100_blocks[i]=averageSlice(r100_throws, i*L, L);
		r210_blocks[i]=averageSlice(r210_throws, i*L, L);
	}

	//Average and error
	for(int i=0; i<N; i++){
		r100_avg[i]=averageSlice(r100_blocks,0,i+1); //Calcola la media dei primi i+1 esperimenti
		r100_avg2[i]=averageQuadSlice(r100_blocks,0,i+1); //Calcola la media dei quadrati dei primi i+1 esperimenti
		r100_err[i]=error(r100_avg,r100_avg2,i); //Calcola la varianza della media per i primi i+1 esperimenti

		r210_avg[i]=averageSlice(r210_blocks,0,i+1); //Calcola la media dei primi i+1 esperimenti
		r210_avg2[i]=averageQuadSlice(r210_blocks,0,i+1); //Calcola la media dei quadrati dei primi i+1 esperimenti
		r210_err[i]=error(r210_avg,r210_avg2,i); //Calcola la varianza della media per i primi i+1 esperimenti
	}



	//Output dei risultati
	ofstream fileout;
	fileout.open("out_r100.txt");
	if(fileout.is_open()){
		fileout<<"N,\tAverage,\tError"<<endl;
		for(int i=0; i<N; i++){
			fileout<<i+1<<",\t"<<r100_avg[i]<<",\t"<<r100_err[i]<<endl;
		}
	}else cerr<< "Unable to open output file"<<endl;

	fileout.close();

	fileout.open("out_r210.txt");
	if(fileout.is_open()){
		fileout<<"N,\tAverage,\tError"<<endl;
		for(int i=0; i<N; i++){
			fileout<<i+1<<",\t"<<r210_avg[i]<<",\t"<<r210_err[i]<<endl;
		}
	}else cerr<< "Unable to open output file"<<endl;
	fileout.close();

	cout<<"Output files for points coordinates: 100.points, 210.points"<<endl;
	cout<<"Output files for <r>: out_r100.txt, out_r210.txt"<<endl;
	return 0;
}
