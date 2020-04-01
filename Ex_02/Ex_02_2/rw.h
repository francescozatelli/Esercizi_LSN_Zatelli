#ifndef __RandomWalk__
#define __RandomWalk__

#include "random.h"
#include <cmath>


using namespace std;

template <typename T> class RandomWalk{
private:
	T m_x,m_y,m_z;
	Random *m_rnd;
	T m_step;

public:
	RandomWalk(Random *rnd, T step) {
		m_x=0;
		m_y=0;
		m_z=0;
		m_rnd=rnd;
		m_step=step;
	}
	~RandomWalk(){};

	void makeDiscreteStep(){
		//Lancio un dado (estraggo intero tra 1 e 6 con uguale probabilità)
		int dice = ceil(6*m_rnd->Rannyu());
		//Se esce 1 mi sposto lungo x in direzione positiva
		//Se esce 2 mi sposto lungo x in direzione negativa
		//Ecc.
		switch (dice) {
			case 1:
				m_x+=m_step;
				break;
			case 2:
				m_x-=m_step;
				break;
			case 3:
				m_y+=m_step;
				break;
			case 4:
				m_y-=m_step;
				break;
			case 5:
				m_z+=m_step;
				break;
			case 6:
				m_z-=m_step;
				break;
		}
	}

	void makeContinuousStep(){
		// Estrazione gaussiana
		/*
		double x_rnd=m_rnd->Gauss(0,1);
		double y_rnd=m_rnd->Gauss(0,1);
		double z_rnd=m_rnd->Gauss(0,1);

		double norm=sqrt(x_rnd*x_rnd+y_rnd*y_rnd+z_rnd*z_rnd);
		*/
		//Rigetto
		//Estraggo un punto casuale nella "palla" di raggio unitario
		double x_rnd=0;
		double y_rnd=0;
		double z_rnd=0;

		double norm=0;
		do{
			x_rnd=m_rnd->Rannyu(-1,1);
			y_rnd=m_rnd->Rannyu(-1,1);
			z_rnd=m_rnd->Rannyu(-1,1);
			norm=sqrt(x_rnd*x_rnd+y_rnd*y_rnd+z_rnd*z_rnd);
		}while(norm>1);

		

		//Normalizzo a 1 il punto ottenuto (così si trova sulla sfera)
		x_rnd=x_rnd/norm;
		y_rnd=y_rnd/norm;
		z_rnd=z_rnd/norm;

		//Sposta il random walker
		m_x+=m_step*x_rnd;
		m_y+=m_step*y_rnd;
		m_z+=m_step*z_rnd;
	}

	T getX() const{return m_x;}
	T getY() const{return m_y;}
	T getZ() const{return m_z;}

	void setPosition(T x, T y, T z){
		m_x=x;
		m_y=y;
		m_z=z;
	}

	vector<T> getPosition(){
		vector<T> position{m_x,m_y,m_z};
		return position;
	}
	double getDistance(){return sqrt(m_x*m_x+m_y*m_y+m_z*m_z);};
};



#endif // __RandomWalk
