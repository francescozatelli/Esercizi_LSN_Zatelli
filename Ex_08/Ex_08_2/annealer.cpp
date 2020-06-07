#include "annealer.h"


//Returns the boltzmann weight for a given energy
double Annealer :: BoltzmannWeight(double beta, double energy){
    return exp(-beta*energy);
}

double Annealer :: CalculateEnergy(double mu, double sigma){
  ProbabilityDistribution p(mu,sigma); //Initialize functions of interest with the new values for mu and sigma
  LocalEnergy eloc(mu,sigma);

  double sum_num=0;
  double sum_den=0;
  for(int i=0; i<m_metrosteps;i++){
    m_energies_new[i]=eloc.Eval(m_points[i]);
    m_weights_new[i]=p.Eval(m_points[i])/m_pevals[i];
    sum_den+=m_weights_new[i];
    sum_num+=m_energies_new[i]*m_weights_new[i];
	}
  return sum_num/sum_den;
}


Annealer :: Annealer(double mu0, double sigma0, double metrowidth, int metrosteps, Random *rnd){
  m_rnd=rnd;
  m_mu=mu0;
  m_sigma=sigma0;
  m_metrowidth=metrowidth;
  m_metrosteps=metrosteps;

  m_energies=vector<double>(metrosteps);
  m_energies_new=vector<double>(metrosteps);
  m_points=vector<double>(metrosteps);
  m_pevals=vector<double>(metrosteps);
  m_weights=vector<double>(metrosteps);
  m_weights_new=vector<double>(metrosteps);

  m_accepted=0;
  m_attempted=0;

  samplePath(mu0,sigma0);
  m_energy=CalculateEnergy(mu0, sigma0);
  cout<<"Starting energy: "<<m_energy<<endl;
}


void Annealer :: annealingStep(double beta, int nsteps, double mu_width, double sigma_width, double metrowidth){
  m_metrowidth=metrowidth; //step width for energy calculation
  m_accepted=0;
  m_attempted=0;
  samplePath(m_mu, m_sigma);
  for(int i=0; i<nsteps; i++){
    double mu_new = m_rnd->Rannyu(m_mu-mu_width, m_mu+mu_width);
    double sigma_new = m_rnd->Rannyu(m_sigma-sigma_width, m_sigma+sigma_width);
    double energy_new = CalculateEnergy(mu_new, sigma_new);
    double q = BoltzmannWeight(beta, energy_new)/BoltzmannWeight(beta, m_energy);
      if(q >= m_rnd->Rannyu()){
        m_energy=energy_new;
        m_mu=mu_new;
        m_sigma=sigma_new;
        m_weights=m_weights_new;
        m_energies=m_energies_new;
        m_accepted++;
      } // else use the same point
      m_attempted++;
  }
}

void Annealer :: samplePath(double mu, double sigma){
  m_psample.setMu(mu); //Probability distribution used to sample the path
  m_psample.setSigma(sigma);

  LogarithmicDerivativeP dlogp(mu,sigma);
  LocalEnergy eloc(mu,sigma);
  Metropolis metro(mu, &m_psample, &dlogp, m_rnd); //mu as starting point

  for(int i=0; i<m_metrosteps;i++){
    metro.StepSmart(m_metrowidth);
		//metro.StepUnif(m_metrowidth);
    m_points[i]=metro.getX();
    m_pevals[i]=m_psample.Eval(m_points[i]); //save also the values of p(x) in order not to evaluate it again in every cycle to calculate the weights
    m_energies[i]=eloc.Eval(m_points[i]);
    m_weights[i]=1;
  }

    //cout<<"Sample acceptance rate: "<<metro.getAcceptanceRate()<<endl;
}

tuple<double,double> Annealer :: dataBlockingEnergy(int block_size){
  int N=m_metrosteps/block_size;
  vector<double> energy_blocks(N);

  //Calcolo i valori per gli N blocchi (le N "misure")
  double sum_num=0;
  double sum_den=0;
	for(int i=0; i<N; i++){
    double block_avg=0;
    for(int j=0; j<block_size; j++){
      sum_num+=m_energies[j+i*block_size]*m_weights[j+i*block_size];
      sum_den+=m_weights[j+i*block_size];
    }
		energy_blocks[i]=sum_num/sum_den;
  }

  double energy_avg=0;
  double energy_avg2=0;
  double energy_err=0;

  for(int i =0; i<N; i++){
    energy_avg = i/double(i+1)*energy_avg+1./double(i+1)*energy_blocks[i];
    energy_avg2 = i/double(i+1)*energy_avg2+1./double(i+1)*energy_blocks[i]*energy_blocks[i];
  }
  energy_err=sqrt((energy_avg2-energy_avg*energy_avg)/double(N-1));

  return make_tuple(energy_avg, energy_err);
}