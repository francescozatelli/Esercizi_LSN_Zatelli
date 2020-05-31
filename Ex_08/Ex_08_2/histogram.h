#ifndef __Histogram__
#define __Histogram__

#include<vector>
#include<cmath>
#include<iostream>

using namespace std;

class Histogram {

private:
  vector<int> m_histo;
  int m_nbins;
  double m_binlength;
  double m_start;
  double m_end;
  double m_interval;
  int m_ndata;

public:
  Histogram(int nbins, double start, double end);
	~Histogram(){};

  void fill(double point);
  void fillWithArray(const vector<double>& data);
  int getBin(int nbin) const{ return m_histo[nbin];};
  double getBinLength() const{ return m_binlength;};
  int getNumberOfBins() const{ return m_nbins;};
  int getNumberOfData() const{ return m_ndata;};
  void reset();
};

#endif // __Histogram__
