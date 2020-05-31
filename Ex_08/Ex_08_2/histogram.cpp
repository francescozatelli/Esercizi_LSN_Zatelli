#include "histogram.h"
#include<iostream>

Histogram :: Histogram(int nbins, double start, double end){
  m_nbins=nbins;
  m_start=start;
  m_end=end;
  m_interval=end-start;
  m_binlength=m_interval/double(nbins);
  m_histo = vector<int>(m_nbins);
  m_ndata=0;
}

void Histogram :: fill(double point){
  
  int bin_index=floor((point-m_start)/m_binlength);
  if(bin_index<m_nbins and bin_index>=0){
    m_histo[bin_index]++;
    m_ndata++;
  }
}

void Histogram :: fillWithArray(const vector<double>& data){
    for(auto& data_i : data){
      fill(data_i);
    }
}

void Histogram :: reset(){
  for(auto& bin : m_histo){
    bin=0;
  }
  m_ndata=0;
}
