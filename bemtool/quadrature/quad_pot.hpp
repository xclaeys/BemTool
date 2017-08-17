#ifndef QUADPOT_H
#define QUADPOT_H

#include "quad.hpp"
#include <fstream>

template <int> class QuadPot;

template<> class QuadPot<1>{
  
public:
  typedef Real qp_t;
  typedef array<1,Real> R1;

private:
  vector<R1>   x;
  vector<Real> w;

public:
  const vector<R1>&   GetPoints () const {return x;}
  const vector<Real>& GetWeights() const {return w;}  
  
  QuadPot<1>(const int& order){    
    vector<Real> w0, x0;
    int n = (int) ceil((order+1)/2.0);
    Quad1D(n,x0,w0);
    for(int j=0; j<x0.size(); j++){
      w.push_back(w0[j]);
      x.push_back(x0[j]);
    }
  }
  
};





template<> class QuadPot<2>{
  
public:
  typedef R2   qp_t;
  
private:
  vector<R2>   x;
  vector<Real> w;

public:
  const vector<R2>&   GetPoints () const {return x;}
  const vector<Real>& GetWeights() const {return w;}  

  QuadPot<2>(const int& order){      
    Quad2D(order,x,w);}

};





#endif
