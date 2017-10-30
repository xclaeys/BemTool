//===================================================================
//
//  Copyright 2017 Xavier Claeys
//
//  bemtool is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  bemtool is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with bemtool.  If not, see <http://www.gnu.org/licenses/>
//
//====================================================================
#ifndef BEMTOOL_QUADRATURE_QUADPOT_HPP
#define BEMTOOL_QUADRATURE_QUADPOT_HPP

#include "quad.hpp"
#include <fstream>


namespace bemtool {

template <int> class QuadPot;

template<> class QuadPot<1>{

public:
  typedef Real qp_t;
  typedef array<1,Real> R1;

private:
  std::vector<R1>   x;
  std::vector<Real> w;

public:
  const std::vector<R1>&   GetPoints () const {return x;}
  const std::vector<Real>& GetWeights() const {return w;}

  QuadPot<1>(const int& order){
    std::vector<Real> w0, x0;
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
  std::vector<R2>   x;
  std::vector<Real> w;

public:
  const std::vector<R2>&   GetPoints () const {return x;}
  const std::vector<Real>& GetWeights() const {return w;}

  QuadPot<2>(const int& order){
    Quad2D(order,x,w);}

};

}



#endif
