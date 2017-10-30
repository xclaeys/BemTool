//===================================================================
//
//  Copyright 2017 Xavier Claeys
//
//  bemtool is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  bemtool is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//  (see ../LICENSE.txt)
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with bemtool.  If not, see <http://www.gnu.org/licenses/>
//
//====================================================================
#ifndef BEMTOOL_MISC_COORDINATES_HPP
#define BEMTOOL_MISC_COORDINATES_HPP

#include <cmath>
#include "../calculus/calculus.hpp"

namespace bemtool{

Real Atan2(const Real& y, const Real& x){
  Real theta = atan2(y,x);
  if(theta>=0){return theta;}
  return 2*pi+theta;
}


/*========================
  COORDONNEES CYLINDRIQUES
  ========================*/

class Cyl{

private:

  R3   x;
  Real r;
  Real theta;
  Real z;

  void ComputeCoordinates(){
    r = sqrt(x[0]*x[0]+x[1]*x[1]);
    z = x[2]; theta = Atan2(x[1],x[0]);
  }

public:

  Cyl(const R3&  p=0.): x(p) {
    r = sqrt(x[0]*x[0]+x[1]*x[1]);
    z = x[2]; theta = Atan2(x[1],x[0]);
  }

  Cyl(const Cyl& p): x(p.x), r(p.r), theta(p.theta), z(p.z) {};

  template <typename r_t> void operator=(const r_t& r){
    x=r;ComputeCoordinates();}

  const Real& operator[] (const int& j) const {return x[j];}
  const Real& R()     const {return r;    }
  const Real& Theta() const {return theta;}
  const Real& Z()     const {return z;    }

  static Real R(const R3& p)    {return sqrt(p[0]*p[0]+p[1]*p[1]);}
  static Real Theta(const R3& p){return Atan2(p[1],p[0]);}
  static Real Z(const R3& p)    {return p[2];}

};




/*======================
  COORDONNEES SPHERIQUES
  ======================*/

class Sph{

private:

  R3   x;
  Real rho;
  Real theta;
  Real phi;

  void ComputeCoordinates(){
    rho   = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    theta = acos(x[2]/rho);
    phi   = Atan2(x[1],x[0]);
  }

public:

  Sph(const R3&  p=0.): x(p) {
    rho   = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    theta = acos(x[2]/rho);
    if(sin(theta)<1e-8){phi=0.;}
    else{phi = Atan2(x[1],x[0]);}
  }

  Sph(const Sph& p): x(p.x),rho(p.rho),theta(p.theta),phi(p.phi) {};

  template <typename r_t> void operator=(const r_t& r){
    x=r; ComputeCoordinates();}

  const Real& operator[]   (const int& j) const {return x[j];}
  const Real& Rho()   const {return rho;  }
  const Real& Theta() const {return theta;}
  const Real& Phi()   const {return phi;  }

  static Real Rho(const R3& p){
    return sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);}
  static Real Theta(const R3& p){
    return acos(p[2]/Sph::Rho(p));}
  static Real Phi(const R3& p){
    if(p[0]*p[0]+p[1]*p[1]<1e-8){return 0.;}
    return Atan2(p[1],p[0]); }

};


}



#endif
