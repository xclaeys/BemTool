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
#ifndef BEMTOOL_FEM_FEMP1_HPP
#define BEMTOOL_FEM_FEMP1_HPP

#include "../mesh/element.hpp"

namespace bemtool{

  inline R2x2 StiffP1(const Elt1D& e){
    Real l = Vol(e);
    R2x2 K;
    K(0,0) =  1./l;   K(1,1) =  1./l;
    K(0,1) = -1./l;   K(1,0) = -1./l;
    return K;
  }

  inline R3x3 StiffP1(const Elt2D& e){
    R3 U[3];
    U[0] = e[2]-e[1];
    U[1] = e[0]-e[2];
    U[2] = e[1]-e[0];
    Real area = Vol(e);

    R3x3 K;
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
	K(j,k) = (U[j],U[k])/(4.*area);
      }
    }
    return K;
  }


  inline R2x2 MassP1(const Elt1D& e){
    Real l = Vol(e);
    R2x2 M;
    M(0,0) = l/3.; M(1,1) = l/3.;
    M(0,1) = l/6.; M(1,0) = l/6.;
    return M;
  }


  inline R3x3 MassP1(const Elt2D& e){
    Real a = Vol(e);
    R3x3 M;
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
	if(j==k){ M(j,k) = a/6. ; }
	else    { M(j,k) = a/12.; }
      }
    }
    return M;
  }


  inline R4x4 MassP1(const Elt3D& e){
    Real v = Vol(e);
    R4x4 M;
    for(int j=0; j<4; j++){
      for(int k=0; k<4; k++){
	if(j==k){ M(j,k) = v/10.; }
	else    { M(j,k) = v/20.; }
      }
    }
    return M;
  }




}
#endif
