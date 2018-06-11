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
#ifndef BEMTOOL_FEM_INTERPOLATION_HPP
#define BEMTOOL_FEM_INTERPOLATION_HPP
#include "dof.hpp"

namespace bemtool {

template<typename> class Interpolator;

template<int D>
class Interpolator< BasisFct<P0,D> >{

  typedef BasisFct<P0,D>  ShapeFct;
  typedef array<1,Cplx>   ReturnType;

 private:
  const Mesh<D>&  mesh;

 public:
  Interpolator< BasisFct<P0,D> >(const Mesh<D>& m): mesh(m) {};

  template <typename fct_t>
    ReturnType operator()(const fct_t& f, const int& j){
    ReturnType res; res[0] = f( Ctr(mesh[j]) );
    return res; }

};



template<int D>
class Interpolator< BasisFct<P1,D> >{

  typedef BasisFct<P1,D>  ShapeFct;
  typedef array<D+1,Cplx> ReturnType;

 private:
  const Mesh<D>&  mesh;

 public:
  Interpolator< BasisFct<P1,D> >(const Mesh<D>& m): mesh(m) {};

  template <typename fct_t>
    ReturnType operator()(const fct_t& f, const int& j){
    ReturnType res;
    for(int k=0; k<D+1; k++){
      res[k] = f(mesh[j][k]);}
    return res;}

};



template<>
class Interpolator<RT0_2D>{

  typedef RT0_2D        ShapeFct;
  typedef array<3,Cplx> ReturnType;

 private:
  const Mesh2D&   mesh;
  ShapeFct&       phi;
  
 public:
  Interpolator<RT0_2D>(ShapeFct& psi): phi(psi), mesh(MeshOf(psi)) {}

  template <typename fct_t>
    ReturnType operator()(const fct_t& f, const int& j){
    ReturnType res;
    const Elt2D& e = mesh[j];

    phi.Assign(j);    
    const R3& orientation = OrientationOf(phi);
    
    array<3,R3> n_;
    NormalToFaces(e,n_);
    array<3,Elt1D> edge = EdgesOf(e);
    for(int k=0; k<3; k++){      
      R3   x  = Ctr(edge[k]);
      Real l  = Vol(edge[k]);
      res[k]  = orientation[k]*l*(n_[k],f(x));
    }
    return res;}

};


}


#endif
