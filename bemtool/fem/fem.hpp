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
#ifndef BEMTOOL_FEM_FEM_HPP
#define BEMTOOL_FEM_FEM_HPP
#include "../quadrature/quad.hpp"
#include "../mesh/mesh.hpp"


namespace bemtool {

template <typename T>
Real ScalarProd(const T&, const T&);


template <>
Real ScalarProd(const Real& u, const Real& v){
  return u*v;}

template <int D>
Real ScalarProd(const array<D,Real>& u, const array<D,Real>& v){
  return (u,v);}

template <typename PhiX,typename PhiY = PhiX>
class LocalMatrix{

  static const int dim      = PhiX::dim;
  static const int nb_dof_x = PhiX::nb_dof_loc;
  static const int nb_dof_y = PhiY::nb_dof_loc;
  typedef mat<nb_dof_x,nb_dof_y,Real>  ReturnType;
  typedef Mesh<dim>                    MeshType;
  typedef typename PhiX::Rd            Rd;

private:
  const MeshType&     mesh;
  ReturnType          inter;
  PhiX                phix;
  PhiY                phiy;
  Real                h;
  QuadFEM<dim>        qr;
  const std::vector<Rd>&   t;
  const std::vector<Real>& dt;

public:
  LocalMatrix<PhiX,PhiY>(const MeshType& m):
  mesh(m), phix(m), phiy(m),
    qr(5), t(points(qr)), dt(weights(qr)) {};

  const ReturnType& operator()(const int& l) {
    const Elt<dim>& e = mesh[l];
    phix.Assign(l);
    phiy.Assign(l);
    inter = 0.;
    h     = DetJac(e);
    for(int j=0; j<t.size(); j++){
      for(int kx = 0; kx<nb_dof_x; kx++){
	for(int ky = 0; ky<nb_dof_y; ky++){
	  inter(kx,ky)+= ScalarProd( phix(kx,t[j]),phiy(ky,t[j]) )*h*dt[j];
	}
      }
    }
    return inter;
  }

};

}



#endif
