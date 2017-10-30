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
#ifndef BEMTOOL_FEM_SHAPEFCT_HPP
#define BEMTOOL_FEM_SHAPEFCT_HPP

#include "../calculus/calculus.hpp"
#include "../mesh/element.hpp"
#include "../mesh/mesh.hpp"


namespace bemtool {


enum Valuation {ScalarValued, VectorValued};
enum FctType   {P0, P1, P2, RT0};

template <int, int> class BasisFct;
template <typename> class GradBasisFct;
template <typename> class DivBasisFct;

/*============
  FONCTIONS P0
  ============*/

template <int D>
class BasisFct<P0,D>{

public:
  static const int dim             = D;
  static const int nb_dof_loc      = 1;
  typedef array<D,Real>            Rd;
  static const Valuation valuation = ScalarValued;

private:
  const Mesh<D>&   mesh;

public:
  BasisFct<P0,D>(const Mesh<D>& m): mesh(m) {};

  void Assign(const int& j){}

  inline Real operator()(const int& j, const Rd& x = 0.) const {
    return 1.;}

};

typedef BasisFct<P0,1> P0_1D;
typedef BasisFct<P0,2> P0_2D;
typedef BasisFct<P0,3> P0_3D;


/*============
  FONCTIONS P1
  ============*/

template <int D>
class BasisFct<P1,D>{

public:
  static const int dim             = D;
  static const int nb_dof_loc      = D+1;
  typedef array<D,Real>            Rd;
  static const Valuation valuation = ScalarValued;

private:
  const Mesh<D>&   mesh;

public:

  BasisFct<P1,D>(const Mesh<D>& m): mesh(m) {};

  void Assign(const int& j){}

  inline Real operator()(const int& j, const Rd& x) const {
    assert(j>=0 && j<nb_dof_loc);
    if(j==0){return 1.-(Rd(1.),x);} return x[j-1];}

};

typedef BasisFct<P1,1> P1_1D;
typedef BasisFct<P1,2> P1_2D;
typedef BasisFct<P1,3> P1_3D;


/*======================
  GRADIENTS FONCTIONS P1
  ======================*/

template <int D>
class GradBasisFct< BasisFct<P1,D> >{

public:
  static const int dim             = D;
  static const int nb_dof_loc      = D+1;
  typedef array<D,Real>            Rd;
  typedef array<D+1,R3>            ArrayGrad;
  static const Valuation valuation = VectorValued;

private:
  const Mesh<D>&      mesh;
  std::vector<ArrayGrad>   grad;
  int                 num_elt;

public:

  GradBasisFct< BasisFct<P1,D> >(const Mesh<D>& m): mesh(m) {
    int nb_elt = NbElt(mesh);
    grad.resize(nb_elt);
    ArrayGrad n;

    for(int j=0; j<nb_elt; j++){
      const Elt<D>& e = mesh[j];
      NormalToFaces(e,n);
      for(int k=0; k<D+1; k++){
	Real z = -1./( n[k],(e[k+1]-e[k]) );
	grad[j][k] = z*n[k];
      }
    }
  }

  void Assign(const int& j){num_elt = j;}

  const R3& operator()(const int& j, const Rd& x = 0.) const {
    assert(j>=0 && j<nb_dof_loc);
    return grad[num_elt][j];}

};

typedef GradBasisFct<P1_1D> Grad_P1_1D;
typedef GradBasisFct<P1_2D> Grad_P1_2D;
typedef GradBasisFct<P1_3D> Grad_P1_3D;



/*============
  FONCTIONS P2
  ============*/

template<>
class BasisFct<P2,1>{

public:
  static const int dim         = 1;
  static const int nb_dof_loc  = 3;
  typedef array<1,Real>        Rd;
  static const Valuation valuation = ScalarValued;

private:
  const Mesh1D&     mesh;

public:
  BasisFct<P2,1>(const Mesh1D& m): mesh(m) {};

  void Assign(const int& j){}

  Real lambda(const int& j, const Rd& x) const {
    if(j==0){return 1.-x[0];} return x[0];}

  Real operator()(const int& j, const Rd& x) const {
    assert(j>=0 && j<nb_dof_loc);
    if(j<2){Real lj = lambda(j,x); return lj*(2.*lj-1.);}
    return 4.*x[0]*(1.-x[0]);}
};


template<>
class BasisFct<P2,2>{

public:
  static const int dim         = 2;
  static const int nb_dof_loc  = 6;
  typedef array<2,Real>        Rd;
  static const Valuation valuation = ScalarValued;

private:
  const Mesh2D&     mesh;

public:
  BasisFct<P2,2>(const Mesh2D& m): mesh(m) {};

  void Assign(const int& j){}

  Real lambda(const int& j, const Rd& x) const {
    if(j==0){return 1.-(Rd(1.),x);} return x[j-1];}

  Real operator()(const int& j, const Rd& x) const {
    assert(j>=0 && j<nb_dof_loc);
    if(j<3){Real lj = lambda(j,x); return lj*(2.*lj-1.);}
    return 4.*lambda((j-1)%3,x)*lambda((j-2)%3,x);}
};


typedef BasisFct<P2,1> P2_1D;
typedef BasisFct<P2,2> P2_2D;




/*======================
  GRADIENTS FONCTIONS P2
  ======================*/

template <>
class GradBasisFct<P2_1D>{

public:
  static const int dim             = 1;
  static const int nb_dof_loc      = 3;
  typedef array<1,Real>            Rd;
  typedef array<2,R3>              ArrayGrad;
  static const Valuation valuation = VectorValued;

private:
  const Mesh1D&       mesh;
  std::vector<ArrayGrad>   grad;
  int                 num_elt;

public:

  GradBasisFct<P2_1D>(const Mesh1D& m): mesh(m) {
    int nb_elt = NbElt(mesh);
    grad.resize(nb_elt);
    ArrayGrad n;

    for(int j=0; j<nb_elt; j++){
      const Elt1D& e = mesh[j];
      NormalToFaces(e,n);
      for(int k=0; k<2; k++){
	Real z = -1./( n[k],(e[k+1]-e[k]) );
	grad[j][k] = z*n[k];
      }
    }
  }

  void Assign(const int& j){num_elt = j;}

  Real lambda(const int& j, const Rd& x) const {
    if(j==0){return 1.-x[0];} return x[0];}

  R3 operator()(const int& j, const Rd& x) const {
    assert(j>=0 && j<nb_dof_loc);
    if(j<2){Real lj = lambda(j,x); return (4.*lj-1.)*grad[num_elt][j];}
    return 4.*(1.-2.*x[0])*grad[num_elt][1];}

};


template <>
class GradBasisFct<P2_2D>{

public:
  static const int dim             = 2;
  static const int nb_dof_loc      = 6;
  typedef array<2,Real>            Rd;
  typedef array<3,R3>              ArrayGrad;
  static const Valuation valuation = VectorValued;

private:
  const Mesh2D&      mesh;
  std::vector<ArrayGrad>   grad;
  int                 num_elt;

public:

  GradBasisFct<P2_2D>(const Mesh2D& m): mesh(m) {
    int nb_elt = NbElt(mesh);
    grad.resize(nb_elt);
    ArrayGrad n;

    for(int j=0; j<nb_elt; j++){
      const Elt2D& e = mesh[j];
      NormalToFaces(e,n);
      for(int k=0; k<3; k++){
	Real z = -1./( n[k],(e[k+1]-e[k]) );
	grad[j][k] = z*n[k];
      }
    }
  }

  void Assign(const int& j){num_elt = j;}

  Real lambda(const int& j, const Rd& x) const {
    if(j==0){return 1.-(Rd(1.),x);} return x[j-1];
  }

  R3 operator()(const int& j, const Rd& x = 0.) const {
    assert(j>=0 && j<nb_dof_loc);
    if(j<3){Real lj = lambda(j,x); return (4.*lj-1.)*grad[num_elt][j];}
    R3 u = 4.*lambda((j-1)%3,x)*grad[num_elt][(j-2)%3];
    u += 4.*lambda((j-2)%3,x)*grad[num_elt][(j-1)%3];
    return u;
  }


};


typedef GradBasisFct<P2_1D> Grad_P2_1D;
typedef GradBasisFct<P2_2D> Grad_P2_2D;





/*=============================
  RAVIART-THOMAS RT0
  =============================*/

template <>
class BasisFct<RT0,2>{

public:
  static const Valuation valuation = VectorValued;
  static const int dim             = 2;
  static const int nb_dof_loc      = 3;
  typedef array<2,Real>            Rd;
  typedef BasisFct<RT0,2>          this_t;

private:
  const Mesh2D&    mesh;
  std::vector<R3x2>     mat_jac;
  int              num_elt;
  R2               a[3];
  std::vector<R3>       orientation;

public:
  BasisFct<RT0,2>(const Mesh2D& m): mesh(m) {

    int nb_elt = NbElt(mesh);
    orientation.resize(nb_elt);
    a[1][0]=1.; a[2][1]=1.;
    mat_jac.resize(NbElt(mesh));

    const int nb_edge_loc = 3;
    int nb_edge =  0;
    int end     = -1;
    int nb_node = NbNode(GeometryOf(mesh));
    std::vector<int>   begin(nb_node,end);
    std::vector<int>   next;
    std::vector<Elt1D> edge;

    for(int j=0; j<nb_elt; j++){
      array<nb_edge_loc,Elt1D> edge_loc = EdgesOf(mesh[j]);
      for(int k=0; k<nb_edge_loc; k++){
	bool exists=false;
	Elt1D e = edge_loc[k]; Order(e);
	for(int p=begin[Key(e)]; p!=end; p=next[p]){
	  if(e==edge[p]){
	    exists=true;
	    orientation[j][k]=-1.;
	    break;
	  }
	}
	if(!exists){
	  next.push_back(begin[Key(e)]);
	  begin[Key(e)] = nb_edge;
	  edge.push_back(e);
	  orientation[j][k] = +1.;
	  nb_edge++;
	}
      }
    }

    for(int j=0; j<nb_elt; j++){
      const Elt2D& e = mesh[j];
      R3x2 B = mat_(e[1]-e[0],e[2]-e[0]);
      mat_jac[j] = (0.5/Vol(e))*B;
    }

  }

  void Assign(const int& j){num_elt = j;}

  R3 operator()(const int& j, const Rd& x) {
    assert(j>=0 && j<nb_dof_loc);
    R3 u = mat_jac[num_elt]*(x-a[j]);
    return orientation[num_elt][j]*u;
  }

};


typedef BasisFct<RT0,2> RT0_2D;


/*==============
  DIVERGENCE RT0
  ==============*/

template <>
class DivBasisFct<RT0_2D>{

public:
  static const int dim         = 2;
  static const int nb_dof_loc  = 3;
  typedef array<2,Real>        Rd;
  static const Valuation valuation = ScalarValued;

private:
  const Mesh2D&    mesh;
  std::vector<R3>       orientation;
  std::vector<Real>     volume;
  int              num_elt;

public:
  DivBasisFct<RT0_2D>(const Mesh2D& m): mesh(m){

    int nb_elt = NbElt(mesh);
    orientation.resize(nb_elt);
    volume.resize(nb_elt);

    const int nb_edge_loc = 3;
    int nb_edge =  0;
    int end     = -1;
    int nb_node = NbNode(GeometryOf(mesh));
    std::vector<int>   begin(nb_node,end);
    std::vector<int>   next;
    std::vector<Elt1D> edge;

    for(int j=0; j<nb_elt; j++){
      volume[j] = Vol(mesh[j]);
      array<nb_edge_loc,Elt1D> edge_loc = EdgesOf(mesh[j]);
      for(int k=0; k<nb_edge_loc; k++){
	bool exists=false;
	Elt1D e = edge_loc[k]; Order(e);
	for(int p=begin[Key(e)]; p!=end; p=next[p]){
	  if(e==edge[p]){
	    exists=true;
	    orientation[j][k]=-1.;
	    break;
	  }
	}
	if(!exists){
	  next.push_back(begin[Key(e)]);
	  begin[Key(e)] = nb_edge;
	  edge.push_back(e);
	  orientation[j][k] = +1.;
	  nb_edge++;
	}
      }
    }
  }

  void Assign(const int& j){num_elt = j;}

  Real operator()(const int& j, const Rd& x = 0.) const {
    assert(j>=0 && j<nb_dof_loc);
    return orientation[num_elt][j]/volume[num_elt];}

};

typedef DivBasisFct<RT0_2D> Div_RT0_2D;

} // namespace bemtool


#endif
