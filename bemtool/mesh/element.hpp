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
#ifndef BEMTOOL_MESH_ELT_HPP
#define BEMTOOL_MESH_ELT_HPP

#include <iostream>
#include <cassert>
#include "../calculus/calculus.hpp"


namespace bemtool{

template<int d, typename ValType = R3>
class Elt{

public:
  typedef ValType         v_t;
  typedef array<d+1,int>  a_t;
  typedef Elt<d,v_t>      e_t;
  typedef Elt<d-1,v_t>    f_t;
  static const int dim =    d;

private:
  // Adresse des sommets
  const v_t* v_[d+1];

  // Adresse du sommet origine
  const v_t* v0;

  // Cle de hachage
  int        key;

public:
  // Constructeurs
  Elt<d,v_t>(): v0(0){};

  Elt<d,v_t>(const e_t& e): v0(e.v0), key(e.key) {
    for(int j=0; j<d+1; j++){v_[j]=(e.v_)[j];}}

  template<typename r_t> Elt<d,v_t>(const r_t& r,
				    const v_t* r0 = 0){
    v0=r0; key=&r[0]-v0;
    for(int j=0; j<d+1; j++){
      v_[j]=&r[j];
      if(key>&r[j]-v0){key=&r[j]-v0;}
    }
  }

  // Affectation
  template<typename r_t> void operator=(const r_t& r){
    key=&r[0]-v0;
    for(int j=0; j<d+1; j++){
      v_[j]=&r[j];
      if(key>&r[j]-v0){key=&r[j]-v0;}
    }
  }

  void operator=(const e_t& e){v0=e.v0; key=e.key;
    for(int j=0; j<d+1; j++){v_[j]=(e.v_)[j];} }

  // Operateur acces
  const v_t& operator[](const int& j) const {return *v_[j%(d+1)];}
  friend const v_t*& Origin(      e_t& e){return e.v0; }
  friend const v_t*  Origin(const e_t& e){return e.v0; }
  friend const int&  Key   (const e_t& e){return e.key;}
  friend const a_t   Num   (const e_t& e){
    if(e.v0==0){std::cout << "error: element origin not intitialised\n";
      exit(EXIT_FAILURE);}
    a_t I; for(int j=0; j<d+1; j++){I[j]=e.v_[j]-e.v0;} return I; }

  // Operateur acces generalise
  template <class i_t>
  subarray<const e_t,i_t> operator[] (const i_t& i_) const {
    return subarray<const e_t,i_t>(*this,i_);}

  // Echange de deux entrees de l'elt
  friend void Swap(e_t& e, const int& j, const int& k){
    const v_t* p=&e[j];   e.v_[j]=e.v_[k]; e.v_[k]=p; }

  // Operateur de flux
  friend std::ostream& operator<<(std::ostream& o, const e_t& e_){
    for(int j=0; j<d+1; j++){o<<e_[j]<<std::endl;} return o;}

  friend bool operator==(const e_t& l_, const e_t& r_){
    bool test=true;
    for(int j=0; j<d+1; j++){
      if( &l_[j] != &r_[j] ){test=false;}}
    return test;}

  friend int Dim(const e_t& e_){return d;}

  friend bemtool::array<d+1,f_t>
    FacesOf(const e_t& e){
    bemtool::array<d+1,f_t> f_;
    bemtool::array<d,int>   k_;
    for(int j=0;j<d;j++)  {k_[j]=j+1;}
    for(int j=0;j<d+1;j++){
      Origin(f_[j]) = Origin(e);
      f_[j]=e[j+k_];}
    return f_;
  }

};


//=======================//
// Declarations de type  //
//=======================//
typedef Elt<0> Elt0D; // noeud
typedef Elt<1> Elt1D; // arrete
typedef Elt<2> Elt2D; // triangle
typedef Elt<3> Elt3D; // tetrahedre


//====================//
// Calcul des aretes  //
//====================//
array<1,Elt1D> EdgesOf(const Elt1D& e){
  array<1,Elt1D> edge;
  Origin(edge[0]) = Origin(e);
  edge[0] = e;
  return edge;}

array<3,Elt1D> EdgesOf(const Elt2D& e){
  array<3,Elt1D> edge; N2 I = N2_(1,2);
  for(int j=0; j<3; j++){
    Origin(edge[j])=Origin(e);
    edge[j]=e[I+j];
  }
  return edge;
}

array<6,Elt1D> EdgesOf(const Elt3D& e){
  array<6,Elt1D> edge; N2 I = N2_(1,2);
  for(int j=0; j<6; j++){Origin(edge[j])=Origin(e);}
  edge[0]=e[I  ];
  edge[1]=e[I+1];
  edge[2]=e[I+2];
  edge[3]=e[I+3];
  I[0]=0;I[1]=2;
  edge[4]=e[I  ];
  edge[5]=e[I+1];
  return edge;
}


//==========================//
// Tri croissant des noeuds //
//==========================//
template <int D>
inline void Order(Elt<D>& e){
  for(int j=1; j<D+1; j++){
    for(int k=0; k<D+1-j; k++){
      if(&e[k]-&e[k+1]>0){Swap(e,k,k+1);}
    }
  }
}


//========//
// Volume //
//========//
inline Real Vol(const Elt1D& e){return norm2(e[1]-e[0]);}
inline Real Vol(const Elt2D& e){return norm2(vprod(e[1]-e[0],e[2]-e[0]))/2.; }
inline Real Vol(const Elt3D& e){return (e[3]-e[0],vprod(e[2]-e[0],e[1]-e[0]))/6.;}


//==========//
// Diametre //
//==========//
inline Real Diam(const Elt1D& e){return norm2(e[1]-e[0]);}
inline Real Diam(const Elt2D& e){
  Real a = norm2(e[1]-e[0]);
  Real b = norm2(e[2]-e[1]); if(a<b){a=b;}
  b = norm2(e[2]-e[0]); if(a<b){a=b;}
  return a;}


//=================//
// Transfo elt ref //
//=================//
template <int D> struct Jac{typedef mat<3,D,Real> Type;};
inline R1x1 MatJac(const Elt<1,R1>& e){return mat_(e[1]-e[0]);}
inline R2x2 MatJac(const Elt<2,R2>& e){return mat_(e[1]-e[0],e[2]-e[0]);}
inline R3x1 MatJac(const Elt1D& e){return mat_(e[1]-e[0]);}
inline R3x2 MatJac(const Elt2D& e){return mat_(e[1]-e[0],e[2]-e[0]);}
inline R3x3 MatJac(const Elt3D& e){return mat_(e[1]-e[0],e[2]-e[0], e[3]-e[0]);}
inline Real DetJac(const Elt1D& e){return norm2(e[1]-e[0]);}
inline Real DetJac(const Elt2D& e){return norm2(vprod(e[1]-e[0],e[2]-e[0]));}
inline Real DetJac(const Elt3D& e){return std::abs((vprod(e[1]-e[0],e[2]-e[0]),e[3]-e[0]));}


//============//
// Barycentre //
//============//
inline R3 Ctr(const Elt1D& e){return (1./2.)*(e[0]+e[1]);}
inline R3 Ctr(const Elt2D& e){return (1./3.)*(e[0]+e[1]+e[2]);}
inline R3 Ctr(const Elt3D& e){return (1./4.)*(e[0]+e[1]+e[2]+e[3]);}


//============================//
// Nombre de noeuds en commun //
//============================//
template <int D>
inline int Adj(const Elt<D>& l, const Elt<D>& r){

  Elt<D> e0=l, e1=r; int p=0;
  for(int j0=0; j0<D+1; j0++){
    for(int j1=p; j1<D+1; j1++){
      if(&e0[j0]==&e1[j1]){
	Swap(e0,p,j0);
	Swap(e1,p,j1);
	p++; break;
      }
    }
  }
  return p;
}


//=============================//
// Comparaison de l'orientation
// de deux elements VOISINS:
// 1 = orientation identique
// 0 = orientation differente
//=============================//
inline bool Comp(const Elt1D& e0, const Elt1D& e1){
  if( e0 == e1 ){ return true;  }
  if( (&e0[0]==&e1[0]) || (&e0[1]==&e1[1]) ){ return false; }
  if( (&e0[0]==&e1[1]) || (&e0[1]==&e1[0]) ){ return true;  }
  std::cout << "\nelt.h: comparaison d'elements non voisins" << std::endl;
  std::abort();
}

inline bool Comp(const Elt2D& e0, const Elt2D& e1){
  if( e0 == e1 ){ return true;  }
  int jj[2],kk[2], n=0;
  for(int j=0; j<3; j++){
    for(int k=0; k<3; k++){
      if( &e0[j]==&e1[k] ){
	jj[n]=j; kk[n]=k; n++;}
    }
  }
  if(n!=2){std::cout << "\nelt.h: comparaison d'elements non voisins" << std::endl; abort();}
  if( (3+jj[1]-jj[0])%3 != (3+kk[1]-kk[0])%3 ){ return true;}
  return false;

}


//=====================================//
// Calcul du vecteur normal a un element
//=====================================//
inline R3 NormalTo(const Elt1D& e){
  R3 e2; e2[2]=1.; return vprod(e2, e[1]-e[0]);}

inline R3 NormalTo(const Elt2D& e){
  return vprod(e[1]-e[0],e[2]-e[1]);}



//===========================//
// Puissance d'un element
//===========================//
inline Real SolidAngle(const R3& p, const Elt1D& e){
  R2x2 M;
  for(int j=0; j<2; j++){
    for(int k=0;k<2;k++){
      M(j,k) = e[k][j] - p[j];
    }
  }
  return det(M)/2;
}

inline Real SolidAngle(const R3& p, const Elt2D& e){
  R3x3 M;
  for(int j=0; j<3; j++){
    for(int k=0;k<3;k++){
      M(j,k) = e[k][j] - p[j];
    }
  }
  return -det(M)/6;
}


//==========================//
//     Calcul normale
//        aux faces
//==========================//

inline void NormalToFaces(const Elt1D& e, array<2,R3>& n_){
  n_[0] = e[1]-e[0];
  normalize(n_[0]);
  n_[1] = (-1.)*n_[0];
}

inline void NormalToFaces(const Elt2D& e, array<3,R3>& n_){

  R3 u1,u2;
  for(int k=0; k<3; k++){

    N2 I;
    I[0]=(k+1)%3;
    I[1]=(k+2)%3;
    Elt1D f = e[I];

    u1 = f[1]-f[0];
    normalize(u1);

    n_[k] = f[0]-e[k];
    n_[k] = n_[k] - (n_[k],u1)*u1;
    normalize(n_[k]);

  }

}

inline void NormalToFaces(const Elt3D& e, array<4,R3>& n_){

  R3 u1,u2;
  for(int k=0; k<4; k++){

    N3 I;
    I[0]=(k+1)%4;
    I[1]=(k+2)%4;
    I[2]=(k+3)%4;
    Elt2D f = e[I];

    u1 = f[1]-f[0];
    normalize(u1);
    u2 = f[2]-f[0];
    u2 = u2 - (u2,u1)*u1;
    normalize(u2);

    n_[k] = f[0]-e[k];
    n_[k] = n_[k] -(n_[k],u1)*u1 -(n_[k],u2)*u2;
    normalize(n_[k]);

  }

}


} // namespace bemtool

#endif
