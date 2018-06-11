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
#ifndef BEMTOOL_POTENTIAL_POTENTIAL_HPP
#define BEMTOOL_POTENTIAL_POTENTIAL_HPP

#include "../quadrature/quad_pot.hpp"
#include "../fem/shapefct.hpp"
#include "../equations.hpp"
#include "../miscellaneous/specialfct.hpp"


namespace bemtool {


/*=====================
  Type de potentiels
  =====================*/

enum PotKernelEnum {
  CST_POT, // noyau constant
   SL_POT, // simple couche
   DL_POT  // double couche
};

template <typename KernelType> class Potential{

public:

  static  const int dimy = KernelType::Trait::dimy;
  typedef typename KernelType::Trait          KernelTypeTrait;
  typedef typename KernelTypeTrait::MeshY     MeshY;
  typedef typename KernelTypeTrait::Rdy       RdY;
  typedef typename KernelTypeTrait::EltY      EltY;
  typedef typename KernelTypeTrait::MatType   MatType;
  typedef QuadPot<dimy>                       QuadType;

private:

  const MeshY&     meshy;
  const Geometry&  nodey;
  KernelType       ker;
  QuadType         qr;
  MatType          mat;
  Cplx             val,val2;

public:

  Potential<KernelType>(const MeshY& my, const Real& k):
  ker(my,k), meshy(my), nodey(GeometryOf(my)), qr(10) {};


  const MatType& operator()(const R3& x, const int& jy){
    const std::vector<RdY>&  t  = qr.GetPoints();
    const std::vector<Real>& w  = qr.GetWeights();
    mat=0; ker.Assign(x,jy);
    for(int j=0; j<w.size(); j++){
      mat += w[j]*ker(x,t[j]);}
    return mat;
  }


  const Cplx& operator()(const R3& x, const N2& Iy){
    const std::vector<RdY>&  t  = qr.GetPoints();
    const std::vector<Real>& w  = qr.GetWeights();
    val=0; ker.Assign(x,Iy[0]);
    for(int j=0; j<w.size(); j++){
      val += w[j]*ker(x,t[j],Iy[1]);}
    return val;
  }

  const Cplx& operator()(const R3& x, const std::vector<N2>& vy){
    val2=0.;
    for(int iy=0; iy<vy.size(); iy++){
      val2 += (*this)(x,vy[iy]);}
    return val2;
  }

};


/*==============================
  GENERALITES SUR LES NOYAUX
  =============================*/

template <int,int,int,typename>
class PotKernel;

  template <typename PhiY, int Valuation=1>
class PotKernelTraits{

public:
  static const int dimy     = PhiY::dim;
  static const int nb_dof_y = PhiY::nb_dof_loc;

  typedef typename Jac<dimy>::Type      JacY;
  typedef typename PhiY::Rd             Rdy;
  typedef Mesh<dimy>                    MeshY;
  typedef Elt<dimy>                     EltY;
  typedef GradBasisFct<PhiY>            GradPhiY;
  typedef DivBasisFct<PhiY>             DivPhiY;
  typedef mat<Valuation,nb_dof_y,Cplx>  MatType;

};


/*=======================
        NOYAU CONSTANT
  =======================*/

template <int D, typename PhiY>
class PotKernel<CST,CST_POT,D,PhiY>{

public:
  typedef PotKernelTraits<PhiY> Trait;

private:
  const typename Trait::MeshY&    meshy;
        typename Trait::MatType   mat;
                        PhiY      phiy;
                        Real      ker;
                        Cplx      val;

public:

  PotKernel<CST,CST_POT,D,PhiY>(const typename Trait::MeshY& my,
				const Real& k=0.): meshy(my), phiy(my) {};


  void Assign(const R3& x, const int& iy){
    const typename Trait::EltY& ey=meshy[iy];
    ker = DetJac(ey);}


  const typename Trait::MatType&
  operator()(const R3& x,const typename Trait::Rdy& tj){
    for(int k=0; k<Trait::nb_dof_y; k++){
      mat(0,k) = ker*phiy(k,tj);}
    return mat;}


  const Cplx&
  operator()(const R3& x,const typename Trait::Rdy& tj, const int& ky){
    return val = ker*phiy(ky,tj);}

};

typedef PotKernel<CST,CST_POT,1,P0_1D> CST_1D_P0;
typedef PotKernel<CST,CST_POT,1,P1_1D> CST_1D_P1;
typedef PotKernel<CST,CST_POT,2,P0_2D> CST_2D_P0;
typedef PotKernel<CST,CST_POT,2,P1_2D> CST_2D_P1;

}
#endif
