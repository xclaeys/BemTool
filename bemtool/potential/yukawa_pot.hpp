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
#ifndef BEMTOOL_POTENTIAL_YUKAWA_HPP
#define BEMTOOL_POTENTIAL_YUKAWA_HPP

#include "potential.hpp"


namespace bemtool {

/*==========================
  SIMPLE COUCHE YUKAWA EN 2D
  ==========================*/

template <typename PhiY>
class PotKernel<YU,SL_POT,2,PhiY>{

public:
  typedef PotKernelTraits<PhiY> Trait;

private:
  const typename Trait::MeshY&   meshy;
        typename Trait::MatType  mat;
        typename Trait::JacY     dy;
                        PhiY     phiy;
  const                 Real     kappa;
                        R3       x0_y0,x_y;
                        Real     h,r;
                        Cplx     ker;

public:

  PotKernel<YU,SL_POT,2,PhiY>(const typename Trait::MeshY& my,
			      const Real& k):
  meshy(my), phiy(my), kappa(k) {};

  void Assign(const R3& x, const int& iy){
    const typename Trait::EltY& ey=meshy[iy];
    x0_y0 = x-ey[0];
    dy    = MatJac(ey);
    h     = DetJac(ey);}

  const typename Trait::MatType&
  operator()(const R3& x,const typename Trait::Rdy& tj){
    x_y  = x0_y0-dy*tj;
    r   = norm2(x_y);
    ker = h*(0.5/pi)*Modified_BesselK0(kappa*r);
    for(int k=0; k<Trait::nb_dof_y; k++){
      mat(0,k) = ker*phiy(k,tj);}
    return mat;}

  const Cplx&
  operator()(const R3& x, const typename Trait::Rdy& tj, const int& ky){
    x_y  = x0_y0-dy*tj;
    r   = norm2(x_y);
    ker = h*(0.5/pi)*Modified_BesselK0(kappa*r);
    return ker *= phiy(ky,tj);
  }

};

  typedef PotKernel<YU,SL_POT,2,P0_1D> YU_SL_2D_P0;
  typedef PotKernel<YU,SL_POT,2,P1_1D> YU_SL_2D_P1;



/*==========================
  SIMPLE COUCHE YUKAWA EN 3D
  ==========================*/

template <typename PhiY>
class PotKernel<YU,SL_POT,3,PhiY>{

public:
  typedef PotKernelTraits<PhiY> Trait;

private:
  const typename Trait::MeshY&   meshy;
        typename Trait::MatType  mat;
        typename Trait::JacY     dy;
                        PhiY     phiy;
  const                 Real     kappa;
                        R3       x0_y0,x_y;
                        Real     h,r;
                        Cplx     ker;

public:

  PotKernel<YU,SL_POT,3,PhiY>(const typename Trait::MeshY& my,
			      const Real& k):
  meshy(my), phiy(my), kappa(k) {};

  void Assign(const R3& x, const int& iy){
    const typename Trait::EltY& ey=meshy[iy];
    x0_y0 = x-ey[0];
    dy    = MatJac(ey);
    h     = DetJac(ey);}

  const typename Trait::MatType&
  operator()(const R3& x,const typename Trait::Rdy& tj){
    x_y  = x0_y0-dy*tj;
    r   = norm2(x_y);
    ker = h*exp(-kappa*r)/(4*pi*r);
    for(int k=0; k<Trait::nb_dof_y; k++){
      mat(0,k) = ker*phiy(k,tj);}
    return mat;}

  const Cplx&
  operator()(const R3& x, const typename Trait::Rdy& tj, const int& ky){
    x_y  = x0_y0-dy*tj;
    r   = norm2(x_y);
    ker = h*exp(-kappa*r)/(4*pi*r);
    return ker *= phiy(ky,tj);
  }

};

typedef PotKernel<YU,SL_POT,3,P0_2D> YU_SL_3D_P0;
typedef PotKernel<YU,SL_POT,3,P1_2D> YU_SL_3D_P1;



/*==========================
  DOUBLE COUCHE YUKAWA EN 2D
  ==========================*/

template <typename PhiY>
class PotKernel<YU,DL_POT,2,PhiY>{

public:
  typedef PotKernelTraits<PhiY> Trait;

private:
  const typename Trait::MeshY&        meshy;
        typename Trait::MatType       mat;
        typename Trait::JacY          dy;
                        PhiY          phiy;
  const                 std::vector<R3>&   normaly;
  const                 Real          kappa;
                        R3            x0_y0,x_y,ny;
                        Real          h,r;
                        Cplx          ker;

public:

  PotKernel<YU,DL_POT,2,PhiY>(const typename Trait::MeshY& my,
			      const Real& k):
  meshy(my), phiy(my), kappa(k), normaly(NormalTo(my)) {};

  void Assign(const R3& x, const int& iy){
    const typename Trait::EltY& ey=meshy[iy];
    x0_y0 = x-ey[0];
    dy    = MatJac(ey);
    h     = DetJac(ey);
    ny    = normaly[iy];
  }

  const typename Trait::MatType&
  operator()(const R3& x,const typename Trait::Rdy& tj){
    x_y  = x0_y0-dy*tj;
    r   = norm2(x_y);
    ker = -h*(ny,x_y)*(1./r)*(0.5/pi)*kappa*Modified_BesselK1(kappa*r);
    for(int k=0; k<Trait::nb_dof_y; k++){
      mat(0,k) = ker*phiy(k,tj);}
    return mat;
  }

  const Cplx&
  operator()(const R3& x, const typename Trait::Rdy& tj, const int& ky){
    x_y  = x0_y0-dy*tj;
    r   = norm2(x_y);
    ker = -h*(ny,x_y)*(1./r)*(0.5/pi)*kappa*Modified_BesselK1(kappa*r);
    return ker *= phiy(ky,tj);
  }

};


typedef PotKernel<YU,DL_POT,2,P0_1D> YU_DL_2D_P0;
typedef PotKernel<YU,DL_POT,2,P1_1D> YU_DL_2D_P1;




/*==========================
  DOUBLE COUCHE YUKAWA EN 3D
  ==========================*/

template <typename PhiY>
class PotKernel<YU,DL_POT,3,PhiY>{

public:
  typedef PotKernelTraits<PhiY> Trait;

private:
  const typename Trait::MeshY&        meshy;
        typename Trait::MatType       mat;
        typename Trait::JacY          dy;
                        PhiY          phiy;
  const                 std::vector<R3>&   normaly;
  const                 Real          kappa;
                        R3            x0_y0,x_y,ny;
                        Real          h,r,r3;
                        Cplx          ker;

public:

  PotKernel<YU,DL_POT,3,PhiY>(const typename Trait::MeshY& my,
			      const Real& k):
  meshy(my), phiy(my), kappa(k), normaly(NormalTo(my)) {};

  void Assign(const R3& x, const int& iy){
    const typename Trait::EltY& ey=meshy[iy];
    x0_y0 = x-ey[0];
    dy    = MatJac(ey);
    h     = DetJac(ey);
    ny    = normaly[iy];
  }

  const typename Trait::MatType&
  operator()(const R3& x,const typename Trait::Rdy& tj){
    x_y  = x0_y0-dy*tj;
    r   = norm2(x_y);
    r3  = r*r*r;
    ker = -h*(ny,x_y)*(kappa*r+1.)*exp(-kappa*r)/(4*pi*r3);
    for(int k=0; k<Trait::nb_dof_y; k++){
      mat(0,k) = ker*phiy(k,tj);}
    return mat;}

  const Cplx&
  operator()(const R3& x, const typename Trait::Rdy& tj, const int& ky){
    x_y  = x0_y0-dy*tj;
    r   = norm2(x_y);
    r3  = r*r*r;
    ker = -h*(ny,x_y)*(kappa*r+1.)*exp(-kappa*r)/(4*pi*r3);
    return ker *= phiy(ky,tj);
  }

};

typedef PotKernel<YU,DL_POT,3,P0_2D> YU_DL_3D_P0;
typedef PotKernel<YU,DL_POT,3,P1_2D> YU_DL_3D_P1;
}

#endif
