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

#ifndef BEMTOOL_POTENTIAL_MAXWELL_HPP
#define BEMTOOL_POTENTIAL_MAXWELL_HPP

#include "potential.hpp"
namespace bemtool {
  
  /*===========================
    SIMPLE COUCHE MAXWELL EN 3D
    ===========================*/
  
  template <typename PhiY>
  class PotKernel<MA,SL_POT,3,PhiY>{
    
  public:
    typedef PotKernelTraits<PhiY,3>      Trait;
    
  private:
    const typename Trait::MeshY&   meshy;
          typename Trait::MatType  inter;
          typename Trait::JacY     dy;
          typename Trait::DivPhiY  div_phiy;
                          PhiY     phiy;
    const                 Real     kappa, kappa2;
                          R3       x0_y0,x_y,nx,ny;
                          Real     h,r;
                          Cplx     ker,val,val2;
                          C3       pot;
    
  public:
    PotKernel<MA,SL_POT,3,PhiY>(const typename Trait::MeshY& my,
				const Real& k):
    meshy(my), phiy(my), div_phiy(my), kappa(k), kappa2(kappa*kappa) {};

    void Assign(const R3& x, const int& iy){
      const typename Trait::EltY& ey=meshy[iy];
      phiy.Assign(iy);
      div_phiy.Assign(iy);
      x0_y0 = x-ey[0];
      dy    = MatJac(ey);
      h     = DetJac(ey);      
    }
    
    const typename Trait::MatType& operator()(const R3& x,const typename Trait::Rdy& ty){
      x_y  = x0_y0-dy*ty;
      r   = norm2(x_y);
      ker = h*exp(iu*kappa*r)/(4*pi*r);
      val = (iu*kappa*r-1)/(kappa2*r*r);
      for(int k=0; k<Trait::nb_dof_y; k++){
	pot  = ker*(phiy(k,ty)+div_phiy(k,ty)*val*x_y);
	inter(0,k) = pot[0];
	inter(1,k) = pot[1];
	inter(2,k) = pot[2];
      }            
      return inter;
    }
    
    const C3&
    operator()(const R3& x, const typename Trait::Rdy& ty, const int& ky){
      x_y  = x0_y0-dy*ty;
      r   = norm2(x_y);
      ker = h*exp(iu*kappa*r)/(4*pi*r);
      val = (iu*kappa*r-1)/(kappa2*r*r);
      pot  = ker*( phiy(ky,ty) + div_phiy(ky,ty)*val*x_y );
      return pot;
    }
    
  };

  typedef PotKernel<MA,SL_POT,3,RT0_2D>  EFIE_RT0;

  /*===========================
    DOUBLE COUCHE MAXWELL EN 3D
    ===========================*/
  
  template <typename PhiY>
  class PotKernel<MA,DL_POT,3,PhiY>{
    
  public:
    typedef PotKernelTraits<PhiY,3>      Trait;

    
  private:
    const typename Trait::MeshY&   meshy;
          typename Trait::MatType  inter;
          typename Trait::JacY     dy;
          typename Trait::DivPhiY  div_phiy;
                          PhiY     phiy;
    const                 Real     kappa, kappa2;
                          R3       x0_y0,x_y,nx,ny;
                          Real     h,r;
                          Cplx     ker,val,val2;
                          C3       pot;
    
  public:
    PotKernel<MA,DL_POT,3,PhiY>(const typename Trait::MeshY& my,
				const Real& k):
    meshy(my), phiy(my), div_phiy(my), kappa(k), kappa2(kappa*kappa) {};

    void Assign(const R3& x, const int& iy){
      const typename Trait::EltY& ey=meshy[iy];
      phiy.Assign(iy);
      div_phiy.Assign(iy);
      x0_y0 = x-ey[0];
      dy    = MatJac(ey);
      h     = DetJac(ey);      
    }
    
    const typename Trait::MatType& operator()(const R3& x,const typename Trait::Rdy& ty){
      x_y  = x0_y0-dy*ty;
      r   = norm2(x_y);
      ker = h*exp(iu*kappa*r)/(4*pi*r);
      val = (iu*kappa*r-1)/(r*r);
      for(int k=0; k<Trait::nb_dof_y; k++){
	pot  = ker*val*vprod(x_y,phiy(k,ty));
	inter(0,k) = pot[0];
	inter(1,k) = pot[1];
	inter(2,k) = pot[2];
      }      
      return inter;
    }
    
    const C3&
    operator()(const R3& x, const typename Trait::Rdy& ty, const int& ky){
      x_y  = x0_y0-dy*ty;
      r   = norm2(x_y);
      ker = h*exp(iu*kappa*r)/(4*pi*r);
      val = (iu*kappa*r-1)/(r*r);
      pot = ker*val*vprod(x_y,phiy(ky,ty));
      return pot;
    }
    
  };

  typedef PotKernel<MA,DL_POT,3,RT0_2D>  MFIE_RT0;

  
}
#endif
