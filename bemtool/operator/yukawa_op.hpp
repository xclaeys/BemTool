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
#ifndef BEMTOOL_OPERATOR_YUKAWA_HPP
#define BEMTOOL_OPERATOR_YUKAWA_HPP

#include "operator.hpp"

namespace bemtool {
/*==========================
  SIMPLE COUCHE YUKAWA EN 2D
  ==========================*/

template <typename PhiX, typename PhiY>
class BIOpKernel<YU,SL_OP,2,PhiX,PhiY>{

public:
  typedef BIOpKernelTraits<YU,SL_OP,2,PhiX,PhiY> Trait;
  static const int EqType = YU;
  static const int OpType = SL_OP;
  static const int Dim    = 2;

private:
  const typename Trait::MeshX&   meshx;
  const typename Trait::MeshY&   meshy;
        typename Trait::MatType  inter;
        typename Trait::JacX     dx;
        typename Trait::JacY     dy;
                        PhiX     phix;
                        PhiY     phiy;
  const                 Real     kappa;
                        R3       x0_y0,x_y;
                        Real     h,r;
                        Cplx     ker;

public:
  BIOpKernel<YU,SL_OP,2,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k):
  meshx(mx), phix(mx), meshy(my), phiy(my), kappa(k) {};

  void Assign(const int& ix, const int& iy){
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];
    x0_y0 = ex[0]-ey[0];
    dx    = MatJac(ex);
    dy    = MatJac(ey);
    h     = DetJac(ex)*DetJac(ey);
  }

  const typename Trait::MatType& operator()(const typename Trait::Rdx& tx,
					    const typename Trait::Rdy& ty){
    x_y  = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    ker = h*(0.5/pi)*Modified_BesselK0(kappa*r);
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
	inter(j,k) = ker*phix(j,tx)*phiy(k,ty);
      }
    }
    return inter;
  }

  const Cplx& operator()(const typename Trait::Rdx& tx,
			 const typename Trait::Rdy& ty,
			 const int& kx, const int& ky){
    x_y  = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    ker = h*(0.5/pi)*Modified_BesselK0(kappa*r);
    return ker *= phix(kx,tx)*phiy(ky,ty);
  }

};

typedef BIOpKernel<YU,SL_OP,2,P0_1D,P0_1D> YU_SL_2D_P0xP0;
typedef BIOpKernel<YU,SL_OP,2,P1_1D,P1_1D> YU_SL_2D_P1xP1;
typedef BIOpKernel<YU,SL_OP,2,P2_1D,P2_1D> YU_SL_2D_P2xP2;
typedef BIOpKernel<YU,SL_OP,2,P0_1D,P1_1D> YU_SL_2D_P0xP1;
typedef BIOpKernel<YU,SL_OP,2,P1_1D,P0_1D> YU_SL_2D_P1xP0;



/*==========================
  SIMPLE COUCHE YUKAWA EN 3D
  ==========================*/

template <typename PhiX, typename PhiY>
class BIOpKernel<YU,SL_OP,3,PhiX,PhiY>{

public:
  typedef BIOpKernelTraits<YU,SL_OP,3,PhiX,PhiY> Trait;
  static const int EqType = YU;
  static const int OpType = SL_OP;
  static const int Dim    = 3;

private:
  const typename Trait::MeshX&   meshx;
  const typename Trait::MeshY&   meshy;
        typename Trait::MatType  inter;
        typename Trait::JacX     dx;
        typename Trait::JacY     dy;
                        PhiX     phix;
                        PhiY     phiy;
  const                 Real     kappa;
                        R3       x0_y0,x_y;
                        Real     h,r;
                        Cplx     ker;

public:
  BIOpKernel<YU,SL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k):
  meshx(mx), phix(mx), meshy(my), phiy(my), kappa(k) {};

  void Assign(const int& ix, const int& iy){
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];
    x0_y0 = ex[0]-ey[0];
    dx    = MatJac(ex);
    dy    = MatJac(ey);
    h     = DetJac(ex)*DetJac(ey);
  }

  const typename Trait::MatType& operator()(const typename Trait::Rdx& tx,
					    const typename Trait::Rdy& ty){
    x_y = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    ker = h*exp(-kappa*r)/(4*pi*r);
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
	inter(j,k) = ker*phix(j,tx)*phiy(k,ty);
      }
    }
    return inter;
  }

  const Cplx& operator()(const typename Trait::Rdx& tx,
			 const typename Trait::Rdy& ty,
			 const int& kx, const int& ky){
    x_y = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    ker = h*exp(-kappa*r)/(4*pi*r);
    return ker *= phix(kx,tx)*phiy(ky,ty);
  }


};

typedef BIOpKernel<YU,SL_OP,3,P0_2D,P0_2D> YU_SL_3D_P0xP0;
typedef BIOpKernel<YU,SL_OP,3,P1_2D,P1_2D> YU_SL_3D_P1xP1;
typedef BIOpKernel<YU,SL_OP,3,P2_2D,P2_2D> YU_SL_3D_P2xP2;
typedef BIOpKernel<YU,SL_OP,3,P0_2D,P1_2D> YU_SL_3D_P0xP1;
typedef BIOpKernel<YU,SL_OP,3,P1_2D,P0_2D> YU_SL_3D_P1xP0;


/*==========================
  DOUBLE COUCHE YUKAWA EN 2D
  ==========================*/

template <typename PhiX, typename PhiY>
class BIOpKernel<YU,DL_OP,2,PhiX,PhiY>{

public:
  typedef BIOpKernelTraits<YU,DL_OP,2,PhiX,PhiY> Trait;
  static const int EqType = YU;
  static const int OpType = DL_OP;
  static const int Dim    = 2;

private:
  const typename Trait::MeshX&   meshx;
  const typename Trait::MeshY&   meshy;
  const            std::vector<R3>&   normaly;
        typename Trait::MatType  inter;
        typename Trait::JacX     dx;
        typename Trait::JacY     dy;
                        PhiX     phix;
                        PhiY     phiy;
  const                 Real     kappa;
                        R3       x0_y0,x_y,ny;
                        Real     h,r;
                        Cplx     ker;

public:
  BIOpKernel<YU,DL_OP,2,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k):
  meshx(mx), phix(mx), meshy(my), phiy(my), kappa(k),
    normaly(NormalTo(my)) {};

  void Assign(const int& ix, const int& iy){
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];
    x0_y0 = ex[0]-ey[0];
    dx    = MatJac(ex);
    dy    = MatJac(ey);
    h     = DetJac(ex)*DetJac(ey);
    ny    = normaly[iy];
  }

  const typename Trait::MatType& operator()(const typename Trait::Rdx& tx,
					    const typename Trait::Rdy& ty){
    x_y = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    ker = -h*(ny,x_y)*(1./r)*(0.5/pi)*kappa*Modified_BesselK1(kappa*r);
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
	inter(j,k) = ker*phix(j,tx)*phiy(k,ty);
      }
    }
    return inter;
  }

  const Cplx& operator()(const typename Trait::Rdx& tx,
			 const typename Trait::Rdy& ty,
			 const int& kx, const int& ky){
    x_y = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    ker = -h*(ny,x_y)*(1./r)*(0.5/pi)*kappa*Modified_BesselK1(kappa*r);
    return ker *= phix(kx,tx)*phiy(ky,ty);
  }


};


typedef BIOpKernel<YU,DL_OP,2,P0_1D,P0_1D> YU_DL_2D_P0xP0;
typedef BIOpKernel<YU,DL_OP,2,P1_1D,P1_1D> YU_DL_2D_P1xP1;
typedef BIOpKernel<YU,DL_OP,2,P2_1D,P2_1D> YU_DL_2D_P2xP2;
typedef BIOpKernel<YU,DL_OP,2,P0_1D,P1_1D> YU_DL_2D_P0xP1;
typedef BIOpKernel<YU,DL_OP,2,P1_1D,P0_1D> YU_DL_2D_P1xP0;



/*==========================
  DOUBLE COUCHE YUKAWA EN 3D
  ==========================*/

template <typename PhiX, typename PhiY>
class BIOpKernel<YU,DL_OP,3,PhiX,PhiY>{

public:
  typedef BIOpKernelTraits<YU,DL_OP,3,PhiX,PhiY> Trait;
  static const int EqType = YU;
  static const int OpType = DL_OP;
  static const int Dim    = 3;

private:
  const typename Trait::MeshX&        meshx;
  const typename Trait::MeshY&        meshy;
        typename Trait::MatType       inter;
        typename Trait::JacX          dx;
        typename Trait::JacY          dy;
                        PhiX          phix;
                        PhiY          phiy;
  const                 std::vector<R3>&   normaly;
  const                 Real          kappa;
                        R3            x0_y0,x_y,ny;
                        Real          h,r,r3;
                        Cplx          ker;

public:
  BIOpKernel<YU,DL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k):
  meshx(mx), phix(mx), meshy(my), phiy(my), kappa(k),
    normaly(NormalTo(my)) {};

  void Assign(const int& ix, const int& iy){
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];
    x0_y0 = ex[0]-ey[0];
    dx    = MatJac(ex);
    dy    = MatJac(ey);
    h     = DetJac(ex)*DetJac(ey);
    ny    = normaly[iy];
  }

  const typename Trait::MatType& operator()(const typename Trait::Rdx& tx,
					    const typename Trait::Rdy& ty){
    x_y = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    r3  = r*r*r;
    ker = -h*(ny,x_y)*(kappa*r+1.)*exp(-kappa*r)/(4*pi*r3);
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
	inter(j,k) = ker*phix(j,tx)*phiy(k,ty);
      }
    }
    return inter;
  }

  const Cplx& operator()(const typename Trait::Rdx& tx,
			 const typename Trait::Rdy& ty,
			 const int& kx, const int& ky){
    x_y = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    r3  = r*r*r;
    ker = -h*(ny,x_y)*(kappa*r+1.)*exp(-kappa*r)/(4*pi*r3);
    return ker *= phix(kx,tx)*phiy(ky,ty);
  }

};

typedef BIOpKernel<YU,DL_OP,3,P0_2D,P0_2D> YU_DL_3D_P0xP0;
typedef BIOpKernel<YU,DL_OP,3,P1_2D,P1_2D> YU_DL_3D_P1xP1;
typedef BIOpKernel<YU,DL_OP,3,P2_2D,P2_2D> YU_DL_3D_P2xP2;
typedef BIOpKernel<YU,DL_OP,3,P0_2D,P1_2D> YU_DL_3D_P0xP1;
typedef BIOpKernel<YU,DL_OP,3,P1_2D,P0_2D> YU_DL_3D_P1xP0;



/*==================================
  DOUBLE COUCHE ADJOINT YUKAWA EN 2D
  ==================================*/

template <typename PhiX, typename PhiY>
class BIOpKernel<YU,TDL_OP,2,PhiX,PhiY>{

public:
  typedef BIOpKernelTraits<YU,TDL_OP,2,PhiX,PhiY> Trait;
  static const int EqType = YU;
  static const int OpType = TDL_OP;
  static const int Dim    = 2;

private:
  const typename Trait::MeshX&        meshx;
  const typename Trait::MeshY&        meshy;
        typename Trait::MatType       inter;
        typename Trait::JacX          dx;
        typename Trait::JacY          dy;
                        PhiX          phix;
                        PhiY          phiy;
  const                 std::vector<R3>&   normalx;
  const                 Real          kappa;
                        R3            x0_y0,x_y,nx;
                        Real          h,r;
                        Cplx          ker;

public:
  BIOpKernel<YU,TDL_OP,2,PhiX,PhiY>(const typename Trait::MeshX& mx,
				    const typename Trait::MeshY& my,
				    const Real& k):
  meshx(mx), phix(mx), meshy(my), phiy(my), kappa(k),
    normalx(NormalTo(mx)) {};

  void Assign(const int& ix, const int& iy){
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];
    x0_y0 = ex[0]-ey[0];
    dx    = MatJac(ex);
    dy    = MatJac(ey);
    h     = DetJac(ex)*DetJac(ey);
    nx    = normalx[ix];
  }

  const typename Trait::MatType& operator()(const typename Trait::Rdx& tx,
					    const typename Trait::Rdy& ty){
    x_y = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    ker = -h*(nx,x_y)*(1./r)*(0.5/pi)*kappa*Modified_BesselK1(kappa*r);
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
	inter(j,k) = ker*phix(j,tx)*phiy(k,ty);
      }
    }
    return inter;
  }

  const Cplx& operator()(const typename Trait::Rdx& tx,
			 const typename Trait::Rdy& ty,
			 const int& kx, const int& ky){
    x_y = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    ker = -h*(nx,x_y)*(1./r)*(0.5/pi)*kappa*Modified_BesselK1(kappa*r);
    return ker *= phix(kx,tx)*phiy(ky,ty);
  }


};

typedef BIOpKernel<YU,TDL_OP,2,P0_1D,P0_1D> YU_TDL_2D_P0xP0;
typedef BIOpKernel<YU,TDL_OP,2,P1_1D,P1_1D> YU_TDL_2D_P1xP1;
typedef BIOpKernel<YU,TDL_OP,2,P2_1D,P2_1D> YU_TDL_2D_P2xP2;
typedef BIOpKernel<YU,TDL_OP,2,P1_1D,P0_1D> YU_TDL_2D_P1xP0;
typedef BIOpKernel<YU,TDL_OP,2,P0_1D,P1_1D> YU_TDL_2D_P0xP1;



/*==================================
  DOUBLE COUCHE ADJOINT YUKAWA EN 3D
  ==================================*/

template <typename PhiX, typename PhiY>
class BIOpKernel<YU,TDL_OP,3,PhiX,PhiY>{

public:
  typedef BIOpKernelTraits<YU,TDL_OP,3,PhiX,PhiY> Trait;
  static const int EqType = YU;
  static const int OpType = TDL_OP;
  static const int Dim    = 3;

private:
  const typename Trait::MeshX&        meshx;
  const typename Trait::MeshY&        meshy;
        typename Trait::MatType       inter;
        typename Trait::JacX          dx;
        typename Trait::JacY          dy;
                        PhiX          phix;
                        PhiY          phiy;
  const                 std::vector<R3>&   normalx;
  const                 Real          kappa;
                        R3            x0_y0,x_y,nx;
                        Real          h,r,r3;
                        Cplx          ker;

public:
  BIOpKernel<YU,TDL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				    const typename Trait::MeshY& my,
				    const Real& k):
  meshx(mx), phix(mx), meshy(my), phiy(my), kappa(k),
    normalx(NormalTo(mx)) {};

  void Assign(const int& ix, const int& iy){
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];
    x0_y0 = ex[0]-ey[0];
    dx    = MatJac(ex);
    dy    = MatJac(ey);
    h     = DetJac(ex)*DetJac(ey);
    nx    = normalx[ix];
  }

  const typename Trait::MatType& operator()(const typename Trait::Rdx& tx,
					    const typename Trait::Rdy& ty){
    x_y = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    r3  = r*r*r;
    ker = -h*(nx,x_y)*(kappa*r+1.)*exp(-kappa*r)/(4*pi*r3);
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
	inter(j,k) = ker*phix(j,tx)*phiy(k,ty);
      }
    }
    return inter;
  }


  const Cplx& operator()(const typename Trait::Rdx& tx,
			 const typename Trait::Rdy& ty,
			 const int& kx, const int& ky){
    x_y = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    r3  = r*r*r;
    ker = -h*(nx,x_y)*(kappa*r+1.)*exp(-kappa*r)/(4*pi*r3);
    return ker *= phix(kx,tx)*phiy(ky,ty);
  }


};

typedef BIOpKernel<YU,TDL_OP,3,P0_2D,P0_2D> YU_TDL_3D_P0xP0;
typedef BIOpKernel<YU,TDL_OP,3,P1_2D,P1_2D> YU_TDL_3D_P1xP1;
typedef BIOpKernel<YU,TDL_OP,3,P2_2D,P2_2D> YU_TDL_3D_P2xP2;
typedef BIOpKernel<YU,TDL_OP,3,P0_2D,P1_2D> YU_TDL_3D_P0xP1;
typedef BIOpKernel<YU,TDL_OP,3,P1_2D,P0_2D> YU_TDL_3D_P1xP0;






/*===========================
  HYPERSINGULIER YUKAWA EN 2D
  ===========================*/

template <typename PhiX, typename PhiY>
class BIOpKernel<YU,HS_OP,2,PhiX,PhiY>{

public:
  typedef BIOpKernelTraits<YU,HS_OP,2,PhiX,PhiY> Trait;
  static const int EqType = YU;
  static const int OpType = HS_OP;
  static const int Dim    = 2;

private:
  const typename Trait::MeshX&       meshx;
  const typename Trait::MeshY&       meshy;
        typename Trait::MatType      inter;
        typename Trait::JacX         dx;
        typename Trait::JacY         dy;
        typename Trait::GradPhiX     grad_phix;
        typename Trait::GradPhiY     grad_phiy;
                        PhiX         phix;
                        PhiY         phiy;
  const                 std::vector<R3>&  normalx;
  const                 std::vector<R3>&  normaly;
  const                 Real         kappa,kappa2;
                        R3           x0_y0,x_y,nx,ny;
                        Real         h,r;
                        Cplx         ker,val,val2;

public:
  BIOpKernel<YU,HS_OP,2,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k):
  meshx(mx), phix(mx), grad_phix(mx), normalx(NormalTo(mx)),
    meshy(my), phiy(my), grad_phiy(my), normaly(NormalTo(my)),
    kappa(k), kappa2(kappa*kappa) {};

  void Assign(const int& ix, const int& iy){
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];
    grad_phix.Assign(ix);
    grad_phiy.Assign(iy);
    h     = DetJac(ex)*DetJac(ey);
    x0_y0 = ex[0]-ey[0];
    dx    = MatJac(ex);
    dy    = MatJac(ey);
    nx    = normalx[ix];
    ny    = normaly[iy];
  }

  const typename Trait::MatType& operator()(const typename Trait::Rdx& tx,
					    const typename Trait::Rdy& ty){
    x_y = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    ker = h*(0.5/pi)*Modified_BesselK0(kappa*r);
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
	val = ( vprod(nx,grad_phix(j,tx)),vprod(ny,grad_phiy(k,ty)) )*ker;
	inter(j,k) = val + kappa2*(nx,ny)*phix(j,tx)*phiy(k,ty)*ker;
      }
    }
    return inter;
  }

  const Cplx& operator()(const typename Trait::Rdx& tx,
			 const typename Trait::Rdy& ty,
			 const int& kx, const int& ky){
    x_y = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    ker = h*(0.5/pi)*Modified_BesselK0(kappa*r);
    val = ( vprod(nx,grad_phix(kx,tx)),vprod(ny,grad_phiy(ky,ty)) )*ker;
    return val2 = val+kappa2*(nx,ny)*phix(kx,tx)*phiy(ky,ty)*ker;
  }


};


typedef BIOpKernel<YU,HS_OP,2,P1_1D,P1_1D> YU_HS_2D_P1xP1;
typedef BIOpKernel<YU,HS_OP,2,P2_1D,P2_1D> YU_HS_2D_P2xP2;


/*===========================
  HYPERSINGULIER YUKAWA EN 3D
  ===========================*/

template <typename PhiX, typename PhiY>
class BIOpKernel<YU,HS_OP,3,PhiX,PhiY>{

public:
  typedef BIOpKernelTraits<YU,HS_OP,3,PhiX,PhiY> Trait;

private:
  const typename Trait::MeshX&       meshx;
  const typename Trait::MeshY&       meshy;
        typename Trait::MatType      inter;
        typename Trait::JacX         dx;
        typename Trait::JacY         dy;
        typename Trait::GradPhiX     grad_phix;
        typename Trait::GradPhiY     grad_phiy;
                        PhiX         phix;
                        PhiY         phiy;
  const                 std::vector<R3>&  normalx;
  const                 std::vector<R3>&  normaly;
  const                 Real         kappa,kappa2;
                        R3           x0_y0,x_y,nx,ny;
                        Real         h,r;
                        Cplx         ker,val,val2;

public:
  BIOpKernel<YU,HS_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k):
  meshx(mx), phix(mx), grad_phix(mx), normalx(NormalTo(mx)),
    meshy(my), phiy(my), grad_phiy(my), normaly(NormalTo(my)),
    kappa(k), kappa2(kappa*kappa) {};

  void Assign(const int& ix, const int& iy){
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];
    grad_phix.Assign(ix);
    grad_phiy.Assign(iy);
    h     = DetJac(ex)*DetJac(ey);
    x0_y0 = ex[0]-ey[0];
    dx    = MatJac(ex);
    dy    = MatJac(ey);
    nx    = normalx[ix];
    ny    = normaly[iy];
  }

  const typename Trait::MatType& operator()(const typename Trait::Rdx& tx,
					    const typename Trait::Rdy& ty){
    x_y = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    ker = h*exp(-kappa*r)/(4*pi*r);
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
	val = ( vprod(nx,grad_phix(j,tx)),vprod(ny,grad_phiy(k,ty)) )*ker;
	inter(j,k) = val + kappa2*(nx,ny)*phix(j,tx)*phiy(k,ty)*ker;
      }
    }
    return inter;
  }

  const Cplx& operator()(const typename Trait::Rdx& tx,
			 const typename Trait::Rdy& ty,
			 const int& kx, const int& ky){
    x_y = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    ker = h*exp(-kappa*r)/(4*pi*r);
    val = ( vprod(nx,grad_phix(kx,tx)),vprod(ny,grad_phiy(ky,ty)) )*ker;
    return val2 = val + kappa2*(nx,ny)*phix(kx,tx)*phiy(ky,ty)*ker;
  }


};

typedef BIOpKernel<YU,HS_OP,3,P1_2D,P1_2D> YU_HS_3D_P1xP1;
typedef BIOpKernel<YU,HS_OP,3,P2_2D,P2_2D> YU_HS_3D_P2xP2;

}



#endif
