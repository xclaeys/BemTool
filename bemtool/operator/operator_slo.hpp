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
#ifndef BEMTOOL_OPERATOR_OPERATOR_SLO_HPP
#define BEMTOOL_OPERATOR_OPERATOR_SLO_HPP

#include <fstream>
#include "../quadrature/quad_bem.hpp"
#include "../fem/shapefct.hpp"
#include "../equations.hpp"
#include "../miscellaneous/specialfct.hpp"
#include "block_op.hpp"

namespace bemtool {


template <int D, typename PhiX, typename PhiY> class BIOp_SLO{

public:

  typedef Mesh<D>                   MeshX;
  typedef Mesh<D>                   MeshY;
  typedef typename PhiX::Rd            RdX;
  typedef typename PhiY::Rd            RdY;
  typedef Elt<D>                    EltX;
  typedef Elt<D>                    EltY;
  typedef BlockMat  MatType;

  static  const int dimx = D;
  static  const int dimy = D;
  static const int nb_dof_x = PhiX::nb_dof_loc;
  static const int nb_dof_y = PhiY::nb_dof_loc;
  typedef mat<dimx,dimx,Real>                JacX_ref;
  typedef mat<dimy,dimy,Real>                JacY_ref;
  typedef typename Jac<dimx>::Type     JacX;
  typedef typename Jac<dimy>::Type     JacY;
  typedef QuadBEM<dimx,dimy>                 QuadType;

private:

  const MeshX&    meshx;
  const MeshY&    meshy;
  const Geometry& nodex;
  const Geometry& nodey;
  PhiX          phix;
  PhiY          phiy;

  QuadType     qr;
  MatType      inter;
  int          rule;
  JacX_ref         dx_ref;
  JacY_ref         dy_ref;
  JacX         dx;
  JacY         dy;
  RdX          x0_ref,x_ref,ax[dimx+1];
  RdY          y0_ref,y_ref,ay[dimy+1];
  R3           x0_y0,x_y;
  Real         h,r,r_D;
  std::vector<std::vector<int>> slo_num_ref;
  std::vector<int> slo_num;

public:
  BIOp_SLO(const MeshX& mx, const MeshY& my):
  meshx(mx), nodex(GeometryOf(mx)),
    meshy(my), nodey(GeometryOf(my)),phix(mx), phiy(my),qr(6) {
    for(int j=0; j<dimx; j++){ax[j+1][j]=1;}
    for(int j=0; j<dimy; j++){ay[j+1][j]=1;}

    for (int j=0;j<D+2;j++){
      // std::cout << nb_dof_x+nb_dof_y<<std::endl;
      std::vector<int> temp(nb_dof_x+nb_dof_y);
      for (int k=0;k<nb_dof_x;k++){
        temp[k]=k;
      }
      if (j==0){ // Disjoint
        for (int k=nb_dof_x;k<nb_dof_x+nb_dof_y;k++){
          temp[k]=k;
        }
      }
      else if (j==1){ // 1 node
        temp[nb_dof_x]=0;
        for (int k=nb_dof_x+1;k<nb_dof_x+nb_dof_y;k++){
          temp[k]=k-1;
        }
      }
      else if (j==2){ // 2 nodes
        temp[nb_dof_x]=0;temp[nb_dof_x+1]=1;
        for (int k=nb_dof_x+2;k<nb_dof_x+nb_dof_y;k++){
          temp[k]=k-2;
        }
      }
      else if (j==3){ // identical
        for (int k=0;k<nb_dof_y;k++){
          temp[k+nb_dof_x]=k;
        }
      }
      slo_num_ref.push_back(temp);
    }

  };

  void ChooseQuad(const int& jx, const int& jy){
    rule = 0; slo_num.clear();
    EltX ex = meshx[jx];   EltY ey = meshy[jy];
    Elt<dimx,RdX> bx(ax);  Elt<dimy,RdY> by(ay);
    std::vector<std::pair<int,int>> permutation_x;
    std::vector<std::pair<int,int>> permutation_y;

    for(int p=0; p<dimx+1; p++){
      for(int q=rule; q<dimy+1; q++){
        if( &ex[p]==&ey[q] ){
          Swap(ex,rule,p); Swap(bx,rule,p);
          Swap(ey,rule,q); Swap(by,rule,q);
          permutation_x.emplace_back(rule,p);
          permutation_y.emplace_back(rule,q);
          rule++; break;
        }
      }
    }
    dx_ref = MatJac(bx); x0_ref=bx[0];
    dy_ref = MatJac(by); y0_ref=by[0];

    slo_num=slo_num_ref[rule];
    for (int j=permutation_x.size()-1;j>=0;j--){
      int temp;
      temp=slo_num[permutation_x[j].first];
      slo_num[permutation_x[j].first]=slo_num[permutation_x[j].second];
      slo_num[permutation_x[j].second]=temp;
    }

    for (int j=permutation_y.size()-1;j>=0;j--){
      int temp;
      temp=slo_num[permutation_y[j].first+nb_dof_x];
      slo_num[permutation_y[j].first+nb_dof_x]=slo_num[permutation_y[j].second+nb_dof_x];
      slo_num[permutation_y[j].second+nb_dof_x]=temp;
    }
  }

  const std::vector<int>& get_slo_num() const {return slo_num;}


  const MatType& operator()(const EltX& ex, const EltY& ey){
    return (*this)(nodex[ex][meshx],nodey[ey][meshy]);}


  const MatType& operator()(const int& jx, const int& jy){
    ChooseQuad(jx,jy);
    const std::vector<RdX>&  s = qr.x(rule);
    const std::vector<RdY>&  t = qr.y(rule);
    const std::vector<Real>& w = qr.w(rule);
    inter.Clear();
    inter.Resize(nb_dof_x+nb_dof_y-rule,nb_dof_x+nb_dof_y-rule);

    // Geometry
    const EltX& ex=meshx[jx];
    const EltY& ey=meshy[jy];
    h     = DetJac(ex)*DetJac(ey);
    x0_y0 = ex[0]-ey[0];
    dx    = MatJac(ex);
    dy    = MatJac(ey);
    // std::cout << "WARNING"<<std::endl;
    // Quadrature
    for(int j=0; j<w.size(); j++){
      // std::cout << "quad "<<j<<std::endl;
      x_ref = x0_ref + dx_ref*s[j];
      y_ref = y0_ref + dy_ref*t[j];
      BlockMat test = w[j]*ker(x_ref,y_ref);
      inter += w[j]*ker(x_ref,y_ref);
    }
    return inter;
  }

  MatType ker(const RdX& tx, const RdY& ty){
    BlockMat inter_loc(NbRow(inter),NbCol(inter));
    // std::cout << size(x_y)<<" "<<size(x0_y0)<<" "<<NbRow(dx)<<
    x_y = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    r_D = std::pow(r,D+1);


    for(int j=0; j<nb_dof_x; j++){
      for(int k=0; k<nb_dof_x; k++){
	      inter_loc(slo_num[j],slo_num[k])+=phix(j,tx)*phix(k,tx);
      }
    }

    for(int j=0; j<nb_dof_y; j++){
      for(int k=0; k<nb_dof_y; k++){
	      inter_loc(slo_num[j+nb_dof_x],slo_num[k+nb_dof_x])+=phiy(j,ty)*phiy(k,ty);
      }
    }

    for(int j=0; j<nb_dof_x; j++){
      for(int k=0; k<nb_dof_y; k++){
	      inter_loc(slo_num[j],slo_num[k+nb_dof_x])+= - phix(j,tx)*phiy(k,ty);
      }
    }

    for(int j=0; j<nb_dof_x; j++){
      for(int k=0; k<nb_dof_y; k++){
	      inter_loc(slo_num[j+nb_dof_x],slo_num[k])+= - phiy(j,ty)*phix(k,tx);
      }
    }
    return (h/r_D)*inter_loc;
  }
  

  // const Cplx& operator()(const N2& Ix, const N2& Iy){
  //   ChooseQuad(Ix[0],Iy[0]);
  //   const std::vector<RdX>&  s = qr.x(rule);
  //   const std::vector<RdY>&  t = qr.y(rule);
  //   const std::vector<Real>& w = qr.w(rule);
  //   ker.Assign(Ix[0],Iy[0]);
  //   for(int j=0; j<w.size(); j++){
  //     x = x0 + dx_ref*s[j]; y = y0 + dy_ref*t[j];
  //     val += w[j]*ker(x,y,Ix[1],Iy[1]);}
  //   return val;
  // }


  // const Cplx& operator()(const std::vector<N2>& vx,
	// 		 const std::vector<N2>& vy){
  //   val2=0.;
  //   for(int ix=0; ix<vx.size(); ix++){
  //     for(int iy=0; iy<vy.size(); iy++){
	// val2 += (*this)(vx[ix],vy[iy]);} }
  //   return val2;
  // }


};

}

#endif
