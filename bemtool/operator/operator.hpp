#ifndef OPERATOR_HPP
#define OPERATOR_HPP

#include <fstream>
#include "../mesh/mesh.hpp"
#include "../quadrature/quad_bem.hpp"

using namespace std;
using bemtool::array;

template <typename KernelType> class BIOp{
  
public:

  typedef typename KernelType::Trait         KernelTypeTrait;
  typedef typename KernelTypeTrait::MeshX    MeshX;
  typedef typename KernelTypeTrait::MeshY    MeshY;  
  typedef typename KernelTypeTrait::Rdx      RdX;
  typedef typename KernelTypeTrait::Rdy      RdY;    
  typedef typename KernelTypeTrait::EltX     EltX;
  typedef typename KernelTypeTrait::EltY     EltY;
  typedef typename KernelTypeTrait::MatType  MatType; 

  static  const int dimx = KernelTypeTrait::dimx;
  static  const int dimy = KernelTypeTrait::dimy;   
  typedef mat<dimx,dimx,Real>                JacX;
  typedef mat<dimy,dimy,Real>                JacY;
  typedef QuadBEM<dimx,dimy>                 QuadType;
    
private: 
  
  const MeshX&    meshx;
  const MeshY&    meshy;  
  const Geometry& nodex;
  const Geometry& nodey;
  
  KernelType   ker;
  QuadType     qr;
  MatType      inter;
  int          rule;
  JacX         dx;
  JacY         dy; 
  RdX          x0,x,ax[dimx+1];
  RdY          y0,y,ay[dimy+1];  
  Cplx         val,val2;
  
public:  
  BIOp<KernelType>(const MeshX& mx, const MeshY& my, const Real& k):
  meshx(mx), nodex(GeometryOf(mx)), ker(mx,my,k),
    meshy(my), nodey(GeometryOf(my)), qr(6) { 
    for(int j=0; j<dimx; j++){ax[j+1][j]=1;}
    for(int j=0; j<dimy; j++){ay[j+1][j]=1;}
  };
  
  void ChooseQuad(const int& jx, const int& jy){
    rule = 0; inter=0.; val=0.;
    EltX ex = meshx[jx];   EltY ey = meshx[jy];
    Elt<dimx,RdX> bx(ax);  Elt<dimy,RdY> by(ay);
    
    for(int p=0; p<dimx+1; p++){
      for(int q=rule; q<dimy+1; q++){
	if( &ex[p]==&ey[q] ){
	  Swap(ex,rule,p); Swap(bx,rule,p);  
	  Swap(ey,rule,q); Swap(by,rule,q);  
	  rule++; break;
	}
      }
    }
    dx = MatJac(bx); x0=bx[0];
    dy = MatJac(by); y0=by[0];
  }
  
  
  const MatType& operator()(const EltX& ex, const EltY& ey){        
    return (*this)(nodex[ex][meshx],nodey[ey][meshy]);}
  
  
  const MatType& operator()(const int& jx, const int& jy){        
    ChooseQuad(jx,jy);    
    const vector<RdX>&  s = qr.x(rule);
    const vector<RdY>&  t = qr.y(rule);    
    const vector<Real>& w = qr.w(rule);    
    ker.Assign(jx,jy);
    for(int j=0; j<w.size(); j++){
      x = x0 + dx*s[j];
      y = y0 + dy*t[j];
      inter += w[j]*ker(x,y);
    }
    return inter;
  }


  const Cplx& operator()(const N2& Ix, const N2& Iy){        
    ChooseQuad(Ix[0],Iy[0]);    
    const vector<RdX>&  s = qr.x(rule);
    const vector<RdY>&  t = qr.y(rule);    
    const vector<Real>& w = qr.w(rule);    
    ker.Assign(Ix[0],Iy[0]);
    for(int j=0; j<w.size(); j++){
      x = x0 + dx*s[j]; y = y0 + dy*t[j];
      val += w[j]*ker(x,y,Ix[1],Iy[1]);}
    return val;
  }
  
  
  const Cplx& operator()(const vector<N2>& vx,
			 const vector<N2>& vy){        
    val2=0.;
    for(int ix=0; ix<vx.size(); ix++){
      for(int iy=0; iy<vy.size(); iy++){
	val2 += (*this)(vx[ix],vy[iy]);} }
    return val2;
  }
  
  
};


/*==============================
  GENERALITES SUR LES NOYAUX
  =============================*/

template <int,int,int,typename,typename>
class BIOpKernel;

template <int EqT, int OpT, int dim, typename PhiX, typename PhiY>
class BIOpKernelTraits{

public:
  static const int EqType   = EqT;
  static const int OpType   = OpT;
  static const int Dim      = dim;
  static const int dimx     = PhiX::dim;
  static const int dimy     = PhiY::dim;
  static const int nb_dof_x = PhiX::nb_dof_loc;
  static const int nb_dof_y = PhiY::nb_dof_loc;
  
  typedef typename Jac<dimx>::Type     JacX;
  typedef typename Jac<dimy>::Type     JacY;
  typedef typename PhiX::Rd            Rdx;
  typedef typename PhiY::Rd            Rdy;  
  typedef Mesh<dimx>                   MeshX;
  typedef Mesh<dimy>                   MeshY;
  typedef Elt<dimx>                    EltX;
  typedef Elt<dimy>                    EltY;
  typedef GradBasisFct<PhiX>           GradPhiX;
  typedef GradBasisFct<PhiY>           GradPhiY;
  typedef DivBasisFct<PhiX>            DivPhiX;
  typedef DivBasisFct<PhiY>            DivPhiY;  
  typedef mat<nb_dof_x,nb_dof_y,Cplx>  MatType;    
  
};


/*=======================
        NOYAU CONSTANT
  =======================*/

template <int D, typename PhiX, typename PhiY>
class BIOpKernel<CST,CST_OP,D,PhiX,PhiY>{

public: 
  typedef BIOpKernelTraits<CST,CST_OP,D,PhiX,PhiY> Trait;
  
private:
  const typename Trait::MeshX&    meshx;
  const typename Trait::MeshY&    meshy;          
        typename Trait::MatType   inter;
                        PhiX      phix;  
                        PhiY      phiy;
                        Real      ker;
                        Cplx      val;
  
public:
  BIOpKernel<CST,CST_OP,D,PhiX,PhiY>(const typename Trait::MeshX& mx,
				 const typename Trait::MeshY& my,
				 const Real& k=0.):
  meshx(mx), phix(mx), meshy(my), phiy(my) {};
  
  
  void Assign(const int& ix, const int& iy){
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];   
    ker = DetJac(ex)*DetJac(ey);
  }

  
  const typename Trait::MatType&
  operator()(const typename Trait::Rdx& tx,
	     const typename Trait::Rdy& ty){
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
    return val = ker*phix(kx,tx)*phiy(ky,ty);
  }


};

typedef BIOpKernel<CST,CST_OP,1,P1_1D,P1_1D> CST_1D_P1xP1;
typedef BIOpKernel<CST,CST_OP,1,P0_1D,P0_1D> CST_1D_P0xP0;
typedef BIOpKernel<CST,CST_OP,1,P0_1D,P1_1D> CST_1D_P0xP1;
typedef BIOpKernel<CST,CST_OP,1,P1_1D,P0_1D> CST_1D_P1xP0;
typedef BIOpKernel<CST,CST_OP,2,P1_2D,P1_2D> CST_2D_P1xP1;
typedef BIOpKernel<CST,CST_OP,2,P0_2D,P0_2D> CST_2D_P0xP0;
typedef BIOpKernel<CST,CST_OP,2,P0_2D,P1_2D> CST_2D_P0xP1;
typedef BIOpKernel<CST,CST_OP,2,P1_2D,P0_2D> CST_2D_P1xP0;
typedef BIOpKernel<CST,CST_OP,1,P2_1D,P2_1D> CST_1D_P2xP2;
typedef BIOpKernel<CST,CST_OP,2,P2_2D,P2_2D> CST_2D_P2xP2;


#endif
