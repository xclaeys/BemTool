#ifndef HELMHOLTZ_POT_HPP
#define HELMHOLTZ_POT_HPP

#include "../calculus/calculus.hpp"
#include "../miscellaneous/specialfct.hpp"


/*=============================
  SIMPLE COUCHE HELMHOLTZ EN 2D
  =============================*/

template <typename PhiY>
class PotKernel<HE,SL_POT,2,PhiY>{
  
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
  
  PotKernel<HE,SL_POT,2,PhiY>(const typename Trait::MeshY& my,
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
    ker = h*0.25*iu*Hankel0(kappa*r);
    for(int k=0; k<Trait::nb_dof_y; k++){
      mat(0,k) = ker*phiy(k,tj);}
    return mat;}
  
  const Cplx&
  operator()(const R3& x, const typename Trait::Rdy& tj, const int& ky){    
    x_y  = x0_y0-dy*tj;
    r   = norm2(x_y);
    ker = h*0.25*iu*Hankel0(kappa*r);
    return ker *= phiy(ky,tj);
  }
    
};

typedef PotKernel<HE,SL_POT,2,P0_1D> HE_SL_2D_P0;
typedef PotKernel<HE,SL_POT,2,P1_1D> HE_SL_2D_P1;



/*=============================
  SIMPLE COUCHE HELMHOLTZ EN 3D
  =============================*/

template <typename PhiY>
class PotKernel<HE,SL_POT,3,PhiY>{
  
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
  
  PotKernel<HE,SL_POT,3,PhiY>(const typename Trait::MeshY& my,
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
    ker = h*exp(iu*kappa*r)/(4*pi*r);
    for(int k=0; k<Trait::nb_dof_y; k++){
      mat(0,k) = ker*phiy(k,tj);}
    return mat;}
  
  const Cplx&
  operator()(const R3& x, const typename Trait::Rdy& tj, const int& ky){    
    x_y  = x0_y0-dy*tj;
    r   = norm2(x_y);
    ker = h*exp(iu*kappa*r)/(4*pi*r);
    return ker *= phiy(ky,tj);
  }
    
};

typedef PotKernel<HE,SL_POT,3,P0_2D> HE_SL_3D_P0;
typedef PotKernel<HE,SL_POT,3,P1_2D> HE_SL_3D_P1;




/*=============================
  DOUBLE COUCHE HELMHOLTZ EN 2D
  =============================*/

template <typename PhiY>
class PotKernel<HE,DL_POT,2,PhiY>{
  
public: 
  typedef PotKernelTraits<PhiY> Trait;
  
private:
  const typename Trait::MeshY&        meshy;    
        typename Trait::MatType       mat;
        typename Trait::JacY          dy;
                        PhiY          phiy;
  const                 vector<R3>&   normaly;
  const                 Real          kappa;
                        R3            x0_y0,x_y,ny;
                        Real          h,r;
                        Cplx          ker; 

public:
  
  PotKernel<HE,DL_POT,2,PhiY>(const typename Trait::MeshY& my,
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
    ker = -h*(ny,x_y)*(1./r)*0.25*iu*kappa*Hankel1(kappa*r);
    for(int k=0; k<Trait::nb_dof_y; k++){
      mat(0,k) = ker*phiy(k,tj);}
    return mat;
  }
  
  const Cplx&
  operator()(const R3& x, const typename Trait::Rdy& tj, const int& ky){    
    x_y  = x0_y0-dy*tj;
    r   = norm2(x_y);
    ker = -h*(ny,x_y)*(1./r)*0.25*iu*kappa*Hankel1(kappa*r);
    return ker *= phiy(ky,tj);
  }
    
};


typedef PotKernel<HE,DL_POT,2,P0_1D> HE_DL_2D_P0;
typedef PotKernel<HE,DL_POT,2,P1_1D> HE_DL_2D_P1;



/*=============================
  DOUBLE COUCHE HELMHOLTZ EN 3D
  =============================*/

template <typename PhiY>
class PotKernel<HE,DL_POT,3,PhiY>{
  
public: 
  typedef PotKernelTraits<PhiY> Trait;
  
private:
  const typename Trait::MeshY&        meshy;    
        typename Trait::MatType       mat;
        typename Trait::JacY          dy;
                        PhiY          phiy;
  const                 vector<R3>&   normaly;
  const                 Real          kappa;
                        R3            x0_y0,x_y,ny;
                        Real          h,r,r3;
                        Cplx          ker; 

public:
  
  PotKernel<HE,DL_POT,3,PhiY>(const typename Trait::MeshY& my,
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
    ker = h*(ny,x_y)*(iu*kappa*r-1.)*exp(iu*kappa*r)/(4*pi*r3);
    for(int k=0; k<Trait::nb_dof_y; k++){
      mat(0,k) = ker*phiy(k,tj);}
    return mat;}
  
  const Cplx&
  operator()(const R3& x, const typename Trait::Rdy& tj, const int& ky){    
    x_y  = x0_y0-dy*tj;
    r   = norm2(x_y);
    r3  = r*r*r; 
    ker = h*(ny,x_y)*(iu*kappa*r-1.)*exp(iu*kappa*r)/(4*pi*r3);
    return ker *= phiy(ky,tj);
  }
    
};

typedef PotKernel<HE,DL_POT,3,P0_2D> HE_DL_3D_P0;
typedef PotKernel<HE,DL_POT,3,P1_2D> HE_DL_3D_P1;




#endif
