#ifndef REFSOL_HPP
#define REFSOL_HPP

#include "../calculus/calculus.hpp"
#include "specialfct.hpp"
using namespace std;





template <int, int, int>
struct AnalyticalSolution;


template<typename Op>
struct RefSol{
  
  static inline Cplx Compute(const int& n, const Real& r=1., const Real& k=1.){
    return AnalyticalSolution<Op::Trait::EqType,Op::Trait::OpType,Op::Trait::Dim>::Compute(n,r,k); }
  
  static inline Cplx Compute(const N2& nm, const Real& r=1., const Real& k=1.){
    return AnalyticalSolution<Op::Trait::EqType,Op::Trait::OpType,Op::Trait::Dim>::Compute(nm,r,k);}  
};


/*==========
  LAPLACE 2D
  ==========*/
template <> struct AnalyticalSolution<LA,SL_OP,2>{
  static inline Cplx  
  Compute(const int& n, const Real& r=1., const Real& k=1.){
    return pi*r/abs(n);} };

template <> struct AnalyticalSolution<LA,DL_OP,2>{
  static inline Cplx
  Compute(const int& n, const Real& r=1., const Real& k=1.){
    return 0.;} };

template <> struct AnalyticalSolution<LA,TDL_OP,2>{
  static inline Cplx
  Compute(const int& n, const Real& r=1., const Real& k=1.){
    return 0.;} };

template <> struct AnalyticalSolution<LA,HS_OP,2>{
  static inline Cplx
  Compute(const int& n, const Real& r=1., const Real& k=1.){
    return pi*abs(n)/r;} };

/*==========
  LAPLACE 3D
  ==========*/

template <> struct AnalyticalSolution<LA,SL_OP,3>{
  static inline Cplx
  Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
    return rho/(2.*nm[0]+1);} };

template <> struct AnalyticalSolution<LA,DL_OP,3>{
  static inline Cplx
  Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
    return 0.5/(2.*nm[0]+1);} };

template <> struct AnalyticalSolution<LA,TDL_OP,3>{
  static inline Cplx
  Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
    return -0.5/(2.*nm[0]+1);} };

template <> struct AnalyticalSolution<LA,HS_OP,3>{
  static inline Cplx
  Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
    return (nm[0]/rho)*(nm[0]+1)/(2.*nm[0]+1);} };

/*============
  HELMHOLTZ 2D
  ============*/

template <> struct AnalyticalSolution<HE,SL_OP,2>{
  static inline Cplx
  Compute(const int& n, const Real& r=1., const Real& k=1.){
    Real r2=r*r;
    Real pi2=pi*pi; 
    return iu*r2*pi2*Hankel(abs(n),k*r)*BesselJ(abs(n),k*r);} };

template <> struct AnalyticalSolution<HE,DL_OP,2>{
  static inline Cplx
  Compute(const int& n, const Real& r=1., const Real& k=1.){
    Real r2 =r*r;
    Real pi2=pi*pi;
    Cplx Res;
    Res  = Hankel(abs(n),k*r)*DBesselJ_Dx(abs(n),k*r);
    Res += DHankel_Dx(abs(n),k*r)*BesselJ(abs(n),k*r);    
    Res *= -0.5*iu*k*r2*pi2;
    return Res;} };

template <> struct AnalyticalSolution<HE,TDL_OP,2>{
  static inline Cplx
  Compute(const int& n, const Real& r=1., const Real& k=1.){
    Real r2 =r*r;
    Real pi2=pi*pi;
    Cplx Res;
    Res  = Hankel(abs(n),k*r)*DBesselJ_Dx(abs(n),k*r);
    Res += DHankel_Dx(abs(n),k*r)*BesselJ(abs(n),k*r);    
    Res *= +0.5*iu*k*r2*pi2;
    return Res;} };

template <> struct AnalyticalSolution<HE,HS_OP,2>{
  static inline Cplx
  Compute(const int& n, const Real& r=1., const Real& k=1.){
    Real r2 =r*r;
    Real pi2=pi*pi;
    Real k2 =k*k; 
    return -iu*k2*r2*pi2*DHankel_Dx(abs(n),k*r)*DBesselJ_Dx(abs(n),k*r);} };

/*============
  HELMHOLTZ 3D
  ============*/

template <> struct AnalyticalSolution<HE,SL_OP,3>{
  static inline Cplx
  Compute(const N2& nm, const Real& rho=1., const Real& k=1.){    
    return iu*k*rho*rho*SphBesselJ(nm[0],k*rho)*SphHankel(nm[0],k*rho);} };

template <> struct AnalyticalSolution<HE,DL_OP,3>{
  static inline Cplx
  Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
    Cplx Res;
    Res  = DSphBesselJ_Dx(nm[0],k*rho)*SphHankel(nm[0],k*rho);
    Res += DSphHankel_Dx(nm[0],k*rho)*SphBesselJ(nm[0],k*rho);
    Res *= -0.5*iu*k*k*rho*rho;
    return Res;} };

template <> struct AnalyticalSolution<HE,TDL_OP,3>{
  static inline Cplx
  Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
    Cplx Res;
    Res  = DSphBesselJ_Dx(nm[0],k*rho)*SphHankel(nm[0],k*rho);
    Res += DSphHankel_Dx(nm[0],k*rho)*SphBesselJ(nm[0],k*rho);
    Res *= 0.5*iu*k*k*rho*rho;
    return Res;} };

template <> struct AnalyticalSolution<HE,HS_OP,3>{
  static inline Cplx
  Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
    Real k3=k*k*k;
    return -iu*k3*rho*rho*DSphBesselJ_Dx(nm[0],k*rho)*DSphHankel_Dx(nm[0],k*rho);} };


/*=======
  MAXWELL
  =======*/

template <> struct AnalyticalSolution<MA,SL_OP,3>{
  static inline Cplx
  Compute(const N2& nm, const Real& rho=1., const Real& k=1.){    
    Cplx res = iu/k;
    res *= SphBesselJ(nm[0],k) +k*DSphBesselJ_Dx(nm[0],k);
    res *= SphHankel (nm[0],k) +k*DSphHankel_Dx (nm[0],k);
    //   res *= k*SphBesselJ(nm[0],k);
    //   res *= k*SphHankel(nm[0],k);
    return res;}
};







#endif
