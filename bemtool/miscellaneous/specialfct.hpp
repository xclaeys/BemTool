#ifndef SPECIALFCT_HPP
#define SPECIALFCT_HPP
#include <cmath>
#include "coordinates.hpp"

using std::acos;

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/hankel.hpp>
#include <boost/math/complex/acos.hpp>

/*===================
  FONCTIONS DE BESSEL
  ===================*/

inline Cplx BesselJ0(const Real& x){
  return boost::math::cyl_bessel_j(0,x);}

inline Cplx BesselJ1(const Real& x){
  return boost::math::cyl_bessel_j(1,x);}

inline Cplx BesselJ(const int& n, const Real& x){
  return boost::math::cyl_bessel_j(n,x);}


/*================================
  DERIVEES DES FONCTIONS DE BESSEL
  ================================*/

inline Cplx DBesselJ0_Dx(const Real& x){
  return -BesselJ1(x);}

inline Cplx DBesselJ_Dx(const int& n, const Real& x){
  if(n==0){return -BesselJ1(x);} return 0.5*( BesselJ(n-1,x) - BesselJ(n+1,x) );}


/*===================
  FONCTIONS DE HANKEL
  ===================*/
inline Cplx Hankel0 (const Real& x){
  return boost::math::cyl_hankel_1(0,x); }

inline Cplx Hankel1 (const Real& x){
  return boost::math::cyl_hankel_1(1,x); }

inline Cplx Hankel (const int& n, const Real& x){
  return boost::math::cyl_hankel_1(n,x);}


/*================================
  DERIVEES DES FONCTIONS DE HANKEL
  ================================*/

inline Cplx DHankel0_Dx (const Real& x) {
  return -Hankel1(x);}

inline Cplx DHankel_Dx  (const int& n, const Real& x){
  if(n==0){return -Hankel1(x);} return 0.5*( Hankel(n-1,x) - Hankel(n+1,x) ); }


/*===================
  FONCTIONS DE KELVIN
  ===================*/

inline Cplx Kelvin0 (const Real& x){
  return boost::math::cyl_bessel_k(0,x); }

inline Cplx Kelvin1 (const Real& x){
  return boost::math::cyl_bessel_k(1,x); }

inline Cplx Kelvin (const int& n, const Real& x){
  return boost::math::cyl_bessel_k(n,x);}


/*======================
  HARMONIQUES SPHERIQUES
  ======================*/

inline Cplx SphHarmo(const N2& nm, const R3& x){
  return boost::math::spherical_harmonic(nm[0],nm[1],Sph::Theta(x),Sph::Phi(x));}

inline Cplx SphHarmo(const int& n, const int& m, const R3& x){
  return boost::math::spherical_harmonic(n,m,Sph::Theta(x),Sph::Phi(x));}


/*=========================
  GRADIENTS SURFACIQUES DES
  HARMONIQUES SPHERIQUES 
  =========================*/

// ATTENTION: ROUTINE BUGEE POUR LE MOMENT....

inline C3 Grad_SphHarmo(const int& n, const int& m, const R3& y){
  R3 x            = (1./norm2(y))*y;  
  R3 e_z          = R3_(0.,0.,1.);
  R3 e_phi        = vprod(e_z,x);
  e_phi           = (1./norm2(e_phi))*e_phi;
  R3 e_theta      = vprod(e_phi,x);  
  Real cos_theta  = x[2];
  Real sin_theta  = sqrt(1-x[2]*x[2]);
  Real coef       = sqrt( (2*n+1)*( n-abs(m) )/(2*n-1 ) );
  Cplx u_phi, u_theta;
  u_phi    = ( n*cos_theta*SphHarmo(n,m,x)-coef*SphHarmo(n-1,m,x) )/sin_theta;
  u_theta  = iu*m*SphHarmo(n,m,x)/sin_theta; 
  return u_phi*e_phi + u_theta*e_theta;}

inline C3 Grad_SphHarmo(const N2& nm, const R3& y){
  return Grad_SphHarmo(nm[0],nm[1],y);}


/*===========================
  BESSEL ET HANKEL SPHERIQUES
  ===========================*/

inline Cplx SphBesselJ(const int& n, const Real& x){
  return boost::math::sph_bessel(n,x);}

inline Cplx SphHankel(const int& n, const Real& x){
  return boost::math::sph_bessel(n,x) + iu*boost::math::sph_neumann(n,x);}


/*========================================
  DERIVEES DES BESSEL ET HANKEL SPHERIQUES
  ========================================*/

inline Cplx DSphBesselJ_Dx(const int& n, const Real& x){
  return -SphBesselJ(n+1,x)+(n/x)*SphBesselJ(n,x);}

inline Cplx DSphHankel_Dx(const int& n, const Real& x){
  return -SphHankel(n+1,x)+(n/x)*SphHankel(n,x);}




#endif
