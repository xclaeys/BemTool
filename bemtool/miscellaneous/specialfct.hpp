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
#ifndef BEMTOOL_MISC_SPECIALFCT_HPP
#define BEMTOOL_MISC_SPECIALFCT_HPP
#include <cmath>
#include "coordinates.hpp"



#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/hankel.hpp>
#include <boost/math/complex/acos.hpp>


namespace bemtool {

/*===================
  BESSEL FUNCTIONS
  ===================*/

inline Cplx BesselJ0(const Real& x){
  return boost::math::cyl_bessel_j(0,x);}

inline Cplx BesselJ1(const Real& x){
  return boost::math::cyl_bessel_j(1,x);}

inline Cplx BesselJ(const int& n, const Real& x){
  return boost::math::cyl_bessel_j(n,x);}


/*================================
  DERIVATIVES OF BESSEL FUNCTIONS
  ================================*/

inline Cplx DBesselJ0_Dx(const Real& x){
  return -BesselJ1(x);}

inline Cplx DBesselJ_Dx(const int& n, const Real& x){
  if(n==0){return -BesselJ1(x);} return 0.5*( BesselJ(n-1,x) - BesselJ(n+1,x) );}


/*===================
  HANKEL FUNCTIONS
  ===================*/
inline Cplx Hankel0 (const Real& x){
  return boost::math::cyl_hankel_1(0,x); }

inline Cplx Hankel1 (const Real& x){
  return boost::math::cyl_hankel_1(1,x); }

inline Cplx Hankel (const int& n, const Real& x){
  return boost::math::cyl_hankel_1(n,x);}


/*================================
  DERIVATIVES OF HANKEL FUNCTIONS
  ================================*/

inline Cplx DHankel0_Dx (const Real& x) {
  return -Hankel1(x);}

inline Cplx DHankel_Dx  (const int& n, const Real& x){
  if(n==0){return -Hankel1(x);} return 0.5*( Hankel(n-1,x) - Hankel(n+1,x) ); }


/*===================
  MODIFIED BESSEL FUNCTIONS
  ===================*/

inline Cplx Modified_BesselI0 (const Real& x){
  return boost::math::cyl_bessel_i(0,x); }

inline Cplx Modified_BesselI1 (const Real& x){
  return boost::math::cyl_bessel_i(1,x); }

inline Cplx Modified_BesselI (const int& n, const Real& x){
  return boost::math::cyl_bessel_i(n,x);}

inline Cplx Modified_BesselK0 (const Real& x){
  return boost::math::cyl_bessel_k(0,x); }

inline Cplx Modified_BesselK1 (const Real& x){
  return boost::math::cyl_bessel_k(1,x); }

inline Cplx Modified_BesselK (const int& n, const Real& x){
  return boost::math::cyl_bessel_k(n,x);}

/*===================
  DERIVATIVES OF MODIFIED BESSEL FUNCTIONS
  ===================*/

inline Cplx DModified_BesselI_Dx (const int& n, const Real& x){
  if(n==0){return Modified_BesselI1(x);} return 0.5*( Modified_BesselI(n-1,x) + Modified_BesselI(n+1,x) );}


inline Cplx DModified_BesselK_Dx (const int& n, const Real& x){
  if(n==0){return -Modified_BesselK1(x);} return -0.5*( Modified_BesselK(n-1,x) + Modified_BesselK(n+1,x) );}

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

/*===========================
  BESSEL SPHERIQUES MODIFIEES
  ===========================*/

inline Cplx ModifiedSphBesseli(const int& n, const Real& x){
  return sqrt(pi/(2*x))*boost::math::cyl_bessel_i(n+0.5,x);}

inline Cplx ModifiedSphBesselk(const int& n, const Real& x){
  return sqrt(pi/(2*x))*boost::math::cyl_bessel_k(n+0.5,x);}


/*========================================
  DERIVEES DES BESSEL SPHERIQUES MODIFIEES
  ========================================*/

inline Cplx DModified_Besseli_Dx (const int& n, const Real& x){
  return ModifiedSphBesseli(n+1,x)+(n/x)*ModifiedSphBesseli(n,x);}


inline Cplx DModified_Besselk_Dx (const int& n, const Real& x){
  return -ModifiedSphBesselk(n+1,x)+(n/x)*ModifiedSphBesselk(n,x);}

}


#endif
