// //===================================================================
// //
// //  Copyright 2017 Xavier Claeys
// //
// //  bemtool is free software: you can redistribute it and/or modify
// //  it under the terms of the GNU General Public License as published by
// //  the Free Software Foundation, either version 3 of the License, or
// //  (at your option) any later version.
// //
// //  bemtool is distributed in the hope that it will be useful,
// //  but WITHOUT ANY WARRANTY; without even the implied warranty of
// //  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// //  GNU General Public License for more details.
// //
// //  You should have received a copy of the GNU General Public License
// //  along with bemtool.  If not, see <http://www.gnu.org/licenses/>
// //
// //====================================================================
// #ifndef BEMTOOL_MISC_REFSOL_HPP
// #define BEMTOOL_MISC_REFSOL_HPP
//
// #include "../operator/operator.hpp"
//
//
//
// namespace bemtool {
//
//
// template <int, int, int>
// struct AnalyticalSolution;
//
//
// template<typename Op>
// struct RefSol{
//
//   static inline Cplx Compute(Real r, Real theta, const int& n, const Real& r=1., const Real& k=1.){
//     return AnalyticalSolution<Op::Trait::EqType,Op::Trait::OpType,Op::Trait::Dim>::Compute(r,theta,n,r,k); }
//
//   static inline Cplx Compute(Real r, Real theta, const N2& nm, const Real& r=1., const Real& k=1.){
//     return AnalyticalSolution<Op::Trait::EqType,Op::Trait::OpType,Op::Trait::Dim>::Compute(r,theta,nm,r,k);}
// };
//
//
// /*==========
//   LAPLACE 2D
//   ==========*/
// template <> struct AnalyticalSolution<LA,SL_OP,2>{
//   static inline Cplx
//   Compute(const int& n, const Real& r=1., const Real& k=1.){
//     return pi*r/abs(n);} };
//
// template <> struct AnalyticalSolution<LA,DL_OP,2>{
//   static inline Cplx
//   Compute(const int& n, const Real& r=1., const Real& k=1.){
//     return 0.;} };
//
// template <> struct AnalyticalSolution<LA,TDL_OP,2>{
//   static inline Cplx
//   Compute(const int& n, const Real& r=1., const Real& k=1.){
//     return 0.;} };
//
// template <> struct AnalyticalSolution<LA,HS_OP,2>{
//   static inline Cplx
//   Compute(const int& n, const Real& r=1., const Real& k=1.){
//     return pi*abs(n)/r;} };
//
// /*==========
//   LAPLACE 3D
//   ==========*/
//
// template <> struct AnalyticalSolution<LA,SL_OP,3>{
//   static inline Cplx
//   Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
//     return rho/(2.*nm[0]+1);} };
//
// template <> struct AnalyticalSolution<LA,DL_OP,3>{
//   static inline Cplx
//   Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
//     return 0.5/(2.*nm[0]+1);} };
//
// template <> struct AnalyticalSolution<LA,TDL_OP,3>{
//   static inline Cplx
//   Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
//     return -0.5/(2.*nm[0]+1);} };
//
// template <> struct AnalyticalSolution<LA,HS_OP,3>{
//   static inline Cplx
//   Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
//     return (nm[0]/rho)*(nm[0]+1)/(2.*nm[0]+1);} };
//
// /*============
//   HELMHOLTZ 2D
//   ============*/
//
// template <> struct AnalyticalSolution<HE,SL_OP,2>{
//   static inline Cplx
//   Compute(const int& n, const Real& r=1., const Real& k=1.){
//     Real r2=r*r;
//     Real pi2=pi*pi;
//     return iu*r2*pi2*Hankel(abs(n),k*r)*BesselJ(abs(n),k*r);} };
//
// template <> struct AnalyticalSolution<HE,DL_OP,2>{
//   static inline Cplx
//   Compute(const int& n, const Real& r=1., const Real& k=1.){
//     Real r2 =r*r;
//     Real pi2=pi*pi;
//     Cplx Res;
//     Res  = Hankel(abs(n),k*r)*DBesselJ_Dx(abs(n),k*r);
//     Res += DHankel_Dx(abs(n),k*r)*BesselJ(abs(n),k*r);
//     Res *= -0.5*iu*k*r2*pi2;
//     return Res;} };
//
// template <> struct AnalyticalSolution<HE,TDL_OP,2>{
//   static inline Cplx
//   Compute(const int& n, const Real& r=1., const Real& k=1.){
//     Real r2 =r*r;
//     Real pi2=pi*pi;
//     Cplx Res;
//     Res  = Hankel(abs(n),k*r)*DBesselJ_Dx(abs(n),k*r);
//     Res += DHankel_Dx(abs(n),k*r)*BesselJ(abs(n),k*r);
//     Res *= +0.5*iu*k*r2*pi2;
//     return Res;} };
//
// template <> struct AnalyticalSolution<HE,HS_OP,2>{
//   static inline Cplx
//   Compute(const int& n, const Real& r=1., const Real& k=1.){
//     Real r2 =r*r;
//     Real pi2=pi*pi;
//     Real k2 =k*k;
//     return -iu*k2*r2*pi2*DHankel_Dx(abs(n),k*r)*DBesselJ_Dx(abs(n),k*r);} };
//
// /*============
//   HELMHOLTZ 3D
//   ============*/
//
// template <> struct AnalyticalSolution<HE,SL_OP,3>{
//   static inline Cplx
//   Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
//     return iu*k*rho*rho*SphBesselJ(nm[0],k*rho)*SphHankel(nm[0],k*rho);} };
//
// template <> struct AnalyticalSolution<HE,DL_OP,3>{
//   static inline Cplx
//   Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
//     Cplx Res;
//     Res  = DSphBesselJ_Dx(nm[0],k*rho)*SphHankel(nm[0],k*rho);
//     Res += DSphHankel_Dx(nm[0],k*rho)*SphBesselJ(nm[0],k*rho);
//     Res *= -0.5*iu*k*k*rho*rho;
//     return Res;} };
//
// template <> struct AnalyticalSolution<HE,TDL_OP,3>{
//   static inline Cplx
//   Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
//     Cplx Res;
//     Res  = DSphBesselJ_Dx(nm[0],k*rho)*SphHankel(nm[0],k*rho);
//     Res += DSphHankel_Dx(nm[0],k*rho)*SphBesselJ(nm[0],k*rho);
//     Res *= 0.5*iu*k*k*rho*rho;
//     return Res;} };
//
// template <> struct AnalyticalSolution<HE,HS_OP,3>{
//   static inline Cplx
//   Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
//     Real k3=k*k*k;
//     return -iu*k3*rho*rho*DSphBesselJ_Dx(nm[0],k*rho)*DSphHankel_Dx(nm[0],k*rho);} };
//
// /*============
//   Modified Helmholtz 2D
//   ============*/
//
// template <> struct AnalyticalSolution<YU,SL_OP,2>{
//   static inline Cplx
//   Compute(const int& n, const Real& r=1., const Real& k=1.){
//     Real r2=r*r;
//     return 2*r2*pi*Modified_BesselK(abs(n),k*r)*Modified_BesselI(abs(n),k*r);} };
//
// template <> struct AnalyticalSolution<YU,DL_OP,2>{
//   static inline Cplx
//   Compute(const int& n, const Real& r=1., const Real& k=1.){
//     Real r2 =r*r;
//     Cplx Res;
//     Res  = DModified_BesselK_Dx(abs(n),k*r)*Modified_BesselI(abs(n),k*r);
//     Res += Modified_BesselK(abs(n),k*r)*DModified_BesselI_Dx(abs(n),k*r);
//     Res *= -k*r2*pi;
//     return Res;} };
//
// template <> struct AnalyticalSolution<YU,TDL_OP,2>{
//   static inline Cplx
//   Compute(const int& n, const Real& r=1., const Real& k=1.){
//     Real r2 =r*r;
//     Cplx Res;
//     Res  = DModified_BesselK_Dx(abs(n),k*r)*Modified_BesselI(abs(n),k*r);
//     Res += Modified_BesselK(abs(n),k*r)*DModified_BesselI_Dx(abs(n),k*r);
//     Res *= k*r2*pi;
//     return Res;} };
//
// template <> struct AnalyticalSolution<YU,HS_OP,2>{
//   static inline Cplx
//   Compute(const int& n, const Real& r=1., const Real& k=1.){
//     Real r2 =r*r;
//     Real k2 =k*k;
//     return -2*k2*r2*pi*DModified_BesselI_Dx(abs(n),k*r)*DModified_BesselK_Dx(abs(n),k*r);} };
//
// /*============
//   Modified Helmholtz 3D
//   ============*/
//
// template <> struct AnalyticalSolution<YU,SL_OP,3>{
//   static inline Cplx
//   Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
//     return iu*k*rho*rho*SphBesselJ(nm[0],k*rho)*SphHankel(nm[0],k*rho);} };
//
// template <> struct AnalyticalSolution<YU,DL_OP,3>{
//   static inline Cplx
//   Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
//     Cplx Res;
//     Res  = DSphBesselJ_Dx(nm[0],k*rho)*SphHankel(nm[0],k*rho);
//     Res += DSphHankel_Dx(nm[0],k*rho)*SphBesselJ(nm[0],k*rho);
//     Res *= -0.5*iu*k*k*rho*rho;
//     return Res;} };
//
// template <> struct AnalyticalSolution<YU,TDL_OP,3>{
//   static inline Cplx
//   Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
//     Cplx Res;
//     Res  = DSphBesselJ_Dx(nm[0],k*rho)*SphHankel(nm[0],k*rho);
//     Res += DSphHankel_Dx(nm[0],k*rho)*SphBesselJ(nm[0],k*rho);
//     Res *= 0.5*iu*k*k*rho*rho;
//     return Res;} };
//
// template <> struct AnalyticalSolution<YU,HS_OP,3>{
//   static inline Cplx
//   Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
//     Real k3=k*k*k;
//     return -iu*k3*rho*rho*DSphBesselJ_Dx(nm[0],k*rho)*DSphHankel_Dx(nm[0],k*rho);} };
//
//
// /*=======
//   MAXWELL
//   =======*/
//
// template <> struct AnalyticalSolution<MA,SL_OP,3>{
//   static inline Cplx
//   Compute(const N2& nm, const Real& rho=1., const Real& k=1.){
//     Cplx res = iu/k;
//     res *= SphBesselJ(nm[0],k) +k*DSphBesselJ_Dx(nm[0],k);
//     res *= SphHankel (nm[0],k) +k*DSphHankel_Dx (nm[0],k);
//     //   res *= k*SphBesselJ(nm[0],k);
//     //   res *= k*SphHankel(nm[0],k);
//     return res;}
// };
//
//
// }
//
//
//
//
// #endif
