//===================================================================
//
//  Copyright 2017 Xavier Claeys
//
//  bemtool is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  bemtool is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//  (see ../LICENSE.txt)
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with bemtool.  If not, see <http://www.gnu.org/licenses/>
//
//====================================================================
#ifndef BEMTOOL_CALCULUS_EXPRESSION_HPP
#define BEMTOOL_CALCULUS_EXPRESSION_HPP


#include <complex>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <vector>

namespace bemtool{

//==========================//
//      Types de base       //
//==========================//

typedef double        Real;
typedef std::complex<Real> Cplx;
const   Cplx iu(0.,1.);
const   Real pi = 3.14159265358979;

inline Cplx operator+ (const int& n,  const Cplx& z){ return (double)n + z;}
inline Cplx operator+ (const Cplx& z, const int& n) { return z + (double)n;}

inline Cplx operator- (const int& n,  const Cplx& z){ return (double)n - z;}
inline Cplx operator- (const Cplx& z, const int& n) { return z - (double)n;}

inline Cplx operator* (const int& n,  const Cplx& z){ return (double)n * z;}
inline Cplx operator* (const Cplx& z, const int& n) { return z * (double)n;}

inline Real conj(const Real& x){return x;}
inline int  conj(const int& x) {return x;}

template <class U> struct isbase                    {static const bool test = false; };
template <>        struct isbase<bool>              {static const bool test = true;  };
template <>        struct isbase<int>               {static const bool test = true;  };
template <>        struct isbase<double>            {static const bool test = true;  };
template <>        struct isbase<float>             {static const bool test = true;  };
template <class U> struct isbase<std::complex<U> >  {static const bool test = true;  };

// Operateurs d'acces
template <class r_t, bool test = isbase<r_t>::test > struct access_ {
    static inline       typename r_t::v_t& op(      r_t& r_, const int& j)              {return r_[j];}
    static inline const typename r_t::v_t& op(const r_t& r_, const int& j)              {return r_[j];}
    static inline       typename r_t::v_t& op(      r_t& r_, const int& j, const int& k){return r_(j,k);}
    static inline const typename r_t::v_t& op(const r_t& r_, const int& j, const int& k){return r_(j,k);}
};

template <class r_t> struct access_<r_t,true> {
    static inline       r_t& op(      r_t& r_, const int& j)              {return r_;}
    static inline const r_t& op(const r_t& r_, const int& j)              {return r_;}
    static inline       r_t& op(      r_t& r_, const int& j, const int& k){return r_;}
    static inline const r_t& op(const r_t& r_, const int& j, const int& k){return r_;}
};

// Type des entrees
template <class r_t, bool test = isbase<r_t>::test> struct entry           {typedef typename r_t::v_t type; };
template <class r_t>                                struct entry<r_t,true> {typedef r_t               type; };

// Gestion des operations entre type
template <class l_t, class r_t> struct res                          {typedef l_t               type;};
template <>                     struct res< int, double >           {typedef double            type;};
template <>                     struct res< double, int >           {typedef double            type;};
template <class l_t>            struct res< l_t,std::complex<l_t> > {typedef std::complex<l_t> type;};
template <class r_t>            struct res< std::complex<r_t>,r_t > {typedef std::complex<r_t> type;};
template <class l_t>            struct res< std::complex<l_t>,int > {typedef std::complex<l_t> type;};
template <class r_t>            struct res< int,std::complex<r_t> > {typedef std::complex<r_t> type;};


template <class l_t, class r_t> struct resop{
  typedef typename res<typename entry<l_t>::type,typename entry<r_t>::type>::type  type; };


template <class l_t, class r_t, bool l_test = isbase<l_t>::test, bool r_test = isbase<r_t>::test >
struct find_nb_rows{ static const int nr = l_t::nr;};
template <class l_t, class r_t >
struct find_nb_rows<l_t, r_t, true, false>{ static const int nr = r_t::nr;};



//==========================//
//        Boucles           //
//==========================//

//---------------------------//
template <class r_t, int D = r_t::dim> struct ascending_loop{
    static inline void apply(r_t& r_){
        r_[D-1]=D-1; ascending_loop<r_t,D-1>::apply(r_);} };

template <class r_t> struct ascending_loop<r_t,1>{
    static inline void apply(r_t& r_){r_[0]=0;} };

//---------------------------//
template <class r_t, int D = r_t::dim> struct decrement_loop{
    static inline void apply(r_t& r_){
        r_[D-1]--; decrement_loop<r_t,D-1>::apply(r_);} };

template <class r_t> struct decrement_loop<r_t,1>{
    static inline void apply(r_t& r_){r_[0]--;} };

//---------------------------//
template <class r_t, int D = r_t::dim> struct increment_loop{
    static inline void apply(r_t& r_){
        r_[D-1]++; increment_loop<r_t,D-1>::apply(r_);} };

template <class r_t> struct increment_loop<r_t,1>{
    static inline void apply(r_t& r_){r_[0]++;} };

//---------------------------//
template <class r_t, int D = r_t::dim> struct istream_loop{
    static inline void apply(std::istream& is, r_t& r_){
        is >> r_[r_t::dim-D]; istream_loop<r_t,D-1>::apply(is,r_);} };

template <class r_t> struct istream_loop<r_t,1>{
    static inline void apply(std::istream& is, r_t& r_){
        is >> r_[r_t::dim-1];} };

//---------------------------//
template <class r_t, int d = r_t::dim> struct ostream_loop{
    static inline void apply(std::ostream& os, const r_t& r_){
        os<<r_[r_t::dim-d] << "\t"; ostream_loop<r_t,d-1>::apply(os,r_);} };

template <class r_t> struct ostream_loop<r_t,1>{
    static inline void apply(std::ostream& os, const r_t& r_){
        os<<r_[r_t::dim-1];} };

//---------------------------//
template <class l_t, class r_t, int D = l_t::dim> struct assign_loop{
    static void apply(l_t& l_, const r_t& r_){
        l_[D-1] = access_<r_t>::op(r_,D-1); assign_loop<l_t,r_t,D-1>::apply(l_,r_);} };

template <class l_t, class r_t> struct assign_loop<l_t,r_t,1>{
    static void apply(l_t& l_, const r_t& r_){ l_[0] = access_<r_t>::op(r_,0);} };

//---------------------------//
template <class l_t, class r_t, int D = l_t::dim> struct plus_assign_loop{
    static void apply(l_t& l_, const r_t& r_){
        l_[D-1] += access_<r_t>::op(r_,D-1); plus_assign_loop<l_t,r_t,D-1>::apply(l_,r_);} };

template <class l_t, class r_t> struct plus_assign_loop<l_t,r_t,1>{
    static void apply(l_t& l_, const r_t& r_){ l_[0] += access_<r_t>::op(r_,0);} };

//---------------------------//
template <class l_t, int d = l_t::dim> struct construct_loop{
    static void apply(l_t& l_){
        l_[d-1] = typename l_t::v_t();
        construct_loop<l_t,d-1>::apply(l_);} };

template <class l_t> struct construct_loop<l_t,1>{
    static void apply(l_t& l_){
        l_[0] = typename l_t::v_t();} };

//---------------------------//
template <class l_t, int j = l_t::nr, int k = l_t::nc> struct construct_dbloop{
    static inline void apply(l_t& l_){
        l_(j-1,k-1) = typename l_t::v_t(); construct_dbloop<l_t,j,k-1>::apply(l_); }};

template <class l_t, int j> struct construct_dbloop<l_t,j,1>{
    static inline void apply(l_t& l_){
        l_(j-1,0) = typename l_t::v_t(); construct_dbloop<l_t,j-1,l_t::nc>::apply(l_); }};

template <class l_t> struct construct_dbloop<l_t,1,1>{
    static inline void apply(l_t& l_){
        l_(0,0) = typename l_t::v_t(); }};

//---------------------------//
template <class l_t, class r_t, int j = l_t::nr, int k = l_t::nc> struct assign_dbloop{
    static inline void apply(l_t& l_, const r_t& r_){
        l_(j-1,k-1) = access_<r_t>::op(r_,j-1,k-1); assign_dbloop<l_t,r_t,j,k-1>::apply(l_,r_); }};

template <class l_t, class r_t, int j> struct assign_dbloop<l_t,r_t,j,1>{
    static inline void apply(l_t& l_, const r_t& r_){
        l_(j-1,0) = access_<r_t>::op(r_,j-1,0); assign_dbloop<l_t,r_t,j-1,l_t::nc>::apply(l_,r_); }};

template <class l_t, class r_t> struct assign_dbloop<l_t,r_t,1,1>{
    static inline void apply(l_t& l_, const r_t& r_){
        l_(0,0) = access_<r_t>::op(r_,0,0); }};

//---------------------------//
template <class l_t, class r_t, int j = l_t::nr, int k = l_t::nc> struct plus_assign_dbloop{
    static inline void apply(l_t& l_, const r_t& r_){
        l_(j-1,k-1) += access_<r_t>::op(r_,j-1,k-1); plus_assign_dbloop<l_t,r_t,j,k-1>::apply(l_,r_); }};

template <class l_t, class r_t, int j> struct plus_assign_dbloop<l_t,r_t,j,1>{
    static inline void apply(l_t& l_, const r_t& r_){
        l_(j-1,0) += access_<r_t>::op(r_,j-1,0); plus_assign_dbloop<l_t,r_t,j-1,l_t::nc>::apply(l_,r_); }};

template <class l_t, class r_t> struct plus_assign_dbloop<l_t,r_t,1,1>{
    static inline void apply(l_t& l_, const r_t& r_){
        l_(0,0) += access_<r_t>::op(r_,0,0); }};

//---------------------------//
template <class r_t, int nr = r_t::nr, int nc = r_t::nc> struct ostream_dbloop{
    static inline void apply(std::ostream& os, const r_t& r_){
        os<<r_(r_t::nr-nr,r_t::nc-nc) << "\t"; ostream_dbloop<r_t,nr,nc-1>::apply(os,r_);} };

template <class r_t, int nr> struct ostream_dbloop<r_t,nr,1>{
    static inline void apply(std::ostream& os, const r_t& r_){
        os<<r_(r_t::nr-nr,r_t::nc-1) << std::endl; ostream_dbloop<r_t,nr-1,r_t::nc>::apply(os,r_);} };

template <class r_t> struct ostream_dbloop<r_t,1,1>{
    static inline void apply(std::ostream& os, const r_t& r_){
        os<<r_(r_t::nr-1,r_t::nc-1);} };

//---------------------------//
template <class v_t, class l_t, class r_t, int d = l_t::nc> struct mat_mult{

    static v_t apply(const l_t& l_, const r_t& r_, const int& j){
        return access_<l_t>::op(l_,j,d-1)*access_<r_t>::op(r_,d-1) + mat_mult<v_t,l_t,r_t,d-1>::apply(l_,r_,j);}

    static v_t apply(const l_t& l_, const r_t& r_, const int& j, const int& k){
        return access_<l_t>::op(l_,j,d-1)*access_<r_t>::op(r_,d-1,k) + mat_mult<v_t,l_t,r_t,d-1>::apply(l_,r_,j);}

};

template <class v_t, class l_t, class r_t> struct mat_mult<v_t,l_t,r_t,1>{

    static v_t apply(const l_t& l_, const r_t& r_, const int& j){
        return access_<l_t>::op(l_,j,0)*access_<r_t>::op(r_,0);}

    static v_t apply(const l_t& l_, const r_t& r_, const int& j, const int& k){
        return access_<l_t>::op(l_,j,0)*access_<r_t>::op(r_,0,k);}

};

//---------------------------//
template <class l_t,class r_t, int k> struct dprod{
    static typename resop<l_t,r_t>::type apply(const l_t& l_, const r_t& r_){
        return l_[k-1]*conj(r_[k-1]) +  dprod<l_t,r_t,k-1>::apply(l_,r_);}};

template <class l_t, class r_t> struct dprod<l_t,r_t,1>{
    static typename resop<l_t,r_t>::type apply(const l_t& l_, const r_t& r_){
        return l_[0]*conj(r_[0]);}};








//==========================//
//       Operations         //
//==========================//


//%%%%%%%%%%//
// Addition //
//%%%%%%%%%%//

template <class lhs_t, class rhs_t>
struct pp{

    typedef lhs_t                          l_t;
    typedef rhs_t                          r_t;
    typedef typename resop<l_t,r_t>::type  v_t;


    static inline v_t apply(const l_t& l_, const r_t& r_, const int& j){
        return access_<l_t>::op(l_,j) + access_<r_t>::op(r_,j);}

    static inline v_t apply(const l_t& l_, const r_t& r_, const int& j, const int& k){
        return access_<l_t>::op(l_,j,k) + access_<r_t>::op(r_,j,k);}

};


//%%%%%%%%%%%%%%//
// Soustraction //
//%%%%%%%%%%%%%%//

template <class lhs_t, class rhs_t>
struct mm{

    typedef lhs_t                          l_t;
    typedef rhs_t                          r_t;
    typedef typename resop<l_t,r_t>::type  v_t;

    static inline v_t apply(const l_t& l_, const r_t& r_, const int& j){
        return access_<l_t>::op(l_,j) - access_<r_t>::op(r_,j);}

    static inline v_t apply(const l_t& l_, const r_t& r_, const int& j, const int& k){
        return access_<l_t>::op(l_,j,k) - access_<r_t>::op(r_,j,k);}
};


//%%%%%%%%%%%%%%%%//
// Multiplication //
//%%%%%%%%%%%%%%%%//

template <class lhs_t, class rhs_t, bool l_test = isbase<lhs_t>::test, bool r_test = isbase<rhs_t>::test>
struct tt{

    typedef lhs_t                            l_t;
    typedef rhs_t                            r_t;
    typedef typename resop<l_t,r_t>::type    v_t;

    static inline v_t apply(const l_t& l_, const r_t& r_, const int& j){
        return mat_mult<v_t,l_t,r_t>::apply(l_,r_,j);}

    static inline v_t apply(const l_t& l_, const r_t& r_, const int& j, const int& k){
        return mat_mult<v_t,l_t,r_t>::apply(l_,r_,j,k);}

};


template <class lhs_t, class rhs_t>
struct tt<lhs_t,rhs_t, true, false>{

    typedef lhs_t                            l_t;
    typedef rhs_t                            r_t;
    typedef typename resop<l_t,r_t>::type    v_t;

    static inline v_t apply(const l_t& l_, const r_t& r_, const int& j){
        return access_<l_t>::op(l_,j)*access_<r_t>::op(r_,j);}

    static inline v_t apply(const l_t& l_, const r_t& r_, const int& j, const int& k){
        return access_<l_t>::op(l_,j,k)*access_<r_t>::op(r_,j,k);}

};

template <class lhs_t, class rhs_t>
struct tt<lhs_t,rhs_t, false, true>{

    typedef lhs_t                            l_t;
    typedef rhs_t                            r_t;
    typedef typename resop<l_t,r_t>::type    v_t;

    static inline v_t apply(const l_t& l_, const r_t& r_, const int& j){
        return access_<l_t>::op(l_,j)*access_<r_t>::op(r_,j);}

    static inline v_t apply(const l_t& l_, const r_t& r_, const int& j, const int& k){
        return access_<l_t>::op(l_,j,k)*access_<r_t>::op(r_,j,k);}

};



//%%%%%%%%%%%%%%%%//
//  Comparaisons  //
//%%%%%%%%%%%%%%%%//

// template <class l_t, class r_t> bool operator==(const l_t& l_, const r_t& r_){
//  for(int j=0; j<l_t::dim; j++){ if( access_<l_t>::op(l_,j)!=access_<r_t>::op(r_,j) ){return false;}}
//  return true;}

// template <class l_t, class r_t> bool operator!=(const l_t& l_, const r_t& r_){
//  return !(l_==r_);}



//==========================//
//       Expressions        //
//==========================//

template <class op>
class xpr{

public:

    typedef typename op::l_t l_t;
    typedef typename op::r_t r_t;
    typedef typename op::v_t v_t;
    typedef xpr<op>       this_t;

    static const int nr = find_nb_rows<l_t,r_t>::nr;
    static const int dim = nr;

private:

    const l_t& l_;
    const r_t& r_;

public:

  xpr(const l_t& l0, const r_t& r0): l_(l0), r_(r0) {};

  v_t operator[](const int& j){assert(j>=0 && j<nr);  return op::apply(l_,r_,j); }
  const v_t operator[](const int& j) const {assert(j>=0 && j<nr); return op::apply(l_,r_,j); }

  v_t operator()(const int& j, const int& k){assert(j>=0 && j<nr); return op::apply(l_,r_,j,k); }
  const v_t operator()(const int& j, const int& k) const {assert(j>=0 && j<nr); return op::apply(l_,r_,j,k); }

  //==== Addition
  template <class r_t> xpr< pp<this_t,r_t> > operator+(const r_t& r_) const {return xpr< pp<this_t,r_t> >(*this,r_);}
  xpr< pp<this_t, int>  > operator+(const int&  r_) const {return xpr< pp<this_t, int>  >(*this,r_);}
  xpr< pp<this_t, Real> > operator+(const Real& r_) const {return xpr< pp<this_t, Real> >(*this,r_);}
  xpr< pp<this_t, Cplx> > operator+(const Cplx& r_) const {return xpr< pp<this_t, Cplx> >(*this,r_);}

  inline friend xpr< pp<int,this_t>  > operator+(const int&  l_, const this_t& r_){return xpr< pp<int,this_t>  >(l_,r_);}
  inline friend xpr< pp<Real,this_t> > operator+(const Real& l_, const this_t& r_){return xpr< pp<Real,this_t> >(l_,r_);}
  inline friend xpr< pp<Cplx,this_t> > operator+(const Cplx& l_, const this_t& r_){return xpr< pp<Cplx,this_t> >(l_,r_);}

  //==== Soustraction
  template <class r_t> xpr< mm<this_t,r_t> > operator-(const r_t& r_) const {return xpr< mm<this_t,r_t> >(*this,r_);}
  xpr< mm<this_t, int>  > operator-(const int&  r_) const {return xpr< mm<this_t, int>  >(*this,r_);}
  xpr< mm<this_t, Real> > operator-(const Real& r_) const {return xpr< mm<this_t, Real> >(*this,r_);}
  xpr< mm<this_t, Cplx> > operator-(const Cplx& r_) const {return xpr< mm<this_t, Cplx> >(*this,r_);}

  inline friend xpr< mm<int,this_t>  > operator-(const int&  l_, const this_t& r_){return xpr< mm<int,this_t>  >(l_,r_);}
  inline friend xpr< mm<Real,this_t> > operator-(const Real& l_, const this_t& r_){return xpr< mm<Real,this_t> >(l_,r_);}
  inline friend xpr< mm<Cplx,this_t> > operator-(const Cplx& l_, const this_t& r_){return xpr< mm<Cplx,this_t> >(l_,r_);}

  //==== Multiplication
  template <class r_t> xpr< tt<this_t,r_t> > operator*(const r_t& r_) const {return xpr< tt<this_t,r_t> >(*this,r_);}
  xpr< tt<this_t, int>  > operator*(const int&  r_) const {return xpr< tt<this_t, int>  >(*this,r_);}
  xpr< tt<this_t, Real> > operator*(const Real& r_) const {return xpr< tt<this_t, Real> >(*this,r_);}
  xpr< tt<this_t, Cplx> > operator*(const Cplx& r_) const {return xpr< tt<this_t, Cplx> >(*this,r_);}

  inline friend xpr< tt<this_t,int>  > operator*(const int&  l_, const this_t& r_){return xpr< tt<this_t,int>  >(r_,l_);}
  inline friend xpr< tt<this_t,Real> > operator*(const Real& l_, const this_t& r_){return xpr< tt<this_t,Real> >(r_,l_);}
  inline friend xpr< tt<this_t,Cplx> > operator*(const Cplx& l_, const this_t& r_){return xpr< tt<this_t,Cplx> >(r_,l_);}

  //===== Produit scalaire
  template <class r_t>typename resop<this_t,r_t>::type operator,(const r_t& r_) const {
    assert(nr==r_t::nr); return dprod<this_t,r_t,nr>::apply(*this,r_); }

};

template <class t> struct access_<xpr<t>, false> {
  static inline typename xpr<t>::v_t op(xpr<t>& r_, const int& j){return r_[j];}
  static inline const typename xpr<t>::v_t op(const xpr<t>& r_, const int& j){return r_[j];}

  static inline typename xpr<t>::v_t op(xpr<t>& r_, const int& j, const int& k){return r_(j,k);}
  static inline const typename xpr<t>::v_t op(const xpr<t>& r_, const int& j, const int& k){return r_(j,k);}
};


} // namespace bemtool

#endif
