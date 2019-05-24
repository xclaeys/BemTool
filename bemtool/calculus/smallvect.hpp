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
#ifndef BEMTOOL_CALCULUS_SMALLVECT_HPP
#define BEMTOOL_CALCULUS_SMALLVECT_HPP
#include <cassert>
#include "expression.hpp"


namespace bemtool{
//==========================//
//        Sub_Array         //
//==========================//

template <class a_t, class i_t>
class subarray{

public:
    static const int dim  =  i_t::dim;
    static const int nr   =  i_t::dim;
    typedef typename a_t::v_t     v_t;
    typedef subarray<a_t,i_t>  this_t;

private:
    a_t&         a_;
    const i_t&   i_;

public:

  subarray(a_t& a0, const i_t& i0): a_(a0), i_(i0) {};

  template <class r_t> void operator=(const r_t& r_){
    assign_loop<this_t,r_t>::apply(*this,r_);}

  template <class r_t> void operator+=(const r_t& r_){
    plus_assign_loop<this_t,r_t>::apply(*this,r_);}

  v_t& operator[](const int& j){assert(j>=0 && j<nr); return a_[i_[j]];}
  const v_t& operator[](const int& j) const {assert(j>=0 && j<nr); return a_[i_[j]];}

  template <class r_t> bool operator==(const r_t& r_){
    for(int j=0; j<dim; j++){if((*this)[j]!=r_[j]){return false;}}
    return true;}

  friend std::ostream& operator<<(std::ostream& os, const this_t& ar){
    ostream_loop<this_t>::apply(os,ar); return os;}

  friend this_t& operator>>(std::istream& is, this_t& ar){
    istream_loop<this_t>::apply(is,ar); return ar;}

  void operator++(int){increment_loop<this_t>::apply(*this);}

  void operator--(int){decrement_loop<this_t>::apply(*this);}

  inline friend int size(const this_t& ar){return dim;}

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

  inline friend xpr< tt<int,this_t>  > operator*(const int&  l_, const this_t& r_){return xpr< tt<int,this_t>  >(l_,r_);}
  inline friend xpr< tt<Real,this_t> > operator*(const Real& l_, const this_t& r_){return xpr< tt<Real,this_t> >(l_,r_);}
  inline friend xpr< tt<Cplx,this_t> > operator*(const Cplx& l_, const this_t& r_){return xpr< tt<Cplx,this_t> >(l_,r_);}

  //===== Produit scalaire
  template <class r_t>typename resop<this_t,r_t>::type operator,(const r_t& r_) const {
    assert(nr==r_t::nr); return dprod<this_t,r_t,nr>::apply(*this,r_); }


};


//==========================//
//           Array          //
//==========================//

  template <int Dim, class V_t>
  class array{

  public:
    static const int dim =      Dim;
    static const int nr  =      Dim;
    typedef array<Dim,V_t>   this_t;
    typedef V_t                 v_t;

  private:
    v_t  v_[dim];

  public:
    array<Dim,V_t>(){ construct_loop<this_t>::apply(*this); }

    operator v_t* const(){return v_;}

    template <class r_t> array<Dim,V_t>(const r_t& r_){
      assign_loop<this_t,r_t>::apply(*this,r_);}

    template <class r_t> void operator=(const r_t& r_){
      assign_loop<this_t,r_t>::apply(*this,r_);}

    template <class r_t> bool operator==(const r_t& r_) const {
      for(int j=0; j<dim; j++){if((*this)[j]!=r_[j]){return false;}}
      return true;}

    template <class r_t> bool operator<(const r_t& r_) const {
      for(int j=0; j<dim; j++){
	if((*this)[j]>r_[j]){return false;}
	if((*this)[j]<r_[j]){return true;}}
      return false;}
    
    template <class r_t> void operator+=(const r_t& r_){
      plus_assign_loop<this_t,r_t>::apply(*this,r_);}

    v_t& operator[](const int& j){assert(j>=0 && j<nr);  return v_[j];}
    const v_t& operator[](const int& j) const {assert(j>=0 && j<nr); return v_[j];}

    inline friend std::ostream& operator<<(std::ostream& os, const this_t& ar){
      ostream_loop<this_t>::apply(os,ar); return os;}

    inline friend this_t& operator>>(std::istream& is, this_t& ar){
      istream_loop<this_t>::apply(is,ar); return ar;}

    void operator++(int){increment_loop<this_t>::apply(*this);}

    void operator--(int){decrement_loop<this_t>::apply(*this);}

    inline friend int size(const this_t& ar){return dim;}

    template <class i_t> subarray<this_t,i_t> operator[] (const i_t& i_){
      return subarray<this_t,i_t>(*this,i_);}

    template <class i_t> subarray<const this_t,i_t> operator[] (const i_t& i_) const {
      return subarray<const this_t,i_t>(*this,i_);}

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

    inline friend xpr< tt<int,this_t>  > operator*(const int&  l_, const this_t& r_){return xpr< tt<int,this_t>  >(l_,r_);}
    inline friend xpr< tt<Real,this_t> > operator*(const Real& l_, const this_t& r_){return xpr< tt<Real,this_t> >(l_,r_);}
    inline friend xpr< tt<Cplx,this_t> > operator*(const Cplx& l_, const this_t& r_){return xpr< tt<Cplx,this_t> >(l_,r_);}


    //===== Produit scalaire
    template <class r_t>typename resop<this_t,r_t>::type operator,(const r_t& r_) const {
      assert(nr==r_t::nr); return dprod<this_t,r_t,nr>::apply(*this,r_); }


  };




//%%%%%%%%%%%%%%%%%%%%%%%%%%//
//   Definitions de type    //
//%%%%%%%%%%%%%%%%%%%%%%%%%%//

typedef array<1,   int>     N1;
typedef array<2,   int>     N2;
typedef array<3,   int>     N3;
typedef array<4,   int>     N4;
typedef array<5,   int>     N5;
typedef array<6,   int>     N6;
typedef array<7,   int>     N7;
typedef array<8,   int>     N8;
typedef array<9,   int>     N9;
typedef array<10,  int>     N10;

typedef array<1,   Real>    R1;
typedef array<2,   Real>    R2;
typedef array<3,   Real>    R3;
typedef array<4,   Real>    R4;
typedef array<5,   Real>    R5;
typedef array<6,   Real>    R6;
typedef array<7,   Real>    R7;
typedef array<8,   Real>    R8;
typedef array<9,   Real>    R9;
typedef array<10,  Real>    R10;

typedef array<1,   Cplx>    C1;
typedef array<2,   Cplx>    C2;
typedef array<3,   Cplx>    C3;
typedef array<4,   Cplx>    C4;
typedef array<5,   Cplx>    C5;
typedef array<6,   Cplx>    C6;
typedef array<7,   Cplx>    C7;
typedef array<8,   Cplx>    C8;
typedef array<9,   Cplx>    C9;
typedef array<10,  Cplx>    C10;



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//   Construction a la volee    //
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//


inline N1 N1_(const int& i0){
  N1 I; I[0]=i0; return I;}
inline N2 N2_(const int& i0,const int& i1){
  N2 I; I[0]=i0;I[1]=i1; return I;}
inline N3 N3_(const int& i0, const int& i1, const int& i2){
  N3 I; I[0]=i0;I[1]=i1;I[2]=i2; return I;}
inline N4 N4_(const int& i0, const int& i1, const int& i2, const int& i3){
  N4 I; I[0]=i0;I[1]=i1;I[2]=i2;I[3]=i3; return I;}


inline R1 R1_(const Real& i0){
  R1 I; I[0]=i0; return I;}
inline R2 R2_(const Real& i0,const Real& i1){
  R2 I; I[0]=i0;I[1]=i1; return I;}
inline R3 R3_(const Real& i0, const Real& i1, const Real& i2){
  R3 I; I[0]=i0;I[1]=i1;I[2]=i2; return I;}
inline R4 R4_(const Real& i0, const Real& i1, const Real& i2, const Real& i3){
  R4 I; I[0]=i0;I[1]=i1;I[2]=i2;I[3]=i3; return I;}


inline C1 C1_(const Cplx& i0){
  C1 I; I[0]=i0; return I;}
inline C2 C2_(const Cplx& i0,const Cplx& i1){
  C2 I; I[0]=i0;I[1]=i1; return I;}
inline C3 C3_(const Cplx& i0, const Cplx& i1, const Cplx& i2){
  C3 I; I[0]=i0;I[1]=i1;I[2]=i2; return I;}
inline C4 C4_(const Cplx& i0, const Cplx& i1, const Cplx& i2, const Cplx& i3){
  C4 I; I[0]=i0;I[1]=i1;I[2]=i2;I[3]=i3; return I;}


} // namespace bemtool

#endif
