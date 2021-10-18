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
#ifndef BEMTOOL_CALCULUS_SMALLMAT_HPP
#define BEMTOOL_CALCULUS_SMALLMAT_HPP
#include "smallvect.hpp"


namespace bemtool {

//==========================//
//         Sub-Matrix       //
//==========================//

template <class m_t, class ro_t, class co_t>
class submat{

public:

  static const int nr  = ro_t::nr;
  static const int nc  = co_t::nr;
  static const int dim = nr*nc;
  typedef typename m_t::v_t        v_t;
  typedef array<nr,v_t>   a_t;
  typedef submat<m_t,ro_t,co_t>    this_t;

private:

  m_t& m;
  const ro_t& I;
  const co_t& J;

public:

  submat<m_t,ro_t,co_t>(m_t& m0, const ro_t& I0, const co_t& J0): m(m0), I(I0), J(J0){};
  template <class r_t> void operator =(const r_t& r_){ assign_dbloop<this_t,r_t>::apply(*this,r_);      }
  template <class r_t> void operator+=(const r_t& r_){ plus_assign_dbloop<this_t,r_t>::apply(*this,r_); }
  v_t& operator()(const int& j, const int& k){assert(j>=0 && j<nr && k>=0 && k<nc); return m(I[j],J[k]); }
  const v_t& operator()(const int& j, const int& k) const {assert(j>=0 && j<nr && k>=0 && k<nc); return m(I[j],J[k]); }
  inline friend std::ostream& operator<<(std::ostream& os, const this_t& r_){
    ostream_dbloop<this_t>::apply(os,r_); return os;}

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

};




//===================================//
//   Specialisation Sub-Matrix       //
//   Cas d'une seule ligne           //
//===================================//

template <class m_t, class co_t>
class submat<m_t,int,co_t>{

public:

  static const int nr  = 1;
  static const int nc  = co_t::nr;
  static const int dim = nc;
  typedef typename m_t::v_t         v_t;
  typedef array<nr,v_t>    a_t;
  typedef submat<m_t,int,co_t>  this_t;

private:

  m_t&        m;
  const int   I;
  const co_t& J;

public:

  submat<m_t,int,co_t>(m_t& m0, const int& I0, const co_t& J0): m(m0), I(I0), J(J0){};
  template <class r_t> void operator =(const r_t& r_){ assign_dbloop<this_t,r_t>::apply(*this,r_);      }
  template <class r_t> void operator+=(const r_t& r_){ plus_assign_dbloop<this_t,r_t>::apply(*this,r_); }
  v_t& operator()(const int& j, const int& k)       {assert(j==0); return m(I,J[k]); }
  const v_t& operator()(const int& j, const int& k) const {assert(j==0); return m(I,J[k]); }
  inline friend std::ostream& operator<<(std::ostream& os, const this_t& r_){
    ostream_dbloop<this_t>::apply(os,r_); return os;}

};


//===================================//
//   Specialisation Sub-Matrix       //
//   Cas d'une seule colonne         //
//===================================//

template <class m_t, class ro_t>
class submat<m_t,ro_t,int>{

public:

  static const int nr  = ro_t::nr;
  static const int nc  = 1;
  static const int dim = nr;
  typedef typename m_t::v_t         v_t;
  typedef array<nr,v_t>    a_t;
  typedef submat<m_t,ro_t,int>      this_t;

private:

  m_t&        m;
  const ro_t& I;
  const int   J;

public:

  submat<m_t,ro_t,int>(m_t& m0, const ro_t& I0, const int& J0): m(m0), I(I0), J(J0){};
  v_t& operator[](const int& j) {assert(j>=0 && j<nr); return m(I[j],J);}
  const v_t& operator[](const int& j) const {assert(j>=0 && j<nr); return m(I[j],J);}
  template <class r_t> void operator= (const r_t& r_){assign_loop<this_t,r_t>::apply(*this,r_);}
  template <class r_t> void operator+=(const r_t& r_){plus_assign_loop<this_t,r_t>::apply(*this,r_);}
  friend std::ostream& operator<<(std::ostream& os, const this_t& ar){ostream_loop<this_t>::apply(os,ar); return os;}

};


//==========================//
//           Matrix         //
//==========================//

template <int Nr, int Nc, class V_t>
class mat{

public:
  static const int dim =  Nr*Nc;
  static const int nr  =     Nr;
  static const int nc  =     Nc;

  typedef V_t                v_t;
  typedef array<nr,v_t>      a_t;
  typedef mat<nr,nc,v_t>  this_t;

private:

  v_t  v_[dim];

public:

  mat<Nr,Nc,V_t>(){ construct_dbloop<this_t>::apply(*this); }
  template <class r_t> mat<Nr,Nc,V_t> (const r_t& r_){ assign_dbloop<this_t,r_t>::apply(*this,r_);  }
  template <class r_t> void operator= (const r_t& r_){ assign_dbloop<this_t,r_t>::apply(*this,r_);  }
  template <class r_t> void operator+=(const r_t& r_){plus_assign_dbloop<this_t,r_t>::apply(*this,r_); }
  v_t& operator()(const int& j, const int& k){assert(j>=0 && j<nr && k>=0 && k<nc); return v_[j+k*nr]; }
  const v_t& operator()(const int& j, const int& k) const {
    assert(j>=0 && j<nr && k>=0 && k<nc); return v_[j+k*nr]; }

  //==== Sous-matrices

  template <class ro_t, class co_t> submat<this_t,ro_t,co_t>
  operator()(const ro_t& I, const co_t& J){ return submat<this_t,ro_t,co_t>(*this,I,J);}

  template <class ro_t, class co_t> submat< const this_t,ro_t,co_t>
  operator()(const ro_t& I, const co_t& J) const { return submat<this_t,ro_t,co_t>(*this,I,J);}

  //==== Affichage

  inline friend std::ostream& operator<<(std::ostream& os, const this_t& m){
    ostream_dbloop<this_t>::apply(os,m); return os;}

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

};


//%%%%%%%%%%%%%%%%%%%%%%%%%%//
//    Produit scalaire      //
//%%%%%%%%%%%%%%%%%%%%%%%%%%//


template <class a_t> inline typename a_t::v_t norm2(const a_t& a_){return sqrt( (a_,a_) ); }
template <class a_t> inline void normalize(a_t& a_){a_ = (1./norm2(a_))*a_;}


//%%%%%%%%%%%%%%%%%%%%%%%%%%//
//  Inversion de matrices   //
//%%%%%%%%%%%%%%%%%%%%%%%%%%//

template<class m_t>
typename m_t::v_t det(const m_t& M){

    if(m_t::nr== 2){
        return M(0,0)*M(1,1)-M(1,0)*M(0,1);}

    if(m_t::nr== 3){
        return M(0,0)*( M(1,1)*M(2,2) - M(2,1)*M(1,2) )
        - M(0,1)*( M(1,0)*M(2,2)-M(2,0)*M(1,2) )
        + M(0,2)*( M(1,0)*M(2,1)-M(2,0)*M(1,1) ); }

    if(m_t::nr > 3){
        std::cout << "matrice trop grosse" << std::endl;
        std::exit(EXIT_FAILURE);}
}


template<class m_t>
mat<m_t::nr,m_t::nc,typename m_t::v_t>
inv(const m_t& M){

    mat<m_t::nr,m_t::nc,typename m_t::v_t> R;
    if(m_t::nr== 2){
        typename m_t::v_t Delta = M(0,0)*M(1,1)-M(1,0)*M(0,1);
        R(0,0) =  M(1,1)/Delta;
        R(1,1) =  M(0,0)/Delta;
        R(0,1) = -M(0,1)/Delta;
        R(1,0) = -M(1,0)/Delta;
        return R; }

    if(m_t::nr== 3){
        typename m_t::v_t Delta = M(0,0)*( M(1,1)*M(2,2) - M(2,1)*M(1,2) )
        + M(0,1)*( M(1,0)*M(2,2)-M(2,0)*M(1,2) )
        + M(0,2)*( M(1,0)*M(2,1)-M(2,0)*M(1,1) );
        R(0,0) =  ( M(1,1)*M(2,2) - M(2,1)*M(1,2) )/Delta;
        R(1,0) = -( M(1,0)*M(2,2) - M(2,0)*M(1,2) )/Delta;
        R(2,0) =  ( M(1,0)*M(2,1) - M(2,0)*M(1,1) )/Delta;
        R(0,1) = -( M(0,1)*M(2,2) - M(2,1)*M(0,2) )/Delta;
        R(1,1) =  ( M(0,0)*M(2,2) - M(2,0)*M(0,2) )/Delta;
        R(2,1) = -( M(0,0)*M(2,1) - M(2,0)*M(0,1) )/Delta;
        R(0,2) =  ( M(0,1)*M(1,2) - M(1,1)*M(0,2) )/Delta;
        R(1,2) = -( M(0,0)*M(1,2) - M(1,0)*M(0,2) )/Delta;
        R(2,2) =  ( M(0,0)*M(1,1) - M(1,0)*M(0,1) )/Delta;
        return R; }

    if(m_t::nr > 3){
        std::cout << "matrice trop grosse" << std::endl;
        std::exit(EXIT_FAILURE);}

}


template<int nr, int nc, class v_t>
mat<nc,nr,v_t> tr(const mat<nr,nc,v_t>& M){
    mat<nc,nr,v_t> R;
    for(int j=0; j<nr; j++){
        for(int k=0; k<nc; k++){
            R(k,j) = M(j,k);}}
    return R;}


//%%%%%%%%%%%%%%%%%%%%%//
//    Combinatoire     //
//%%%%%%%%%%%%%%%%%%%%%//

template <class array_t>
inline void swap(array_t& ar, const int& j, const int& k){
    typename array_t::v_t val = ar[j];
    ar[j] = ar[k]; ar[k] = val;}

template <class array_t>
inline void sort(array_t& ar){
    const int N = size(ar);
    for(int j=1;j<N;j++){
        for(int k=0;k<N-j;k++){
            if(ar[k]>ar[k+1]){swap(ar,k,k+1);}
        }
    }
}

template <int dim>
inline void init(array<dim,int>& I){
    ascending_loop< array<dim,int> >::apply(I);  }

//%%%%%%%%%%%%%%%%%%%%%%%%%%//
//    Produit vectoriel     //
//%%%%%%%%%%%%%%%%%%%%%%%%%%//

template <class l_t, class r_t>
inline array<3,typename resop<l_t,r_t>::type> vprod(const l_t& u, const r_t& v){
    array<3,typename resop<l_t,r_t>::type> w;
    w[0] = u[1]*v[2] - u[2]*v[1];
    w[1] = u[2]*v[0] - u[0]*v[2];
    w[2] = u[0]*v[1] - u[1]*v[0];
    return w; }


//%%%%%%%%%%%%%%%%%%%%%%%%%%//
//   Definitions de type    //
//%%%%%%%%%%%%%%%%%%%%%%%%%%//

typedef mat<1,1,   Real>    R1x1;
typedef mat<1,2,   Real>    R1x2;
typedef mat<2,1,   Real>    R2x1;
typedef mat<2,2,   Real>    R2x2;
typedef mat<3,3,   Real>    R3x3;
typedef mat<4,4,   Real>    R4x4;
typedef mat<5,5,   Real>    R5x5;
typedef mat<6,6,   Real>    R6x6;
typedef mat<10,10, Real>    R10x10;
typedef mat<3,1,   Real>    R3x1;
typedef mat<3,2,   Real>    R3x2;
typedef mat<1,3,   Real>    R1x3;
typedef mat<2,3,   Real>    R2x3;
typedef mat<5,2,   Real>    R5x2;
typedef mat<5,3,   Real>    R5x3;

typedef mat<1,1,   Cplx>    C1x1;
typedef mat<2,1,   Cplx>    C2x1;
typedef mat<1,2,   Cplx>    C1x2;
typedef mat<1,3,   Cplx>    C1x3;
typedef mat<1,4,   Cplx>    C1x4;
typedef mat<1,5,   Cplx>    C1x5;
typedef mat<2,2,   Cplx>    C2x2;
typedef mat<3,3,   Cplx>    C3x3;
typedef mat<4,4,   Cplx>    C4x4;
typedef mat<6,6,   Cplx>    C6x6;
typedef mat<10,10, Cplx>    C10x10;
typedef mat<3,2,   Cplx>    C3x2;
typedef mat<2,3,   Cplx>    C2x3;
typedef mat<3,1,   Cplx>    C3x1;
  


  template <int nr, int nc, typename v_t>
  mat<nc,nr,v_t> transpose(const mat<nr,nc,v_t>& m0){
    mat<nc,nr,v_t> m;
    for(int j=0; j<nr; j++){
      for(int k=0; k<nc; k++){
	m(k,j)=m0(j,k);
      }
    }
    return m;
  }


template <class x_t>
mat<x_t::nr,1, typename x_t::v_t> mat_(const x_t& x0){
    mat<x_t::nr,1, typename x_t::v_t> M;
    for(int j=0; j<x_t::nr; j++){M(j,0) = x0[j];}
    return M;
}

template <class x_t>
mat<x_t::nr,2, typename x_t::v_t> mat_(const x_t& x0, const x_t& x1){
    mat<x_t::nr,2, typename x_t::v_t> M;
    for(int j=0; j<x_t::nr; j++){
        M(j,0) = x0[j];
        M(j,1) = x1[j];
    }
    return M;
}

template <class x_t>
mat<x_t::nr,3, typename x_t::v_t> mat_(const x_t& x0, const x_t& x1, const x_t& x2){
    mat<x_t::nr,3, typename x_t::v_t> M;
    for(int j=0; j<x_t::nr; j++){
        M(j,0) = x0[j];
        M(j,1) = x1[j];
        M(j,2) = x2[j];
    }
    return M;
}

} // namespace bemtool

#endif
