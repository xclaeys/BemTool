#ifndef EIGEN_WRAP_H
#define EIGEN_WRAP_H

#include <Eigen/Dense>
#include "../calculus/calculus.hpp"
#include "misc.hpp"

using namespace std;

/*===========================
||   Interface avec eigen  ||
===========================*/

template <class EigenMatType,class ValType>
  class EigenMat{

private:
  EigenMatType     mat;
  vector<size_t>   ipvt;
  
public:
  
  typedef EigenMat<EigenMatType,ValType>   m_t;
  typedef ValType                          v_t;
  
  EigenMat(){};
  EigenMat(const int& n, const int& m): mat(n,m), ipvt(n){};
  
  v_t& operator()(const int& i, const int& j){
    return mat(i,j);}

  const v_t& operator()(const int& i, const int& j) const{
    return mat(i,j);}
  
  template <class r_t, class c_t> submat<m_t,r_t,c_t>
  operator()(const r_t& I0,const c_t& J0 ){
    return submat<m_t,r_t,c_t>(*this,I0,J0);}
  
  template <class r_t, class c_t> submat<const m_t,r_t,c_t>
  operator()(const r_t& I0,const c_t& J0 ) const{
    return submat<const m_t,r_t,c_t>(*this,I0,J0);}
  
  template <class c_t> submat<m_t,int,c_t>
  operator()(const int& I0,const c_t& J0 ){
    return submat<m_t,int,c_t>(*this,I0,J0);}

  template <class c_t> submat<const m_t,int,c_t>
  operator()(const int& I0,const c_t& J0 ) const {
    return submat<m_t,int,c_t>(*this,I0,J0);}
  
  template <class r_t> submat<m_t,r_t,int>
  operator()(const r_t& I0,const int& J0 ) {
    return submat<m_t,r_t,int>(*this,I0,J0);}

  template <class r_t> submat<const m_t,r_t,int>
  operator()(const r_t& I0,const int& J0 ) const {
    return submat<m_t,r_t,int>(*this,I0,J0);}
  
  void operator=(const m_t& m){
    for(int j=0; j<m.mat.rows(); j++){
      for(int k=0; k<m.mat.cols(); k++){
	mat(j,k) = m(j,k);
      }
    }
  }
  
  friend m_t operator+(const m_t& m1, const m_t& m2){
    m_t m; 
    for(int j=0; j<m.mat.rows(); j++){
      for(int k=0; k<m.mat.cols(); k++){
	m(j,k) = m1(j,k) + m2(j,k);
      }
    }
    return m;
  }

  friend m_t operator-(const m_t& m1, const m_t& m2){
    m_t m; 
    for(int j=0; j<m.mat.rows(); j++){
      for(int k=0; k<m.mat.cols(); k++){
	m(j,k) = m1(j,k) - m2(j,k);
      }
    }
    return m;
  }
  
  
  friend int NbRows(const m_t& m){
    return m.mat.rows();}
  
  friend int NbCols(const m_t& m){
    return m.mat.cols();}
  
  friend ostream& operator<<(ostream& os, m_t& m){
    for(int j=0; j<NbRows(m); j++){
      for(int k=0; k<NbCols(m); k++){
	os << m(j,k) << "\t\t";}
      os << endl;
    }
    return os;
  }
  
  friend void Clear(m_t& m){
    for(int j=0; j<NbRows(m); j++){
      for(int k=0; k<NbCols(m); k++){
	m(j,k)= 0.;	
      }
    }
  }
    
};


/*===========================
||   Definitions de type   ||
===========================*/

typedef Cplx Field;
static const int Dynamic = Eigen::Dynamic;
typedef Eigen::Matrix<Field, Dynamic, Dynamic>  DenseMatrix;
typedef EigenMat< DenseMatrix ,Field > EigenDense;

#endif
