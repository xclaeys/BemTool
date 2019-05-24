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
#ifndef BEMTOOL_MISC_EIGEN_WRAP_HPP
#define BEMTOOL_MISC_EIGEN_WRAP_HPP

#include <Eigen/Dense>
#include "../calculus/calculus.hpp"
#include "misc.hpp"

namespace bemtool {

  /*===========================
    ||   Interface avec eigen  ||
    ===========================*/

  template <class EigenMatType,class ValType>
  class EigenMat{

  private:
    EigenMatType          mat;
    std::vector<size_t>   ipvt;

  public:

    typedef EigenMat<EigenMatType,ValType>          m_t;
    typedef ValType                                 v_t;
    typedef std::vector<v_t>                        VectorType;
    typedef Eigen::Matrix<v_t,Eigen::Dynamic,1>     EigenVectorType;

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

    friend std::ostream& operator<<(std::ostream& os, m_t& m){
      for(int j=0; j<NbRows(m); j++){
	for(int k=0; k<NbCols(m); k++){
	  os << m(j,k) << "\t\t";}
	os << std::endl;
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

    friend void Write(m_t& m, char const * const name){
      std::ofstream file; file.open(name);
      int nb_row = NbRows(m);
      int nb_col = NbCols(m);

      file << "# name: " << name       << std::endl;
      file << "# type: complex matrix" << std::endl;
      file << "# rows:   \t" << nb_row << std::endl;
      file << "# columns:\t" << nb_col << std::endl;

      int count = 0;
      for(int k=0; k<nb_col; k++){
        for(int j=0; j<nb_row; j++){
      file << m(j,k) << "  ";
      count++;
      if(count==10){
      file << std::endl;
      count = 0;}
        }
      }
      file.close();
    }
    friend int matlab_save(m_t& m, const std::string& file){
        std::ofstream out(file);
        out << std::setprecision(18);
        if(!out) {
            std::cout << "Cannot open file."<<std::endl;
            return 1;
        }
        int rows = NbRows(m);
        int cols = NbCols(m);

        for (int i=0;i<rows;i++){
            for (int j=0;j<cols;j++){
                out<<std::real(m(i,j));
                if (std::imag(m(i,j))<0){
                    out<<std::imag(m(i,j))<<"i\t";
                }
                else if (std::imag(m(i,j))==0){
                    out<<"+"<<0<<"i\t";
                }
                else{
                    out<<"+"<<std::imag(m(i,j))<<"i\t";
                }
            }
            out << std::endl;
        }
        out.close();
        return 0;
    }
    /*=======================
      Singular value decomposition
      =======================*/
    
    std::vector<double> SVD(){
      Eigen::JacobiSVD<EigenMatType> svd(mat);
      const Eigen::VectorXd sv = svd.singularValues();
      std::vector<double> std_sv;
      for(int j=0; j<sv.size(); j++){ std_sv.push_back(sv[j]);}
      return std_sv;
    }

    /*=======================
      Produit matrice-vecteur
      =======================*/

    friend void mv_prod(VectorType& b0, const m_t& A, const VectorType& x0){
      assert( b0.size()==NbRows(A) && x0.size()==NbCols(A) );
      EigenVectorType x(NbCols(A));
      EigenVectorType b(NbRows(A));
      for(int j=0; j<NbCols(A); j++){x[j]=x0[j];} b = A.mat*x;
      for(int j=0; j<NbRows(A); j++){b0[j]=b[j];}
    }

    friend void add_mv_prod(VectorType& b0, const m_t& A, const VectorType& x0){
      assert( b0.size()==NbRows(A) && x0.size()==NbCols(A) );
      EigenVectorType x(NbCols(A));
      EigenVectorType b(NbRows(A));
      for(int j=0; j<NbCols(A); j++){x[j]=x0[j];} b = A.mat*x;
      for(int j=0; j<NbRows(A); j++){b0[j]+=b[j];}
    }


    /*========
      Solveurs
      ========*/

    friend void lu_solve(m_t& A, const VectorType& b0, VectorType& x0){
      assert( b0.size()==NbRows(A) && x0.size()==NbCols(A) );
      EigenVectorType b(NbRows(A));
      for(int j=0; j<NbRows(A); j++){b[j]=b0[j];}
      EigenVectorType x = A.mat.lu().solve(b);
      for(int j=0; j<NbCols(A); j++){x0[j]=x[j];}
    }

  };


  /*===========================
    ||   Definitions de type   ||
    ===========================*/

  typedef Cplx Field;
  static const int Dynamic = Eigen::Dynamic;
  typedef Eigen::Matrix<Field, Dynamic, Dynamic>  DenseMatrix;
  typedef EigenMat< DenseMatrix ,Field >          EigenDense;


}
#endif
