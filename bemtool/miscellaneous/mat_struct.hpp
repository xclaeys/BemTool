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
#ifndef BEMTOOL_MISC_MAT_STRUCT_HPP
#define BEMTOOL_MISC_MAT_STRUCT_HPP
#include "misc.hpp"
#include <map>
#include <vector>
#include <list>

namespace bemtool{

  //=====================//
  //  Matrices creuses   //
  //=====================//
  
  template <typename ValType>
  class SparseMatrix{

  public:
    typedef ValType                         v_t;
    typedef std::vector<v_t>                VectorType;
    typedef SparseMatrix<v_t>               this_t;
    typedef std::map<N2,ValType>            map_t;
    typedef typename map_t::iterator        it_t;
    typedef typename map_t::const_iterator  const_it_t;   
    
  private:    
    map_t data;
    int   nr;
    int   nc;

    std::vector<int> ir;
    std::vector<int> ic;
    
  public:
    SparseMatrix<ValType>(const int& nr0,const int& nc0): nr(nr0), nc(nc0){}
   
    v_t operator()(const int& j, const int& k) const {
      N2 I; I[0]=j; I[1]=k;
      const_it_t it = data.find(I);
      if(it==data.end()){return 0.;}
      return it->second; 
    }
    
    template <typename r_t, typename c_t>
    mat<r_t::dim,c_t::dim,v_t>
    operator()(const r_t& r, const c_t& c) const {
    mat<r_t::dim,c_t::dim,v_t> m;
      for(int j=0; j<r_t::dim; j++){
	for(int k=0; k<c_t::dim; k++){
	  m(j,k) = (*this)(r[j],c[k]);
	}
      }
      return m;
    }
    
    template <typename r_t, typename c_t>
    void Insert(const r_t& r, const c_t& c,
		const mat<r_t::dim,c_t::dim,v_t>& m){
      for(int j=0; j<r_t::dim; j++){
	for(int k=0; k<c_t::dim; k++){
	  this->Insert(r[j],c[k],m(j,k));
	}
      }
    }

    void Insert(const int& j, const int& k, const v_t& v){
      N2 I; I[0]=j; I[1]=k;
      it_t it = data.find(I);      
      if(it==data.end()){data[I]=v;}
      else{it->second += v;}
    }

    void Set(const int& j, const int& k, const v_t& v){
      N2 I; I[0]=j; I[1]=k; data[I]=v; }

    template <typename r_t, typename c_t>
    void Set(const r_t& r, const c_t& c,
	     const mat<r_t::dim,c_t::dim,v_t>& m){
      for(int j=0; j<r_t::dim; j++){
	for(int k=0; k<c_t::dim; k++){
	  this->Set(r[j],c[k],m(j,k));
	}
      }
    }
    
    friend int NbRows(const this_t& m){return m.nr;}
    friend int NbCols(const this_t& m){return m.nc;}    

    friend std::ostream& operator<<(std::ostream& o, const this_t& m){
      for(int j=0; j<NbRows(m); j++){
	for(int k=0; k<NbCols(m); k++){
	  o << m(j,k) << "\t";
	}
	o << "\n";
      }
      return o;
    }

    friend void mv_prod(VectorType& b,
			const this_t& m,
			const VectorType& x){
      assert(x.size()==NbCols(m) );
      b.clear(); b.resize(NbRows(m),v_t());
      const map_t& data = m.data;
      for(const_it_t it=data.begin(); it!=data.end(); it++){
	int  j    = it->first[0];
	int  k    = it->first[1];
	Cplx m_jk = it->second;
	b[j] += m_jk*x[k];
      }
    }
    
    friend void mv_prod(const std::vector<int>& ir,
			const std::vector<int>& ic,
			VectorType& b,
			const this_t& m,
			const VectorType& x){
      assert(ic.size()==NbCols(m) && ir.size()==NbRows(m));
      const map_t& data = m.data;
      for(const_it_t it=data.begin(); it!=data.end(); it++){
	int  j    = it->first[0];
	int  k    = it->first[1];
	Cplx m_jk = it->second;
	b[ir[j]] += m_jk*x[ic[k]];
      }
    }


    void Export(const char* filename) const {
      std::ofstream file;
      file.open(filename);
      for(const_it_t it = data.begin(); it!= data.end(); it++){
	int  j    = it->first[0];
	int  k    = it->first[1];
	Cplx m_jk = it->second;
	
	file << j << "\t" << k << "\t";
	file << m_jk.real() << "\t" << m_jk.imag() << std::endl;
      }
      file.close();
    }
    
    ValType Norm(const VectorType& u){
      ValType res = ValType();
      VectorType v; mv_prod(v,*this,u);
      for(int j=0; j<u.size();j++){res+=v[j]*u[j].conj();}
      return res;
    }
    
    
  };

  
  //=====================//
  //  Matrices denses    //
  //=====================//
  
  template <typename ValType>
  class DenseMatrix{

  public:
    typedef ValType             v_t;
    typedef std::vector<v_t>    VectorType;
    typedef DenseMatrix<v_t>    this_t;

  private:
    VectorType data;
    int        nr;
    int        nc;

    std::vector<int> ir;
    std::vector<int> ic;
        
  public:
    DenseMatrix<ValType>(const int& nr0,
			 const int& nc0): nr(nr0), nc(nc0), data(nr0*nc0,0.){}

    v_t& operator()(const int& j, const int& k){
      assert( 0<=j && j<nr && 0<=k && k<nc);      
      return data[j*nc+k];}

    const v_t& operator()(const int& j, const int& k) const {
      assert( 0<=j && j<nr && 0<=k && k<nc);      
      return data[j*nc+k];}

    template <class r_t, class c_t> submat<this_t,r_t,c_t>
    operator()(const r_t& I0,const c_t& J0 ){
      return submat<this_t,r_t,c_t>(*this,I0,J0);}

    template <class r_t, class c_t> submat<const this_t,r_t,c_t>
    operator()(const r_t& I0,const c_t& J0 ) const{
      return submat<const this_t,r_t,c_t>(*this,I0,J0);}
    
    friend int NbRows(const this_t& m){return m.nr;}
    friend int NbCols(const this_t& m){return m.nc;}    
    
    friend std::ostream& operator<<(std::ostream& o, const this_t& m){
      for(int j=0; j<NbRows(m); j++){
	for(int k=0; k<NbCols(m); k++){
	  o << m(j,k) << "\t";
	}
	o << "\n";
      }
      return o;
    }
    
    friend void mv_prod(VectorType& b,
			const this_t& m,
			const VectorType& x){
      int nr = NbRows(m);
      int nc = NbCols(m);
      
      assert(x.size()==nc);
      b.clear(); b.resize(nr,v_t());
      for(int j=0; j<nr; j++){
	for(int k=0; k<nc; k++){
	  b[j]+= m(j,k)*x[k]; 
	}
      }      
    }



    friend void mv_prod(const std::vector<int>& ir,
			const std::vector<int>& ic,
			VectorType& b,
			const this_t& m,
			const VectorType& x){
      int nr = NbRows(m);
      int nc = NbCols(m);
      assert(ic.size()==nc && ir.size()==nr);
      for(int j=0; j<nr; j++){
	for(int k=0; k<nc; k++){
	  b[ir[j]]+= m(j,k)*x[ic[k]]; 
	}
      }      
    }


    ValType Norm(const VectorType& u){
      ValType res = ValType();
      VectorType v; mv_prod(v,*this,u);
      for(int j=0; j<u.size();j++){res+=v[j]*u[j].conj();}
      return res;
    }
    friend int matlab_save(const DenseMatrix& m, const std::string& file){
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
    
  };
  

 
  //========================//
  //  Matrices heterogenes  //
  //========================//

  template<typename> class BlockMatrix;

  
  template <typename ValType>
  class BlockInserter{

    BlockMatrix<ValType>*    bm;
    const std::vector<int>&  ir;
    const std::vector<int>&  ic;

  public:
    BlockInserter<ValType>(const std::vector<int>&    ir0,
			   const std::vector<int>&    ic0,
			   BlockMatrix<ValType>* bm0): ir(ir0),ic(ic0),bm(bm0){};
    
    void operator+=(const DenseMatrix<ValType>& m){
      assert( ir.size()==NbRows(m) && ic.size()==NbCols(m) );
      bm->dmat.push_back(&m);
      bm->drow.push_back(ir);
      bm->dcol.push_back(ic);	
    }
    
    void operator+=(const SparseMatrix<ValType>& m){
      assert( ir.size()==NbRows(m) && ic.size()==NbCols(m) );
      bm->smat.push_back(&m);
      bm->srow.push_back(ir);
      bm->scol.push_back(ic);	
    }
    
  };


  
  
  template <typename ValType>
  class BlockMatrix{

    typedef ValType             v_t;
    typedef std::vector<v_t>    VectorType;
    typedef BlockMatrix<v_t>    this_t;
    
  private:
    friend class BlockInserter<ValType>;    
    std::vector<const DenseMatrix<Cplx>*  > dmat;
    std::vector< std::vector<int> >         drow;
    std::vector< std::vector<int> >         dcol;    
    std::vector<const SparseMatrix<Cplx>* > smat;
    std::vector< std::vector<int> >         srow;
    std::vector< std::vector<int> >         scol;
    
    int nr;
    int nc;
    
  public:
    BlockMatrix(const int& nr0, const int& nc0): nr(nr0), nc(nc0) {}
    
    friend int NbRows(const this_t& m){return m.nr;}
    friend int NbCols(const this_t& m){return m.nc;}    
    
    BlockInserter<ValType> operator()(const std::vector<int>& ir,
				      const std::vector<int>& ic){
      return BlockInserter<ValType>(ir,ic,this);}
    
    friend void mv_prod(VectorType& b, const this_t& m, const VectorType& x){    
      assert( x.size()== NbCols(m) );
      b.clear(); b.resize(NbRows(m));
      for(int j=0; j<m.dmat.size(); j++){
	mv_prod((m.drow)[j],(m.dcol)[j],b,*(m.dmat[j]),x);}
      for(int j=0; j<m.smat.size(); j++){
	mv_prod((m.srow)[j],(m.scol)[j],b,*(m.smat[j]),x);}
    }
    
  };

  
  

}

#endif
