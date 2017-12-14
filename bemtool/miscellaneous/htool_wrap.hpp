#ifndef BEMTOOL_MISC_HTOOLWRAP_HPP
#define BEMTOOL_MISC_HTOOLWRAP_HPP

#include <htool/htool.hpp>
#include <bemtool/fem/dof.hpp>
#include <bemtool/operator/operator.hpp>

namespace bemtool{

template <typename KernelType, typename Discretization>
class BIO_Generator : public htool::IMatrix<Cplx>{
  BIOp<KernelType>& V;
  Dof<Discretization>& dof;
  SubBIOp<BIOp<KernelType>> subV;

public:
  BIO_Generator(BIOp<KernelType>& V0, Dof<Discretization>& dof0):IMatrix(NbDof(dof0),NbDof(dof0)), V(V0), dof(dof0), subV(V0,dof0,dof0) {}

  Cplx get_coef(const int& i, const int& j) const {
    return V(dof.ToElt(i),dof.ToElt(j));
  }

  htool::SubMatrix<Cplx> get_submatrix(const std::vector<int>& J, const std::vector<int>& K) {
    htool::SubMatrix<Cplx> mat(J,K);
    subV.compute_block(J,K,mat);
    // for (int i=0; i<mat.nb_rows(); i++)
    // for (int j=0; j<mat.nb_cols(); j++)
    //     mat(i,j) = this->get_coef(J[i], K[j]);
    return mat;
  }

};


template <typename KernelType, typename Discretization>
class POT_Generator : public htool::IMatrix<Cplx>{
  Potential<KernelType>& V;
  Dof<Discretization>& dof;
  Geometry& geometry;

public:
  POT_Generator(Potential<KernelType>& V0, Dof<Discretization>& dof0, Geometry& geometry0):IMatrix(NbNode(geometry0),NbDof(dof0)), V(V0), dof(dof0), geometry(geometry0) {}

  Cplx get_coef(const int& i, const int& j) const {
    return V(geometry[i],dof.ToElt(j));
  }

};

}// namespace

#endif
