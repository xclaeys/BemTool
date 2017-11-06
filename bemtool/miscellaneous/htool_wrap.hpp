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

public:
  BIO_Generator(BIOp<KernelType>& V0, Dof<Discretization>& dof0):IMatrix(NbDof(dof0),NbDof(dof0)), V(V0), dof(dof0) {}

  Cplx get_coef(const int& i, const int& j) const {
    return V(dof.ToElt(i),dof.ToElt(j));
  }
};


template <typename KernelType, typename Discretization>
class POT_Generator : public htool::IMatrix<Cplx>{
  Potential<KernelType>& V;
  Dof<Discretization>& dof;
  Geometry& geometry;

public:
  POT_Generator(Potential<KernelType>& V0, Dof<Discretization>& dof0, Geometry& geometry0):IMatrix(NbDof(dof0),NbDof(dof0)), V(V0), dof(dof0), geometry(geometry0) {}

  Cplx get_coef(const int& i, const int& j) const {
    return V(geometry[i],dof.ToElt(j));
  }
};

}// namespace

#endif
