#ifndef BEMTOOL_MISC_HTOOLWRAP_HPP
#define BEMTOOL_MISC_HTOOLWRAP_HPP

#include <htool/htool.hpp>
#include <bemtool/fem/dof.hpp>
#include <bemtool/operator/operator.hpp>
#include <bemtool/operator/block_op.hpp>
#include <bemtool/potential/potential.hpp>

namespace bemtool{

template <typename KernelType, typename Discretization>
class BIO_Generator : public htool::VirtualGenerator<Cplx>{
  Dof<Discretization> dof;
  SubBIOp<BIOp<KernelType>> subV;
  // std::vector<int> boundary;

public:
    BIO_Generator(const Dof<Discretization>& dof0, const double& kappa):VirtualGenerator(NbDof(dof0),NbDof(dof0)), dof(dof0),subV(dof,dof,kappa) {}
    // {boundary=is_boundary_nodes(dof);}

  void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {
    SubBIOp<BIOp<KernelType>> subV_local = subV;
    subV_local.compute_block(M,N,rows,cols,ptr);
  }

};

template<int K>
class BIO_Generator<BIOpKernel<K,HS_OP,2,P0_1D,P0_1D>,P0_1D> : public htool::VirtualGenerator<Cplx>{
  public:
    BIO_Generator(const Dof<P0_1D>& dof0, const double& kappa):VirtualGenerator(0,0) {std::cout << "BemTool error: cannot use P0 discretization with Hyper Singular operator." << std::endl; assert(0);}
    void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {assert(0);}
};
template<int K>
class BIO_Generator<BIOpKernel<K,HS_OP,3,P0_2D,P0_2D>,P0_2D> : public htool::VirtualGenerator<Cplx>{
  public:
    BIO_Generator(const Dof<P0_2D>& dof0, const double& kappa):VirtualGenerator(0,0) {std::cout << "BemTool error: cannot use P0 discretization with Hyper Singular operator." << std::endl; assert(0);}
    void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {assert(0);}
};

template <typename KernelType, typename Discretization>
class BIO_Generator_w_mass : public htool::VirtualGenerator<Cplx>{
  Dof<Discretization> dof;
  SubBIOp<BIOp<KernelType>> subV;
  // std::vector<int> boundary;
  double coef;

public:
    BIO_Generator_w_mass(const Dof<Discretization>& dof0, const double& kappa,const double& coef0):VirtualGenerator(NbDof(dof0),NbDof(dof0)), dof(dof0),subV(dof,dof,kappa),coef(coef0) {}
    // {boundary=is_boundary_nodes(dof);}



  void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {
    SubBIOp<BIOp<KernelType>> subV_local = subV;
    subV_local.compute_block_w_mass(M,N,rows,cols,ptr,coef);

}

};

template<int K>
class BIO_Generator_w_mass<BIOpKernel<K,HS_OP,2,P0_1D,P0_1D>,P0_1D> : public htool::VirtualGenerator<Cplx>{
  public:
    BIO_Generator_w_mass(const Dof<P0_1D>& dof0, const double& kappa,const double& coef0):VirtualGenerator(0,0) {std::cout << "BemTool error: cannot use P0 discretization with Hyper Singular operator." << std::endl; assert(0);}
    void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {assert(0);}
};
template<int K>
class BIO_Generator_w_mass<BIOpKernel<K,HS_OP,3,P0_2D,P0_2D>,P0_2D> : public htool::VirtualGenerator<Cplx>{
  public:
    BIO_Generator_w_mass(const Dof<P0_2D>& dof0, const double& kappa,const double& coef0):VirtualGenerator(0,0) {std::cout << "BemTool error: cannot use P0 discretization with Hyper Singular operator." << std::endl; assert(0);}
    void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {assert(0);}
};

template <typename KernelType1, typename KernelType2, typename Discretization>
class Combined_BIO_Generator : public htool::VirtualGenerator<Cplx>{
  Dof<Discretization> dof;
  SubBIOp<BIOp<KernelType1>> sub1;
  SubBIOp<BIOp<KernelType2>> sub2;
  // std::vector<int> boundary;
  Cplx combined_coef_1,combined_coef_2;
  double mass_coef;

public:
    Combined_BIO_Generator(const Dof<Discretization>& dof0, const double& kappa,const Cplx& coef1,const Cplx& coef2,const double& mass_coef0):VirtualGenerator(NbDof(dof0),NbDof(dof0)), dof(dof0),sub1(dof,dof,kappa),sub2(dof,dof,kappa),combined_coef_1(coef1), combined_coef_2(coef2),mass_coef(mass_coef0) {}
    
    Combined_BIO_Generator(const Dof<Discretization>& dof0, const double& kappa,const Cplx& coef1,const double& mass_coef0):VirtualGenerator(NbDof(dof0),NbDof(dof0)), dof(dof0),sub1(dof,dof,kappa),sub2(dof,dof,kappa),combined_coef_1(coef1), combined_coef_2(1),mass_coef(mass_coef0) {}


  void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {
    SubBIOp<BIOp<KernelType1>> sub1_local = sub1;
    SubBIOp<BIOp<KernelType2>> sub2_local = sub2;
    sub1_local.compute_block(M,N,rows,cols,ptr);
    std::transform(ptr, ptr+M*N, ptr,
               std::bind(std::multiplies<Cplx>(), std::placeholders::_1, combined_coef_1));
    std::vector<Cplx> tmp(M*N,0);
    sub2_local.compute_block_w_mass(M,N,rows,cols,tmp.data(),mass_coef);
    int size = M*N;
    int incx(1), incy(1);
    htool::Blas<Cplx>::axpy(&size, &combined_coef_2, tmp.data(),&incx,ptr,&incy);



}

};

template<int K>
class Combined_BIO_Generator<BIOpKernel<K,HS_OP,2,P0_1D,P0_1D>,BIOpKernel<K,HS_OP,2,P0_1D,P0_1D>,P0_1D> : public htool::VirtualGenerator<Cplx>{
  public:
    Combined_BIO_Generator(const Dof<P0_1D>& dof0, const double& kappa,const Cplx& coef1,const Cplx& coef2,const double& mass_coef0):VirtualGenerator(0,0) {std::cout << "BemTool error: cannot use P0 discretization with Hyper Singular operator." << std::endl; assert(0);}
    Combined_BIO_Generator(const Dof<P0_1D>& dof0, const double& kappa,const Cplx& coef1,const double& mass_coef0):VirtualGenerator(0,0) {std::cout << "BemTool error: cannot use P0 discretization with Hyper Singular operator." << std::endl; assert(0);}
    void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {assert(0);}
};
template<int K, typename KernelType>
class Combined_BIO_Generator<BIOpKernel<K,HS_OP,2,P0_1D,P0_1D>,KernelType,P0_1D> : public htool::VirtualGenerator<Cplx>{
  public:
    Combined_BIO_Generator(const Dof<P0_1D>& dof0, const double& kappa,const Cplx& coef1,const Cplx& coef2,const double& mass_coef0):VirtualGenerator(0,0) {std::cout << "BemTool error: cannot use P0 discretization with Hyper Singular operator." << std::endl; assert(0);}
    Combined_BIO_Generator(const Dof<P0_1D>& dof0, const double& kappa,const Cplx& coef1,const double& mass_coef0):VirtualGenerator(0,0) {std::cout << "BemTool error: cannot use P0 discretization with Hyper Singular operator." << std::endl; assert(0);}
    void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {assert(0);}
};
template<int K, typename KernelType>
class Combined_BIO_Generator<KernelType,BIOpKernel<K,HS_OP,2,P0_1D,P0_1D>,P0_1D> : public htool::VirtualGenerator<Cplx>{
  public:
    Combined_BIO_Generator(const Dof<P0_1D>& dof0, const double& kappa,const Cplx& coef1,const Cplx& coef2,const double& mass_coef0):VirtualGenerator(0,0) {std::cout << "BemTool error: cannot use P0 discretization with Hyper Singular operator." << std::endl; assert(0);}
    Combined_BIO_Generator(const Dof<P0_1D>& dof0, const double& kappa,const Cplx& coef1,const double& mass_coef0):VirtualGenerator(0,0) {std::cout << "BemTool error: cannot use P0 discretization with Hyper Singular operator." << std::endl; assert(0);}
    void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {assert(0);}
};
template<int K>
class Combined_BIO_Generator<BIOpKernel<K,HS_OP,3,P0_2D,P0_2D>,BIOpKernel<K,HS_OP,3,P0_2D,P0_2D>,P0_2D> : public htool::VirtualGenerator<Cplx>{
  public:
    Combined_BIO_Generator(const Dof<P0_2D>& dof0, const double& kappa,const Cplx& coef1,const Cplx& coef2,const double& mass_coef0):VirtualGenerator(0,0) {std::cout << "BemTool error: cannot use P0 discretization with Hyper Singular operator." << std::endl; assert(0);}
    Combined_BIO_Generator(const Dof<P0_2D>& dof0, const double& kappa,const Cplx& coef1,const double& mass_coef0):VirtualGenerator(0,0) {std::cout << "BemTool error: cannot use P0 discretization with Hyper Singular operator." << std::endl; assert(0);}
    void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {assert(0);}
};
template<int K, typename KernelType>
class Combined_BIO_Generator<BIOpKernel<K,HS_OP,3,P0_2D,P0_2D>,KernelType,P0_2D> : public htool::VirtualGenerator<Cplx>{
  public:
    Combined_BIO_Generator(const Dof<P0_2D>& dof0, const double& kappa,const Cplx& coef1,const Cplx& coef2,const double& mass_coef0):VirtualGenerator(0,0) {std::cout << "BemTool error: cannot use P0 discretization with Hyper Singular operator." << std::endl; assert(0);}
    Combined_BIO_Generator(const Dof<P0_2D>& dof0, const double& kappa,const Cplx& coef1,const double& mass_coef0):VirtualGenerator(0,0) {std::cout << "BemTool error: cannot use P0 discretization with Hyper Singular operator." << std::endl; assert(0);}
    void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {assert(0);}
};
template<int K, typename KernelType>
class Combined_BIO_Generator<KernelType,BIOpKernel<K,HS_OP,3,P0_2D,P0_2D>,P0_2D> : public htool::VirtualGenerator<Cplx>{
  public:
    Combined_BIO_Generator(const Dof<P0_2D>& dof0, const double& kappa,const Cplx& coef1,const Cplx& coef2,const double& mass_coef0):VirtualGenerator(0,0) {std::cout << "BemTool error: cannot use P0 discretization with Hyper Singular operator." << std::endl; assert(0);}
    Combined_BIO_Generator(const Dof<P0_2D>& dof0, const double& kappa,const Cplx& coef1,const double& mass_coef0):VirtualGenerator(0,0) {std::cout << "BemTool error: cannot use P0 discretization with Hyper Singular operator." << std::endl; assert(0);}
    void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {assert(0);}
};


template <typename KernelType, typename Discretization>
class SubBIO_Neumann_Generator : public htool::VirtualGenerator<Cplx>{
    typedef typename BIOp<KernelType>::KernelTypeTrait  KernelTrait;
    typedef typename KernelTrait::ShapeFctX     PhiX;
    typedef typename KernelTrait::ShapeFctY     PhiY;
    typedef typename KernelTrait::MatType       MatType;
    typedef Dof<PhiX>                           DofX;
    typedef Dof<PhiY>                           DofY;

    static const int nb_dof_loc_x = KernelTrait::nb_dof_x;
    static const int nb_dof_loc_y = KernelTrait::nb_dof_y;
    typedef NDofLoc<nb_dof_loc_x> Nlocx;
    typedef NDofLoc<nb_dof_loc_y> Nlocy;
    typedef typename std::map<int,Nlocx>::iterator ItTypeX;
    typedef typename std::map<int,Nlocy>::iterator ItTypeY;


    Dof<Discretization> dof;
    SubBIOp<BIOp<KernelType>> subV;
    std::map<int,Nlocx>  Ix;
    std::map<int,Nlocy>  Iy;

    // std::vector<int> boundary;

public:
    SubBIO_Neumann_Generator(const Dof<Discretization>& dof0, const double& kappa, const std::vector<int>& targets0 , const std::vector<int>& sources0):VirtualGenerator(NbDof(dof0),NbDof(dof0)), dof(dof0),subV(dof,dof,kappa) {
        // boundary=is_boundary_nodes(dof);
        for(int k=0; k<targets0.size(); k++){
            const std::vector<N2>& jj = dof.ToElt(targets0[k]);
            for(int l=0; l<jj.size(); l++){
                const N2& j = jj[l]; Ix[j[0]][j[1]] = k;
            }
        }

        for(int k=0; k<sources0.size(); k++){
            const std::vector<N2>& jj = dof.ToElt(sources0[k]);
            for(int l=0; l<jj.size(); l++){
                const N2& j = jj[l]; Iy[j[0]][j[1]] = k;
            }
        }
    }



  void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {

    SubBIOp<BIOp<KernelType>> subV_local = subV;
    subV_local.compute_neumann_block(M,N,rows,cols,ptr,Ix,Iy);
  }

};

// template <typename KernelType, typename Discretization>
// class SubBIO_Generator : public htool::VirtualGenerator<Cplx>{
//   Dof<Discretization> dof;
//   SubBIOp<BIOp<KernelType>> subV;
//   std::vector<int> targets;
//   std::vector<int> sources;
//
// public:
//   SubBIO_Generator(const Dof<Discretization>& dof0, const double& kappa, const std::vector<int>& targets0 , const std::vector<int>& sources0):VirtualGenerator(targets0.size(),sources0.size()),dof(dof0),subV(dof,dof,kappa), targets(targets0),sources(sources0) {}
//
//   Cplx get_coef(const int& i, const int& j) const {
//     // Get the rank of the process
//     int rank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     std::vector<int> J(1,targets[i]);
//     std::vector<int> K(1,sources[j]);
//     if (rank==0){
//         std::cout << i<<" "<<j<<std::endl;
//         std::cout << targets[i]<<" "<<sources[j]<<std::endl;
//         // std::cout << i<<" "<<j<<std::endl;
//     }
//     htool::SubMatrix<Cplx> mat(J,K);
//     SubBIOp<BIOp<KernelType>> subV_local = subV;
//     subV_local.compute_block(J,K,mat);
//     return mat(0,0);
//   }
//
//   htool::SubMatrix<Cplx> get_submatrix(const std::vector<int>& J, const std::vector<int>& K) const{
//     std::vector<int> Jb(J.size(),0),Kb(K.size(),0);
//     for (int i=0; i<Jb.size();i++){
//         Jb[i]=targets[J[i]];
//     }
//     for (int i=0; i<Kb.size();i++){
//         Kb[i]=sources[K[i]];
//     }
//     htool::SubMatrix<Cplx> mat(Jb,Kb);
//     SubBIOp<BIOp<KernelType>> subV_local = subV;
//     subV_local.compute_block(Jb,Kb,mat);
//     return mat;
//   }
//
// };


template <typename KernelType, typename Discretization>
class POT_Generator : public htool::VirtualGenerator<Cplx>{
  Potential<KernelType> V;
  Dof<Discretization> dof;
  Geometry& geometry;

public:
  POT_Generator(Dof<Discretization>& dof0, Geometry& geometry0, double kappa):VirtualGenerator(NbNode(geometry0),NbDof(dof0)), dof(dof0), geometry(geometry0), V(MeshOf(dof0),kappa) {}

  void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {
    Potential<KernelType> V_local = V;
    for (int i=0; i<M; i++)
        for (int j=0; j<N; j++)
        {
            ptr[i+j*M] = V_local(geometry[rows[i]],dof.ToElt(cols[j]));
            }
  }

};

template<int K>
class POT_Generator<BIOpKernel<K,HS_OP,2,P0_1D,P0_1D>,P0_1D> : public htool::VirtualGenerator<Cplx>{
  public:
    POT_Generator(Potential<BIOpKernel<K,HS_OP,2,P0_1D,P0_1D>>& V0, Dof<P0_1D>& dof0, Geometry& geometry0):VirtualGenerator(0,0) {std::cout << "BemTool error: cannot use P0 discretization with Hyper Singular operator." << std::endl; assert(0);}
    void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {assert(0);}
};
template<int K>
class POT_Generator<BIOpKernel<K,HS_OP,3,P0_2D,P0_2D>,P0_2D> : public htool::VirtualGenerator<Cplx>{
  public:
    POT_Generator(Potential<BIOpKernel<K,HS_OP,3,P0_2D,P0_2D>>& V0, Dof<P0_2D>& dof0, Geometry& geometry0):VirtualGenerator(0,0) {std::cout << "BemTool error: cannot use P0 discretization with Hyper Singular operator." << std::endl; assert(0);}
    void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {assert(0);}
};

}// namespace

#endif
