#include <iostream>
#include <vector>
#include <map>
#include <tools.hpp>
#include <fstream>

using namespace bemtool;


template <typename OperatorType>
struct Test{
  static inline void Launch(){

    Real kappa = 1.;
    Real kappa2 = kappa*kappa;
    Geometry node("mesh/sphere.msh");
    Mesh2D mesh; mesh.Load(node,1);
    Orienting(mesh);
    int nb_elt = NbElt(mesh);
    std::cout << "nb_elt:\t" << nb_elt << std::endl;

    BIOp<OperatorType> V(mesh,mesh,kappa);
    Dof<P1_2D> dof(mesh);
    int nb_dof = NbDof(dof);
    std::cout << "nb_dof:\t" << nb_dof << std::endl;
    EigenDense  A(nb_dof,nb_dof); Clear(A);

    const int n = 1;
    const int m = 1;
    const N2 nm = N2_(n,m);
    std::vector<Cplx> En(nb_dof);
    for(int j=0; j<nb_elt; j++){
      const N3&         jdof = dof[j];
      const array<3,R3> xdof = dof(j);
      for(int k=0; k<3; k++){
	En[jdof[k]] = SphHarmo(nm,xdof[k]);}
    }

    progress bar("Assemblage\t",nb_elt);
    for(int j=0; j<nb_elt; j++){
      bar++;
      for(int k=0; k<nb_elt; k++){
	A(dof[j],dof[k]) += V(j,k);
      }
    }
    bar.end();

    Cplx sum =0.;
    for(int j=0; j<nb_dof; j++){
      for(int k=0; k<nb_dof; k++){
	sum+= A(j,k)*conj(En[j])*En[k];
      }
    }

    Cplx refsol = RefSol<OperatorType>::Compute(n,1.,kappa);
    std::cout << "Erreur relative:\t";
    std::cout << 100*abs(sum - refsol)/abs(refsol) << " %" << std::endl;
  }

};



int main(int argc, char* argv[]){

  Test<HE_SL_3D_P1xP1>::Launch();

}
