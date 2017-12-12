#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include "bemtool/tools.hpp"
#include <fstream>

using namespace bemtool;


int main(int argc, char* argv[]){

  Real kappa = 0.5;
  Real kappa2 = kappa*kappa;
  Geometry node("mesh/circle.msh");
  Mesh1D mesh; mesh.Load(node,1);
  Orienting(mesh);
  int nb_elt = NbElt(mesh);
  std::cout << "nb_elt:\t" << nb_elt << std::endl;

  typedef LA_DL_2D_P1xP1 OperatorType;

  BIOp<OperatorType> V(mesh,mesh,kappa);
  Dof<P1_1D> dof(mesh);
  int nb_dof = NbDof(dof);
  std::cout << "nb_dof:\t" << nb_dof << std::endl;
  EigenDense  A(nb_dof,nb_dof); Clear(A);

  const int n = 1;
  std::vector<Cplx> En(nb_dof);
  for(int j=0; j<nb_elt; j++){
    const N2&         jdof = dof[j];
    const array<2,R3> xdof = dof(j);
    for(int k=0; k<2; k++){
      En[jdof[k]] = pow( xdof[k][0]+iu*xdof[k][1], n);}
  }

  /*
  progress bar("Assemblage",nb_elt);
  for(int j=0; j<nb_elt; j++){
    bar++;
    for(int k=0; k<nb_elt; k++){
      A(dof[j],dof[k]) += V(j,k);
    }
  }
  */
  progress bar("Assemblage",nb_dof);
  for(int j=0; j<nb_dof; j++){
    bar++;
    for(int k=0; k<nb_dof; k++){
      A(j,k) += V(dof.ToElt(j),dof.ToElt(k));
    }
  }  
  bar.end();
  
  Cplx sum =0.;
  for(int j=0; j<nb_dof; j++){
    for(int k=0; k<nb_dof; k++){
      sum+= A(j,k)*conj(En[j])*En[k];
    }
  }

  //  Cplx refsol = RefSol<OperatorType>::Compute(n,1.,kappa);
  Cplx refsol = RefEigenvalue<OperatorType>::Compute(n,1.,kappa);
  std::cout << "erreur relative:\t";
  //std::cout << abs(sum) << " %" << std::endl;
  std::cout << 100*sqrt( abs(sum - refsol)/abs(refsol) )<< " %" << std::endl;

}
