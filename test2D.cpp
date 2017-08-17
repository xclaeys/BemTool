#include <iostream>
#include <vector>
#include <map>
#include <fstream> 
#include <iostream>
#include <vector>
#include <map>
#include <tools.hpp>
#include <fstream> 

using namespace std;
using bemtool::array;


int main(int argc, char* argv[]){

  Real kappa = 0.5;
  Real kappa2 = kappa*kappa;
  Geometry node("mesh/circle.msh");
  Mesh1D mesh; mesh.Load(node,1);
  Orienting(mesh);
  int nb_elt = NbElt(mesh);
  cout << "nb_elt:\t" << nb_elt << endl;

  typedef HE_HS_2D_P1xP1 OperatorType;
  
  BIOp<OperatorType> V(mesh,mesh,kappa);
  Dof<P1_1D> dof(mesh);    
  int nb_dof = NbDof(dof);    
  cout << "nb_dof:\t" << nb_dof << endl;
  EigenDense  A(nb_dof,nb_dof); Clear(A);     
  
  const int n = 1;
  vector<Cplx> En(nb_dof);
  for(int j=0; j<nb_elt; j++){
    const N2&         jdof = dof[j];
    const array<2,R3> xdof = dof(j);
    for(int k=0; k<2; k++){
      En[jdof[k]] = pow( xdof[k][0]+iu*xdof[k][1], n);}   
  }
  
  progress bar("Assemblage",nb_elt);
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
  cout << "erreur relative:\t"; 
  //cout << abs(sum) << " %" << endl;
  cout << 100*sqrt( abs(sum - refsol)/abs(refsol) )<< " %" << endl;
  

}
