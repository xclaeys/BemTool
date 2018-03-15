//===================================================================
//
//  Copyright 2017 Xavier Claeys
//
//  bemtool is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  bemtool is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with bemtool.  If not, see <http://www.gnu.org/licenses/>
//
//====================================================================
#ifndef BEMTOOL_MISC_OUTPUT_GMSH_HPP
#define BEMTOOL_MISC_OUTPUT_GMSH_HPP

#include "../mesh/mesh.hpp"
namespace bemtool {

  
  void WriteEltVectGmsh(const Dof<RT0_2D>& dof,
			char const * const name,
			const std::vector<Real>& x)
  {
    
    const Mesh2D&   mesh   = MeshOf(dof);
    const Geometry& node = GeometryOf(mesh);
    RT0_2D phi(mesh);
    int nb_dof   = NbDof(dof);
    int nb_elt   = NbElt(dof);    
    int nb_node  = NbNode(node);
    int elt_type = 2;
        
    std::ofstream file; file.open(name);    
    file << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";
    file << nb_node << std::endl;
    for(int j=0; j<nb_node; j++){
      file << j+1 << "\t" << node[j] << std::endl;}
    file << "$EndNodes\n$Elements\n";
    file << nb_elt << std::endl;
    for(int j=0; j<nb_elt; j++){
      file << j+1 << "\t"<<elt_type<<"\t2\t2\t2\t" << Num(mesh[j])+1 << std::endl;}
    file << "$EndElements\n$ElementData\n";
    file << "1\n\"vector field\"\n1\n0.0\n3\n0\n3\n";
    file << nb_elt << std::endl;
    for(int j=0; j<nb_elt; j++){
      phi.Assign(j);
      const N3& I = dof[j];
      R2 t; t[0]=1./3.,t[1]=1./3.;
      R3 U;
      U += x[I[0]]*phi(0,t);
      U += x[I[1]]*phi(1,t);
      U += x[I[2]]*phi(2,t);      
      file << j+1 << "\t" << U << std::endl;
    }
    file << "$EndElementData\n";
        
  }
  

  
//// Gmsh
template<int dim, typename T>
void WritePointValGmsh(const Dof<BasisFct<P1,dim>>& dof,
		       char const * const name,
		       const std::vector<T>& x){

	int nb_dof = NbDof(dof);
  int nb_elt = NbElt(dof);
  std::ofstream file; file.open(name);
  file << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";
  file << nb_dof << std::endl;
  for(int j=0; j<nb_dof; j++){
    file << j+1 << "\t";
    file << dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][0]<<" ";
    file << dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][1]<<" ";
    file << dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][2]<< std::endl;}
  file << "$EndNodes\n$Elements\n";
  file << nb_elt << std::endl;
  int elt_type;
  switch (dim) {
      case 1:
      elt_type=1;
      break;
      case 2:
      elt_type=2;
      break;
      case 3:
      elt_type=4;
      break;
  }
  for(int j=0; j<nb_elt; j++){
    file << j+1 << "\t"<<elt_type<<"\t2\t2\t2\t" << dof[j]+1 << std::endl;}
  file << "$EndElements\n$NodeData\n";
  file << "1\n\"Normal field\"\n1\n0.0\n3\n0\n1\n";
  file << nb_dof << std::endl;
  for(int j=0; j<nb_dof; j++){file << j+1 << "\t" << x[j] << std::endl;}
  file << "$EndNodeData\n";
}

template<int dim, typename T>
void WriteEltValGmsh(const Dof<BasisFct<P1,dim>>& dof,
		       char const * const name,
		       const std::vector<T>& x){

	int nb_dof = NbDof(dof);
  int nb_elt = NbElt(dof);
  std::ofstream file; file.open(name);
  file << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";
  file << nb_dof << std::endl;
  for(int j=0; j<nb_dof; j++){
    file << j+1 << "\t" << dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][0]<<" "<<dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][1]<<" "<<dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][2]<< std::endl;}
  file << "$EndNodes\n$Elements\n";
  file << nb_elt << std::endl;
  int elt_type;
  switch (dim) {
      case 1:
      elt_type=1;
      break;
      case 2:
      elt_type=2;
      break;
      case 3:
      elt_type=4;
      break;
  }
  for(int j=0; j<nb_elt; j++){
    file << j+1 << "\t"<<elt_type<<"\t2\t2\t2\t" << dof[j]+1 << std::endl;}
  file << "$EndElements\n$NodeData\n";
    file << "1\n\"Normal field\"\n1\n0.0\n3\n0\n1\n";
    file << nb_elt << std::endl;
    for(int j=0; j<nb_elt; j++){file << j+1 << "\t" << x[j] << std::endl;}
    file << "$EndElementData\n";
}

template<typename T,int dim>
void WritePointValGmsh(const Mesh<dim>& mesh,
		       char const * const name,
		       const std::vector<T>& x){

  const Geometry& node = GeometryOf(mesh);
  int nb_node = NbNode(node);
  int nb_elt  = NbElt(mesh);
  std::ofstream file; file.open(name);
  file << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";
  file << NbNode(node) << std::endl;
  for(int j=0; j<nb_node; j++){
    file << j+1 << "\t" << node[j] << std::endl;}
  file << "$EndNodes\n$Elements\n";
  file << nb_elt << std::endl;
  int elt_type;
  switch (dim) {
      case 1:
      elt_type=1;
      break;
      case 2:
      elt_type=2;
      break;
      case 3:
      elt_type=4;
      break;
  }
  for(int j=0; j<nb_elt; j++){
    file << j+1 << "\t"<<elt_type<<"\t2\t2\t2\t" << Num(mesh[j])+1 << std::endl;}
  file << "$EndElements\n$NodeData\n";
  file << "1\n\"Normal field\"\n1\n0.0\n3\n0\n1\n";
  file << nb_node << std::endl;
  for(int j=0; j<nb_node; j++){file << j+1 << "\t" << x[j] << std::endl;}
  file << "$EndNodeData\n";
}

template<typename T,int dim>
void WriteEltValGmsh(const Mesh<dim>& mesh,
		     char const * const name,
		     const std::vector<Real>& x){

  const Geometry& node = GeometryOf(mesh);
  int nb_node = NbNode(node);
  int nb_elt  = NbElt(mesh);
  std::ofstream file; file.open(name);
  file << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";
  file << NbNode(node) << std::endl;
  for(int j=0; j<nb_node; j++){
    file << j+1 << "\t" << node[j] << std::endl;}
  file << "$EndNodes\n$Elements\n";
  file << nb_elt << std::endl;
  int elt_type;
  switch (dim) {
      case 1:
      elt_type=1;
      break;
      case 2:
      elt_type=2;
      break;
      case 3:
      elt_type=4;
      break;
  }
  for(int j=0; j<nb_elt; j++){
    file << j+1 << "\t"<<elt_type<<"\t2\t2\t2\t" << Num(mesh[j])+1 << std::endl;}
  file << "$EndElements\n$ElementData\n";
  file << "1\n\"Normal field\"\n1\n0.0\n3\n0\n1\n";
  file << nb_elt << std::endl;
  for(int j=0; j<nb_elt; j++){file << j+1 << "\t" << x[j] << std::endl;}
  file << "$EndElementData\n";
}




}






#endif
