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
#ifndef BEMTOOL_MISC_OUTPUT_PARAVIEW_HPP
#define BEMTOOL_MISC_OUTPUT_PARAVIEW_HPP

#include "../mesh/mesh.hpp"
namespace bemtool {

//// Paraview -- Ensight gold format
void WriteCaseParaview(char const * const name, char const * const meshname){
  std::ofstream file; file.open(name);
  file << "FORMAT" << std::endl;
  file << "type: ensight gold" << std::endl << std::endl;
  file << "GEOMETRY" << std::endl;
  file << "model: 1 "+std::string(meshname) << std::endl << std::endl;
}

void WriteCaseParaview(char const * const name, char const * const meshname, char const * const valname, char const * const valfile){
  std::ofstream file; file.open(name);
  file << "FORMAT" << std::endl;
  file << "type: ensight gold" << std::endl << std::endl;
  file << "GEOMETRY" << std::endl;
  file << "model: 1 "+std::string(meshname) << std::endl << std::endl;
  file << "VARIABLE"<<std::endl;
  file << "scalar per node: 1 "+std::string(valname)+" "+std::string(valfile)<<std::endl;
}

void WriteCaseParaview(char const * const name, char const * const meshname, char const * const valname, char const * const valfile, const std::vector<double>& times){
  std::ofstream file; file.open(name);
  file << "FORMAT" << std::endl;
  file << "type: ensight gold" << std::endl << std::endl;
  file << "GEOMETRY" << std::endl;
  file << "model: 1 "+std::string(meshname) << std::endl << std::endl;
  file << "VARIABLE"<<std::endl;
  file << "scalar per node: 1 "+std::string(valname)+" "+std::string(valfile)+"*"<<std::endl<<std::endl;
  file << "TIME"<<std::endl;
  file << "time set: 1"<<std::endl;
  file << "number of steps: "<<times.size()<<std::endl;
  file << "filename start number: 0"<<std::endl;
  file << "filename increment: 1"<<std::endl;
  file << "time values: ";
  for (int i=0;i<times.size();i++){
      file <<times[i]<<" ";
  }
  // std::copy(times.begin(),times.end(),std::ostream_iterator<double>(file, " "));
}

template<int dim>
void WriteMeshParaview(const Dof<BasisFct<P1,dim>>& dof,
		       char const * const name){

  int nb_dof = NbDof(dof);
  int nb_elt  = NbElt(dof);
  std::ofstream file; file.open(name);
  file << "description line 1" << std::endl;
  file << "description line 2" << std::endl;
  file << "node id given" << std::endl;
  file << "element id given" << std::endl;
  file << "part" << std::endl;
  file << 1 << std::endl;
  file << "description part 1" << std::endl;
  file << "coordinates" << std::endl;
  file << nb_dof << std::endl;
  for(int j=0; j<nb_dof; j++)
    file << j << std::endl;
  for(int j=0; j<nb_dof; j++)
    file << dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][0] << std::endl;
  for(int j=0; j<nb_dof; j++)
    file << dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][1] << std::endl;
  for(int j=0; j<nb_dof; j++)
    file << dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][2] << std::endl;
  switch (dim) {
      case 1:
      file << "bar2" << std::endl;
      break;
      case 2:
      file << "tria3" << std::endl;
      break;
      case 3:
      file << "tetra4" << std::endl;
      break;
  }

  file << nb_elt << std::endl;
  for(int j=0; j<nb_elt; j++)
    file << j << std::endl;
  for(int j=0; j<nb_elt; j++)
    file << dof[j]+1 << std::endl;
}

template<int dim>
void WritePointValParaview(const Dof<BasisFct<P1,dim>>& dof,
		       char const * const name,
               const std::vector<Real>& x){

  int nb_dof = NbDof(dof);
  int nb_elt  = NbElt(dof);
  std::ofstream file; file.open(name);
  file << "description line 1" << std::endl;
  file << "part" << std::endl;
  file << 1 << std::endl;
  file << "coordinates" << std::endl;
  for (int k=0;k<nb_dof;k++)
          file << x[k] << std::endl;
}

template<int dim>
void WriteMeshParaview(const Mesh<dim>& mesh,
		       char const * const name){

  const Geometry& node = GeometryOf(mesh);
  int nb_node = NbNode(node);
  int nb_elt  = NbElt(mesh);
  std::ofstream file; file.open(name);
  file << "description line 1" << std::endl;
  file << "description line 2" << std::endl;
  file << "node id given" << std::endl;
  file << "element id given" << std::endl;
  file << "part" << std::endl;
  file << 1 << std::endl;
  file << "description part 1" << std::endl;
  file << "coordinates" << std::endl;
  file << NbNode(node) << std::endl;
  for(int j=0; j<nb_node; j++)
    file << j << std::endl;
  for(int j=0; j<nb_node; j++)
    file << node[j][0] << std::endl;
  for(int j=0; j<nb_node; j++)
    file << node[j][1] << std::endl;
  for(int j=0; j<nb_node; j++)
    file << node[j][2] << std::endl;
  switch (dim) {
      case 1:
      file << "bar2" << std::endl;
      break;
      case 2:
      file << "tria3" << std::endl;
      break;
      case 3:
      file << "tetra4" << std::endl;
      break;
  }

  file << nb_elt << std::endl;
  for(int j=0; j<nb_elt; j++)
    file << j << std::endl;
  for(int j=0; j<nb_elt; j++)
    file << Num(mesh[j])+1 << std::endl;
}

template<int dim>
void WritePointValParaview(const Mesh<dim>& mesh,
		       char const * const name,
               const std::vector<Real>& x){

  const Geometry& node = GeometryOf(mesh);
  int nb_node = NbNode(node);
  int nb_elt  = NbElt(mesh);
  std::ofstream file; file.open(name);
  file << "description line 1" << std::endl;
  file << "part" << std::endl;
  file << 1 << std::endl;
  file << "coordinates" << std::endl;
  for (int k=0;k<nb_node;k++)
          file << x[k] << std::endl;
}




}






#endif
