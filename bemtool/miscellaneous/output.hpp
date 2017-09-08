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
#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include "../mesh/mesh.hpp"

namespace bemtool {

void WritePointValGmsh(const Mesh2D& mesh,
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
  for(int j=0; j<nb_elt; j++){
    file << j+1 << "\t2\t2\t2\t2\t" << Num(mesh[j])+1 << std::endl;}
  file << "$EndElements\n$NodeData\n";
  file << "1\n\"Normal field\"\n1\n0.0\n3\n0\n1\n";
  file << nb_elt << std::endl;
  for(int j=0; j<nb_elt; j++){file << j+1 << "\t" << x[j] << std::endl;}
  file << "$EndNodeData\n";
}

void WriteEltValGmsh(const Mesh2D& mesh,
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
  for(int j=0; j<nb_elt; j++){
    file << j+1 << "\t2\t2\t2\t2\t" << Num(mesh[j])+1 << std::endl;}
  file << "$EndElements\n$ElementData\n";
  file << "1\n\"Normal field\"\n1\n0.0\n3\n0\n1\n";
  file << nb_elt << std::endl;
  for(int j=0; j<nb_elt; j++){file << j+1 << "\t" << x[j] << std::endl;}
  file << "$EndElementData\n";
}

void WriteEltValGmsh(const Mesh2D& mesh,
		     char const * const name,
		     const std::vector<R3>& v){

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
  for(int j=0; j<nb_elt; j++){
    file << j+1 << "\t2\t2\t2\t2\t" << Num(mesh[j])+1 << std::endl;}
  file << "$EndElements\n$ElementData\n";
  file << "1\n\"Normal field\"\n1\n0.0\n3\n0\n3\n";
  file << nb_elt << std::endl;
  for(int j=0; j<nb_elt; j++){file << j+1 << "\t" << v[j] << std::endl;}
  file << "$EndElementData\n";
}


}






#endif
