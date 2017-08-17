#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include "../mesh/mesh.hpp"

using namespace std;

void WritePointValGmsh(const Mesh2D& mesh,
		       char const * const name,
		       const vector<Real>& x){

  const Geometry& node = GeometryOf(mesh);
  int nb_node = NbNode(node);
  int nb_elt  = NbElt(mesh);
  ofstream file; file.open(name);        
  file << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";
  file << NbNode(node) << endl;
  for(int j=0; j<nb_node; j++){
    file << j+1 << "\t" << node[j] << endl;}
  file << "$EndNodes\n$Elements\n";
  file << nb_elt << endl;
  for(int j=0; j<nb_elt; j++){
    file << j+1 << "\t2\t2\t2\t2\t" << Num(mesh[j])+1 << endl;}
  file << "$EndElements\n$NodeData\n";
  file << "1\n\"Normal field\"\n1\n0.0\n3\n0\n1\n";
  file << nb_elt << endl;
  for(int j=0; j<nb_elt; j++){file << j+1 << "\t" << x[j] << endl;}
  file << "$EndNodeData\n"; 
}

void WriteEltValGmsh(const Mesh2D& mesh,
		     char const * const name,
		     const vector<Real>& x){

  const Geometry& node = GeometryOf(mesh);
  int nb_node = NbNode(node);
  int nb_elt  = NbElt(mesh);
  ofstream file; file.open(name);        
  file << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";
  file << NbNode(node) << endl;
  for(int j=0; j<nb_node; j++){
    file << j+1 << "\t" << node[j] << endl;}
  file << "$EndNodes\n$Elements\n";
  file << nb_elt << endl;
  for(int j=0; j<nb_elt; j++){
    file << j+1 << "\t2\t2\t2\t2\t" << Num(mesh[j])+1 << endl;}
  file << "$EndElements\n$ElementData\n";
  file << "1\n\"Normal field\"\n1\n0.0\n3\n0\n1\n";
  file << nb_elt << endl;
  for(int j=0; j<nb_elt; j++){file << j+1 << "\t" << x[j] << endl;}
  file << "$EndElementData\n"; 
}

void WriteEltValGmsh(const Mesh2D& mesh,
		     char const * const name,
		     const vector<R3>& v){
  
  const Geometry& node = GeometryOf(mesh);
  int nb_node = NbNode(node);
  int nb_elt  = NbElt(mesh);
  ofstream file; file.open(name);        
  file << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";
  file << NbNode(node) << endl;
  for(int j=0; j<nb_node; j++){
    file << j+1 << "\t" << node[j] << endl;}
  file << "$EndNodes\n$Elements\n";
  file << nb_elt << endl;
  for(int j=0; j<nb_elt; j++){
    file << j+1 << "\t2\t2\t2\t2\t" << Num(mesh[j])+1 << endl;}
  file << "$EndElements\n$ElementData\n";
  file << "1\n\"Normal field\"\n1\n0.0\n3\n0\n3\n";
  file << nb_elt << endl;
  for(int j=0; j<nb_elt; j++){file << j+1 << "\t" << v[j] << endl;}
  file << "$EndElementData\n"; 
}









#endif
