//===================================================================
//
//  Copyright 2017 Xavier Claeys
//
//  bemtool is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  bemtool is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//  (see ../LICENSE.txt)
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with bemtool.  If not, see <http://www.gnu.org/licenses/>
//
//====================================================================
#ifndef BEMTOOL_MESH_NORMAL_H
#define BEMTOOL_MESH_NORMAL_H

#include <vector>
#include "mesh.hpp"
#include "adjacency.hpp"


namespace bemtool{



//==========================//
//       Normale            //
//==========================//

template<int D>
class Normal{

public:
  //______________
  // Alias de type
  static const int dim = D;
  typedef Normal<D>    this_t;
  typedef Geometry     g_t;
  typedef Mesh<D>      m_t;
  typedef Elt<D>       e_t;

private:
  //_________________________
  // Instances pre-existantes
  const g_t&                 node;
  const m_t&                 mesh;
  const std::vector<e_t>&         elt;

  //_______________
  // Donnee membres
  std::vector<bool>         orientation;
  std::vector<R3>           normal;
  static const R3           none;

private:

  R3& operator[](const int& j) { return normal[j];}

public:

  Normal(const m_t&);
  const R3& operator[](const int& j) const {return normal[j];}

  const R3& operator[](const e_t& e) const {
    const int& j = node[e][mesh];
    if(j<0){return none;}
    return normal[j];}

  void set(const int& j, const R3& n) {
    normal[j] = n;
  }

  inline friend const m_t& MeshOf(const this_t& N) {return N.mesh;}
  inline friend void swap(this_t& N, const int& j){N[j] = (-1.)*N[j]; }
  inline friend void swap(this_t& N){
    for(int j=0; j<N.normal.size(); j++){N[j] = (-1.)*N[j];} }

  friend void WriteGmsh(const this_t& N, char const * const name){
    const m_t& mesh = N.mesh;
    const g_t& node = N.node;
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
    for(int j=0; j<nb_elt; j++){file << j+1 << "\t" << N[j] << std::endl;}
    file << "$EndElementData\n";
  }

};

  template<int dim> const R3 Normal<dim>::none = 0.;

  //===========================//
  // Constructeur de la normale
  template<int dim>
  Normal<dim>::Normal(const m_t& m):
    mesh(m), node(GeometryOf(m)), elt(EltOf(m)) {

    int  nbelt  = NbElt(mesh);
    bool ok = true;
    orientation.assign(nbelt,ok);
    std::vector<bool> visited(nbelt,false);

    //====================================
    // Calcul de l'adjacence entre elements
    Adjacency<m_t> adj(mesh);

    //====================================
    // Recherche des composantes connexes
    Connected<m_t> component(mesh);
    int nbc = NbOf(component);
    Real global_orientation[nbc];

    //====================================
    // initialisation  de la recherche
    // d'un point extremal du maillage
    int  Iext = 0;
    Real Ext = norm2(Ctr( mesh[ component[Iext][0] ] ));




    //===============================//
    //   Breadth First Search sur    //
    //   chaque composante connexe   //
    //===============================//
    for(int I=0; I<nbc; I++){

      int nbe = component[I].size();
      int nb_visited = 0;
      std::queue<int> visit;

      int j0 = component[I][0];
      visit.push(j0);
      visited[j0]=true;

      while(nb_visited<nbe){

	j0 = visit.front();
	const e_t& e0 = mesh[j0];
	visit.pop();
	nb_visited++;

	if( norm2(Ctr(e0)) > Ext){
	  Ext = norm2(Ctr(e0));
	  Iext = I; }

	for(int k0=0; k0<dim+1; k0++){
	  const int& j1 = adj[j0][k0];
	  const e_t& e1 = mesh[j1];

	  if(!visited[j1]){
	    bool same = Comp(e1,e0);
	    if(same){orientation[j1]=orientation[j0];}
	    else{orientation[j1]=!orientation[j0];}
	    visited[j1]=true;
	    visit.push(j1);
	  }
	}
      }

      global_orientation[I] = 0.;
      const R3& p = mesh[ component[I][0] ][0];
      for(int j=0; j<nbe; j++){
	j0 = component[I][j];
	const e_t& e = mesh[j0];
	Real r = SolidAngle(p,e);
	if(!orientation[j0]){r = -r;}
	global_orientation[I] += r;
      }

      if(global_orientation[I]>0){
	for(int j=0; j<nbe; j++){
	  j0 = component[I][j];
	  orientation[j0] = !orientation[j0];
	}
      }

    }

    //=====================================
    // Calcul effectif des vecteurs normaux
    normal.resize(nbelt);
    for(int j=0; j<nbelt; j++){
      normal[j]=NormalTo(mesh[j]);
      normalize(normal[j]);
      if(orientation[j]){normal[j] = (-1.)*normal[j];}
    }

    //=====================================
    // Si le domaine est borne la
    // composante exterieure du bord
    // doit etre orientee dans l'autre sens
    if(mesh.IsBounded()){
      for(int j=0; j<component[Iext].size(); j++){
	int jj = component[Iext][j];
	normal[jj] = (-1.)*normal[jj];
      }
    }


  }



  template <int D> void Orienting(Mesh<D>& mesh){
    Normal<D> normal(mesh); mesh.Orienting(normal);}

  typedef Normal<1> Nrml1D;
  typedef Normal<2> Nrml2D;



} // namespace bemtool

#endif
