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
#ifndef BEMTOOL_FEM_DOF_HPP
#define BEMTOOL_FEM_DOF_HPP
#include "shapefct.hpp"



namespace bemtool{

template <int> struct EdgePerElt;
template <>    struct EdgePerElt<1>{static const int nb=1;};
template <>    struct EdgePerElt<2>{static const int nb=3;};
template <>    struct EdgePerElt<3>{static const int nb=6;};

template <typename> class Dof;
template <typename ShapeFct> struct DofTraits{

public:
  static const int dim             = ShapeFct::dim;
  static const int nb_dof_loc      = ShapeFct::nb_dof_loc;
  static const int nb_edge_loc     = EdgePerElt<dim>::nb;
  typedef array<nb_dof_loc,int>    Nloc;
  typedef array<nb_dof_loc,R3>     R3loc;
  typedef Mesh<dim>                mesh_t;

};


/*====================================================
                  NOTE EXPLICATIVE
  ====================================================
  Notons N1 = nb_elt
         N2 = nb_dof_loc et
	 N3 = nb_dof

  - le tableau elt_to_dof modelise une application
  surjective de [0,N1]x[0,N2] -> [0,N3]

  - le tableau dof_to_elt modelise l'application inverse
  a gauche [0,N3] -> [0,N1]x[0,N2]

  ===================================================*/


/*=========================
      Numerotation P0
  ========================*/

template <int Dim> class Dof< BasisFct<P0,Dim> >{

public:
  typedef BasisFct<P0,Dim>     ShapeFct;
  typedef Dof<ShapeFct>        this_t;
  typedef DofTraits<ShapeFct>  Trait;

private:
  const Geometry*                       node_p;
  const typename Trait::mesh_t*         mesh_p;
  std::vector<typename Trait::Nloc>     elt_to_dof;
  std::vector< std::vector<N2> >        dof_to_elt;
  int                                   nb_elt;
  int                                   nb_dof;
  int                                   offset;

public:
  Dof<ShapeFct>(): mesh_p(0), node_p(0), offset(0) {};

  Dof<ShapeFct>(const typename Trait::mesh_t& m):
  mesh_p(&m), node_p(&GeometryOf(m)), nb_elt(NbElt(m)), offset(0) {
    nb_dof = nb_elt;
    elt_to_dof.resize(nb_elt);
    dof_to_elt.resize(nb_dof);
    for(int j=0; j<nb_elt; j++){
      elt_to_dof[j]=j;
      dof_to_elt[j].push_back(N2_(j,0));
    }
  }

  typename Trait::R3loc operator()(const int& j) const {
    Elt<Trait::dim> e = (*mesh_p)[j];
    typename Trait::R3loc P; P[0] = Ctr(e);
    return P;
  }

  const typename Trait::Nloc&
  operator[](const int& j) const {
    return elt_to_dof[j];}

  const typename Trait::Nloc&
  operator[](const Elt<Trait::dim>& e) const {
    return (*this)[ (*node_p)[e][(*mesh_p)] ];}

  friend const int& NbDof(const this_t& d){return d.nb_dof;}
  friend const int& NbElt(const this_t& d){return d.nb_elt;}

  void operator+=(const int& offset0){
    offset += offset0;
    for(int j=0; j<nb_elt; j++){elt_to_dof[j]+=offset0;}}

  const std::vector<N2>& ToElt(const int& j) const {
    assert(j>=offset && j<offset+nb_dof);
    return dof_to_elt[j-offset];}

  friend const typename Trait::mesh_t&
  MeshOf(const this_t& dof){return *(dof.mesh_p);}

};




/*===============
  NUMEROTATION P1
  ===============*/

template <int Dim> class Dof< BasisFct<P1,Dim> >{

public:
  typedef BasisFct<P1,Dim>     ShapeFct;
  typedef Dof<ShapeFct>        this_t;
  typedef DofTraits<ShapeFct>  Trait;

private:
  const Geometry*                      node_p;
  const typename Trait::mesh_t*        mesh_p;
  std::vector<typename Trait::Nloc>    elt_to_dof;
  std::vector< std::vector<N2> >       dof_to_elt;
  int                                  nb_elt;
  int                                  nb_dof;
  int                                  offset;

public:
  Dof<ShapeFct>(): mesh_p(0), node_p(0), offset(0) {};

  Dof<ShapeFct>(const typename Trait::mesh_t& m):
  mesh_p(&m), node_p(&GeometryOf(m)), nb_elt(NbElt(m)), offset(0) {
    const int nb_node = NbNode(*node_p);
    const typename Trait::mesh_t&  mesh = *mesh_p;

    nb_dof = 0;
    elt_to_dof.resize(NbElt(mesh));
    std::vector<int> num(nb_node,-1);
    for(int j=0; j<NbElt(mesh); j++){
      array<Trait::dim+1,int> I = Num(mesh[j]);
      for(int k=0; k<Trait::dim+1; k++){
	if(num[I[k]]==-1){
	  num[I[k]]=nb_dof;
	  dof_to_elt.push_back(std::vector<N2>());
	  nb_dof++;
	}
	elt_to_dof[j][k] = num[I[k]];
	dof_to_elt[ num[I[k]] ].push_back( N2_(j,k) );
      }
    }
  }
  // special constructor for FreeFEM dof numbering
  Dof<ShapeFct>(const typename Trait::mesh_t& m, bool numToFF):
  mesh_p(&m), node_p(&GeometryOf(m)), nb_elt(NbElt(m)), offset(0) {
    
    const int nb_node = NbNode(*node_p);
    const typename Trait::mesh_t&  mesh = *mesh_p;

    nb_dof = nb_node;
    elt_to_dof.resize(NbElt(mesh));
    dof_to_elt.resize(nb_node);
    std::vector<int> num(nb_node,-1);
    for(int j=0; j<NbElt(mesh); j++){
      array<Trait::dim+1,int> I = Num(mesh[j]);
      for(int k=0; k<Trait::dim+1; k++){
		elt_to_dof[j][k] = I[k] ;
	    dof_to_elt[ I[k] ].push_back( N2_(j,k) );
       }
    }
	
  }

  typename Trait::R3loc operator()(const int& j) const {
    Elt<Trait::dim> e = (*mesh_p)[j];
    typename Trait::R3loc x;
    for(int j=0; j<Trait::dim+1; j++){x[j]=e[j];}
    return x;
  }

  const typename Trait::Nloc& operator[](const int& j) const {
    return elt_to_dof[j];}

  const typename Trait::Nloc&  operator[](const Elt<Trait::dim>& e) const {
    return (*this)[ (*node_p)[e][(*mesh_p)] ];}

  friend const int& NbDof(const this_t& d){return d.nb_dof;}
  friend const int& NbElt(const this_t& d){return d.nb_elt;}

  void operator+=(const int& offset0){
    offset += offset0;
    for(int j=0; j<nb_elt; j++){elt_to_dof[j]+=offset0;}}

  const std::vector<N2>& ToElt(const int& j) const {
    assert(j>=offset && j<offset+nb_dof);
    return dof_to_elt[j-offset];}


    const Elt<Trait::dim>& get_elt(int j) const{
      return (*mesh_p)[j];
    }

  friend const typename Trait::mesh_t&
  MeshOf(const this_t& dof){return *(dof.mesh_p);}
};


/*=========================
      Numerotation P2
  ========================*/

template <int Dim>
class Dof< BasisFct<P2,Dim> >{

public:
  typedef BasisFct<P2,Dim>     ShapeFct;
  typedef Dof<ShapeFct>        this_t;
  typedef DofTraits<ShapeFct>  Trait;

private:
  const Geometry*                  node_p;
  const typename Trait::mesh_t*    mesh_p;
  std::vector<typename Trait::Nloc>     elt_to_dof;
  std::vector< std::vector<N2> >             dof_to_elt;
  int                              nb_elt;
  int                              nb_dof;
  int                              offset;

public:
  Dof<ShapeFct>(): mesh_p(0), node_p(0), offset(0) {};

  Dof<ShapeFct>(const typename Trait::mesh_t& m):
  mesh_p(&m), node_p(&GeometryOf(m)), nb_elt(NbElt(m)), offset(0) {

    const int dim = Trait::dim;
    const int nb_edge_loc = Trait::nb_edge_loc;
    const int nb_node = NbNode(*node_p);
    const typename Trait::mesh_t&  mesh = *mesh_p;
    nb_dof = 0;
    elt_to_dof.resize(NbElt(mesh));

    int end = -1;
    std::vector<int>   num  (nb_node,-1);
    std::vector<int>   begin(nb_node,-1);
    std::vector<int>   next;
    std::vector<Elt1D> edge;
    int nb_edge=0;

    //numerotation dofs associes aux noeuds
    for(int j=0; j<NbElt(mesh); j++){
      array<dim+1,int> I = Num(mesh[j]);
      for(int k=0; k<dim+1; k++){
	if(num[I[k]]==-1){
	  num[I[k]]=nb_dof;
	  dof_to_elt.push_back(std::vector<N2>());
	  nb_dof++;}
	elt_to_dof[j][k] = num[I[k]];
	dof_to_elt[num[I[k]]].push_back( N2_(j,k) );
      }
    }

    int nb_nodal_dofs = nb_dof;
    //numerotation dofs associes aux aretes
    for(int j=0; j<NbElt(mesh); j++){
      array<nb_edge_loc,Elt1D> edge_loc = EdgesOf(mesh[j]);
      for(int k=0; k<nb_edge_loc; k++){
	bool exists = false;
	int  num_edge;
	Elt1D e = edge_loc[k]; Order(e);
	for(int p=begin[Key(e)]; p!=end; p=next[p]){
	  if(e==edge[p]){
	    exists=true;
	    int num_dof = nb_nodal_dofs+p;
	    elt_to_dof[j][dim+1+k] = num_dof;
	    dof_to_elt[num_dof].push_back( N2_(j,dim+1+k) );
	    break;
	  }
	}
	if(!exists){
	  next.push_back(begin[Key(e)]);
	  begin[Key(e)]=nb_edge;
	  elt_to_dof[j][dim+1+k]=nb_nodal_dofs+nb_edge;
	  edge.push_back(e);

	  dof_to_elt.push_back(std::vector<N2>());
	  dof_to_elt[nb_dof].push_back( N2_(j,dim+1+k) );

	  nb_edge++;
	  nb_dof++;
	}
      }
    }
  }


  typename Trait::R3loc operator()(const int&) const;

  const typename Trait::Nloc& operator[](const int& j) const {
    return elt_to_dof[j];}

  const typename Trait::Nloc&
  operator[](const Elt<Trait::dim>& e) const {
    return (*this)[ (*node_p)[e][(*mesh_p)] ];}

  friend const int& NbDof(const this_t& d){return d.nb_dof;}
  friend const int& NbElt(const this_t& d){return d.nb_elt;}

  void operator+=(const int& offset0){
    offset+=offset0;
    for(int j=0; j<nb_elt; j++){elt_to_dof[j]+=offset0;}}

  const std::vector<N2>& ToElt(const int& j) const {
    assert(j>=offset && j<offset+nb_dof);
    return dof_to_elt[j-offset];}

  friend const typename Trait::mesh_t&
  MeshOf(const this_t& dof){return *(dof.mesh_p);}

};


template<>
Dof<P2_1D>::Trait::R3loc
Dof<P2_1D>::operator()(const int& j) const {
  Elt<Trait::dim> e = (*mesh_p)[j];
   Dof<P2_1D>::Trait::R3loc x;
  x[0] = e[0];
  x[1] = e[1];
  x[2] = e[2];
  return x;
}


template<>
 Dof<P2_2D>::Trait::R3loc
 Dof<P2_2D>::operator()(const int& j) const {
  Elt<Trait::dim> e = (*mesh_p)[j];
   Dof<P2_2D>::Trait::R3loc x;
  x[0] = e[0];
  x[1] = e[1];
  x[2] = e[2];
  x[3] = 0.5*(x[1]+x[2]);
  x[4] = 0.5*(x[2]+x[0]);
  x[5] = 0.5*(x[0]+x[1]);
  return x;
}



/*================
  NUMEROTATION RT0
  ================*/

template <>
class Dof<RT0_2D>{

public:
  typedef RT0_2D               ShapeFct;
  typedef Dof<ShapeFct>        this_t;
  typedef DofTraits<ShapeFct>  Trait;

private:
  const Geometry*                node_p;
  const Trait::mesh_t*           mesh_p;
  std::vector<Trait::Nloc>       elt_to_dof;
  std::vector<std::vector<N2> >  dof_to_elt;
  int                            nb_elt;
  int                            nb_dof;
  int                            offset;

public:
  Dof<ShapeFct>(): mesh_p(0), node_p(0), offset(0) {};

  Dof<ShapeFct>(const Trait::mesh_t& m):
  mesh_p(&m), node_p(&GeometryOf(m)), nb_elt(NbElt(m)), offset(0) {

    const int dim = Trait::dim;
    const int nb_edge_loc = Trait::nb_edge_loc;
    const int nb_node = NbNode(*node_p);
    const Trait::mesh_t&  mesh = *mesh_p;
    nb_dof = 0;
    elt_to_dof.resize(NbElt(mesh));

    int end = -1;
    std::vector<int>   num  (nb_node,-1);
    std::vector<int>   begin(nb_node,-1);
    std::vector<int>   next;
    std::vector<Elt1D> edge;
    int nb_edge=0;

    //numerotation dofs associes aux aretes
    for(int j=0; j<NbElt(mesh); j++){
      array<nb_edge_loc,Elt1D> edge_loc = EdgesOf(mesh[j]);
      for(int k=0; k<nb_edge_loc; k++){
	bool exists = false;
	int  num_edge;
	Elt1D e = edge_loc[k]; Order(e);
	for(int p=begin[Key(e)]; p!=end; p=next[p]){
	  if(e==edge[p]){
	    exists=true;
	    elt_to_dof[j][k] = p;
	    dof_to_elt[p].push_back( N2_(j,k) );
	    break;
	  }
	}
	if(!exists){
	  next.push_back(begin[Key(e)]);
	  begin[Key(e)]=nb_edge;
	  elt_to_dof[j][k]=nb_edge;
	  edge.push_back(e);

	  dof_to_elt.push_back(std::vector<N2>());
	  dof_to_elt[nb_dof].push_back( N2_(j,k) );

	  nb_edge++;
	  nb_dof++;
	}
      }
    }

  }

  Trait::R3loc operator()(const int& j) const {
    Elt2D e = (*mesh_p)[j];
    Dof<RT0_2D>::Trait::R3loc x;
    x[0] = 0.5*(e[1]+e[2]);
    x[1] = 0.5*(e[2]+e[0]);
    x[2] = 0.5*(e[0]+e[1]);
    return x;
  }

  const  Trait::Nloc& operator[](const int& j) const {
    return elt_to_dof[j];}

  const  Trait::Nloc&
  operator[](const Elt<Trait::dim>& e) const {
    return (*this)[ (*node_p)[e][(*mesh_p)] ];}

  friend const int& NbDof(const this_t& d){return d.nb_dof;}
  friend const int& NbElt(const this_t& d){return d.nb_elt;}

  void operator+=(const int& offset0){
    offset+=offset0;
    for(int j=0; j<nb_elt; j++){elt_to_dof[j]+=offset0;}}

  const std::vector<N2>& ToElt(const int& j) const {
    assert(j>=offset && j<offset+nb_dof);
    return dof_to_elt[j-offset];}

  friend const  Trait::mesh_t&
  MeshOf(const this_t& dof){return *(dof.mesh_p);}

};

}// namespace bemtool


#endif
