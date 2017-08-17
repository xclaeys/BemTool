#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP
#include "dof.hpp"



template<typename> class Interpolator;

template<int D> 
class Interpolator< BasisFct<P0,D> >{

  typedef BasisFct<P0,D>  ShapeFct;
  typedef array<1,Cplx>   ReturnType;
  
 private:
  const Mesh<D>&  mesh;

 public:
  Interpolator< BasisFct<P0,D> >(const Mesh<D>& m): mesh(m) {};
  
  template <typename fct_t>
    ReturnType operator()(const fct_t& f, const int& j){
    ReturnType res; res[0] = f( Ctr(mesh[j]) );
    return res; }
  
};



template<int D>
class Interpolator< BasisFct<P1,D> >{

  typedef BasisFct<P1,D>  ShapeFct;
  typedef array<D+1,Cplx> ReturnType;
  
 private:
  const Mesh<D>&  mesh;

 public:
  Interpolator< BasisFct<P1,D> >(const Mesh<D>& m): mesh(m) {};
  
  template <typename fct_t>
    ReturnType operator()(const fct_t& f, const int& j){
    ReturnType res;
    for(int k=0; k<D+1; k++){
      res[k] = f(mesh[j][k]);}
    return res;}
  
};



template<>
class Interpolator<RT0_2D>{
  
  typedef RT0_2D        ShapeFct;
  typedef array<3,Cplx> ReturnType;
  
 private:
  const Mesh2D&  mesh;
  
 public:
  Interpolator<RT0_2D>(const Mesh2D& m): mesh(m) {}
  
  template <typename fct_t>
    ReturnType operator()(const fct_t& f, const int& j){
    ReturnType res; const Elt2D& e = mesh[j];
    array<3,R3> n_;
    NormalToFaces(mesh[j],n_);
    array<3,Elt1D> edge = EdgesOf(e);    
    for(int k=0; k<3; k++){
      R3   x  = Ctr(edge[k]);
      Real l  = Vol(edge[k]);
      res[k]  = l*(n_[k],f(x));
    }
    return res;}
  
  
  
};





#endif
