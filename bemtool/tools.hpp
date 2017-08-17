#ifndef TOOLS_HPP
#define TOOLS_HPP


/*=====================
  Equations disponibles
  =====================*/

enum EquationEnum {
  CST, // noyau constant
  HE , // Helmholtz
  LA , // Laplace
  YU , // Yukawa
  MA   // Maxwell
};

/*=====================
  Type d'operateurs 
  =====================*/

enum BIOpKernelEnum {
  CST_OP,  // noyau constant
   SL_OP,  // trace dirichlet du simple couche
   DL_OP,  // trace dirichlet du double couche
  TDL_OP,  // trace neumann   du simple couche
   HS_OP   // trace neumann   du double couche
};

/*=====================
  Type de potentiels 
  =====================*/

enum PotKernelEnum {
  CST_POT, // noyau constant
   SL_POT, // simple couche
   DL_POT  // double couche
};



//=====================//

#include "calculus/calculus.hpp"

#include "mesh/element.hpp"
#include "mesh/adjacency.hpp"
#include "mesh/normal.hpp"
#include "mesh/mesh.hpp"

#include "fem/fem.hpp"
#include "fem/femP1.hpp"
#include "fem/dof.hpp"
#include "fem/shapefct.hpp"
#include "fem/interpolation.hpp"

#include "operator/operator.hpp"
#include "operator/helmholtz_op.hpp"
#include "operator/laplace_op.hpp"
#include "operator/yukawa_op.hpp"
#include "operator/maxwell_op.hpp"

#include "potential/potential.hpp"
#include "potential/helmholtz_pot.hpp"
#include "potential/laplace_pot.hpp"
#include "potential/yukawa_pot.hpp"

#include "quadrature/dunavant.hpp" 
#include "quadrature/quad.hpp"
#include "quadrature/quad_bem.hpp"
#include "quadrature/quad_pot.hpp"

#include "miscellaneous/eigen_wrap.hpp"
#include "miscellaneous/output.hpp"
#include "miscellaneous/misc.hpp"
#include "miscellaneous/coordinates.hpp"
#include "miscellaneous/specialfct.hpp"
#include "miscellaneous/refsol.hpp"

#endif
