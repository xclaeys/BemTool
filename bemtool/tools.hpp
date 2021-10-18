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

#ifndef BEMTOOL_TOOLS_HPP
#define BEMTOOL_TOOLS_HPP

//=====================//

#include "equations.hpp"

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
#include "operator/operator_slo.hpp"
#include "operator/helmholtz_op.hpp"
#include "operator/laplace_op.hpp"
#include "operator/yukawa_op.hpp"
#include "operator/maxwell_op.hpp"
#include "operator/block_op.hpp"
#include "operator/block_op_slo.hpp"

#include "potential/potential.hpp"
#include "potential/helmholtz_pot.hpp"
#include "potential/laplace_pot.hpp"
#include "potential/yukawa_pot.hpp"
#include "potential/maxwell_pot.hpp"

#include "quadrature/dunavant.hpp"
#include "quadrature/quad.hpp"
#include "quadrature/quad_bem.hpp"
#include "quadrature/quad_pot.hpp"

#include "miscellaneous/output_gmsh.hpp"
#include "miscellaneous/output_paraview.hpp"
#include "miscellaneous/misc.hpp"
#include "miscellaneous/coordinates.hpp"
#include "miscellaneous/specialfct.hpp"
#include "miscellaneous/refeigenvalue.hpp"
#include "miscellaneous/refsol.hpp"
#include "miscellaneous/overlap.hpp"


#endif
