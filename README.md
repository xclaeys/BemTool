# BemTool
BemTool is a C++ header library that only relies on the STL and Boost::math library.
Its main purpose is the computation of assembly routines for boundary elements. At low level,
it contains a set of elementary classes for low dimensional linear algebra (so as to deal with
geometric computations) and simplicial meshes in 1D,2D and 3D. It contains basic wrappers
toward Gmsh (cf http://gmsh.info/) for mesh generation, and eigen (http://eigen.tuxfamily.org/)
for the treatment of large matrices.


Currently it provides routines for the assembly boundary element matrices for Laplace,
Helmholtz, Yukawa (purely dissipative Hekmholtz) and Maxwell matrices in 2D and 3D.
The available discretisation spaces are P0, continuous P1, continuous P2 and lowest order
Raviart-Thomas elements.
