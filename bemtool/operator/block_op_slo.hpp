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
#ifndef BEMTOOL_SUBOPERATOR_SLO_HPP
#define BEMTOOL_SUBOPERATOR_SLO_HPP

#include "block_op_slo.hpp"
#include "operator_slo.hpp"
#include <iostream>
#include <map>
#include <vector>

namespace bemtool {

template <int D, typename PhiX, typename PhiY> class SubBIOp_slo {

  typedef BlockMat  MatType;
  typedef Dof<PhiX> DofX;
  typedef Dof<PhiY> DofY;

  static const int nb_dof_x = PhiX::nb_dof_loc;
  static const int nb_dof_y = PhiY::nb_dof_loc;
  typedef NDofLoc<nb_dof_x> Nlocx;
  typedef NDofLoc<nb_dof_y> Nlocy;
  typedef typename std::map<int, Nlocx>::iterator ItTypeX;
  typedef typename std::map<int, Nlocy>::iterator ItTypeY;

private:
  BIOp_SLO<D,PhiX,PhiY> biop;
  const DofX &dofx;
  const DofY &dofy;

  BlockMat block_mat;
  MatType elt_mat;
  std::map<int, Nlocx> Ix;
  std::map<int, Nlocy> Iy;

public:
  SubBIOp_slo(const DofX &dofx0, const DofY &dofy0)
      : biop(MeshOf(dofx0), MeshOf(dofy0)), dofx(dofx0), dofy(dofy0) {}

  const BlockMat &operator()(const std::vector<int> &jjx,
                             const std::vector<int> &jjy) {
    block_mat.Resize(jjx.size(), jjy.size());
    elt_mat = 0.;
    Ix.clear();
    Iy.clear();

    for (int k = 0; k < jjx.size(); k++) {
      const std::vector<N2> &jj = dofx.ToElt(jjx[k]);
      for (int l = 0; l < jj.size(); l++) {
        const N2 &j = jj[l];
        Ix[j[0]][j[1]] = k;
      }
    }

    for (int k = 0; k < jjy.size(); k++) {
      const std::vector<N2> &jj = dofy.ToElt(jjy[k]);
      for (int l = 0; l < jj.size(); l++) {
        const N2 &j = jj[l];
        Iy[j[0]][j[1]] = k;
      }
    }

    for (ItTypeX itx = Ix.begin(); itx != Ix.end(); itx++) {
      const int &jx = itx->first;
      const Nlocx &nx = itx->second;

      for (ItTypeY ity = Iy.begin(); ity != Iy.end(); ity++) {
        const int &jy = ity->first;
        const Nlocy &ny = ity->second;

        elt_mat = biop(jx, jy);

        // Adapt to the size and numbering of elt_mat
        std::vector<int> slo_num = biop.get_slo_num();
        std::vector<int> I(NbRow(elt_mat),0);

        for (int l=0;l<size(nx);l++){
          I[slo_num[l]]=nx[l];
          I[slo_num[l+size(nx)]]=ny[l];
        } 


        for (int kx = 0; kx < I.size(); kx++) {
          if (I[kx] != -1) {
            for (int ky = 0; ky < I.size(); ky++) {
              if (I[ky] != -1) {
                block_mat(I[kx], I[ky]) += elt_mat(kx, ky);
              }
            }
          }
        }
      }
    }

    return block_mat;
  }

  template <typename Matrix>
  void compute_block(const std::vector<int> &jjx, const std::vector<int> &jjy,
                     Matrix &mat) {
    elt_mat.Clear();
    Ix.clear();
    Iy.clear();

    for (int k = 0; k < jjx.size(); k++) {
      const std::vector<N2> &jj = dofx.ToElt(jjx[k]);
      for (int l = 0; l < jj.size(); l++) {
        const N2 &j = jj[l];
        Ix[j[0]][j[1]] = k;
      }
    }

    for (int k = 0; k < jjy.size(); k++) {
      const std::vector<N2> &jj = dofy.ToElt(jjy[k]);
      for (int l = 0; l < jj.size(); l++) {
        const N2 &j = jj[l];
        Iy[j[0]][j[1]] = k;
      }
    }
  std::cout << Ix.size()<<" "<<Iy.size()<<std::endl;
    for (ItTypeX itx = Ix.begin(); itx != Ix.end(); itx++) {
      const int &jx = itx->first;
      const Nlocx &nx = itx->second;
      std::cout << "jx : "<<jx<<" "<<nx[0]<<" "<<nx[1]<<std::endl;



      for (ItTypeY ity = Iy.begin(); ity != Iy.end(); ity++) {
        std::cout << "size_x "<<Ix.size()<<" size_z "<<Iy.size()<<std::endl;
        const int &jy = ity->first;
        const Nlocy &ny = ity->second;
        std::cout << "jy : "<< jy<<std::endl;
        elt_mat = biop(jx, jy);
        std::cout << "biop"<<std::endl;
        // Adapt to the size and numbering of elt_mat
        std::vector<int> slo_num = biop.get_slo_num();
        std::cout << "slo_num ";
        for (int i=0;i<slo_num.size();i++){
          std::cout <<slo_num[i]<<" ";
        }
        std::cout << std::endl;
        std::vector<int> I(NbRow(elt_mat),0);
        std::cout << "I ";
        for (int l=0;l<size(nx);l++){
          I[slo_num[l]]=nx[l];
          I[slo_num[l+size(nx)]]=ny[l];
        } 
        for (int i=0;i<I.size();i++){
          std::cout <<I[i]<<" ";
        }
        std::cout << std::endl;
        for (int kx = 0; kx < I.size(); kx++) {
          if (I[kx] != -1) {
            for (int ky = 0; ky < I.size(); ky++) {
              if (I[ky] != -1) {
                std::cout <<mat.get_ir().size()<<" "<<mat.get_ic().size()<<std::endl;
                std::cout << "TEST "<<kx<<" "<<I[kx]<<" "<<ky<<" "<<I[ky]<<std::endl;
                mat(I[kx], I[ky]) += elt_mat(kx, ky);
              }
            }
          }
        }
        std::cout << "mat"<<std::endl;
      }
    }
  }

  template <typename Matrix>
  void compute_neumann_block(const std::vector<int> &jjx,
                             const std::vector<int> &jjy, Matrix &mat) {
    elt_mat.Clear();
    Ix.clear();
    Iy.clear();

    for (int k = 0; k < jjx.size(); k++) {
      const std::vector<N2> &jj = dofx.ToElt(jjx[k]);
      for (int l = 0; l < jj.size(); l++) {
        const N2 &j = jj[l];
        Ix[j[0]][j[1]] = k;
      }
    }

    for (int k = 0; k < jjy.size(); k++) {
      const std::vector<N2> &jj = dofy.ToElt(jjy[k]);
      for (int l = 0; l < jj.size(); l++) {
        const N2 &j = jj[l];
        Iy[j[0]][j[1]] = k;
      }
    }

    for (ItTypeX itx = Ix.begin(); itx != Ix.end(); itx++) {
      const int &jx = itx->first;
      const Nlocx &nx = itx->second;

      for (ItTypeY ity = Iy.begin(); ity != Iy.end(); ity++) {
        const int &jy = ity->first;
        const Nlocy &ny = ity->second;

        bool testx = 1;
        bool testy = 1;
        for (int kx = 0; kx < nb_dof_x; kx++) {
          if (nx[kx] == -1) {
            testx = 0;
          }
        }
        for (int ky = 0; ky < nb_dof_y; ky++) {
          if (ny[ky] == -1) {
            testy = 0;
          }
        }
        if (testx && testy) {
          elt_mat = biop(jx, jy);
          // Adapt to the size and numbering of elt_mat
          std::vector<int> slo_num = biop.get_slo_num();
          std::vector<int> I(NbRow(elt_mat),0);

          for (int l=0;l<size(nx);l++){
            I[slo_num[l]]=nx[l];
            I[slo_num[l+size(nx)]]=ny[l];
          } 
          for (int kx = 0; kx < I.size(); kx++) {
            for (int ky = 0; ky < I.size(); ky++) {
              mat(I[kx], I[ky]) += elt_mat(kx, ky);
            }
          }
        }
      }
    }
  }

//   template <typename Matrix>
//   void compute_neumann_block(const std::vector<int> &jjx,
//                              const std::vector<int> &jjy, Matrix &mat,
//                              const std::map<int, Nlocx> &Ix_g,
//                              const std::map<int, Nlocx> &Iy_g) {
//     elt_mat = 0.;
//     Ix.clear();
//     Iy.clear();

//     for (int k = 0; k < jjx.size(); k++) {
//       const std::vector<N2> &jj = dofx.ToElt(jjx[k]);
//       for (int l = 0; l < jj.size(); l++) {
//         const N2 &j = jj[l];
//         Ix[j[0]][j[1]] = k;
//       }
//     }

//     for (int k = 0; k < jjy.size(); k++) {
//       const std::vector<N2> &jj = dofy.ToElt(jjy[k]);
//       for (int l = 0; l < jj.size(); l++) {
//         const N2 &j = jj[l];
//         Iy[j[0]][j[1]] = k;
//       }
//     }

//     for (auto itx = Ix.begin(); itx != Ix.end(); itx++) {
//       const int &jx = itx->first;
//       const Nlocx &nx = itx->second;
//       const Nlocx &nx_g = Ix_g.at(jx);

//       for (auto ity = Iy.begin(); ity != Iy.end(); ity++) {
//         const int &jy = ity->first;
//         const Nlocy &ny = ity->second;
//         const Nlocy &ny_g = Iy_g.at(jy);

//         bool testx = 1;
//         bool testy = 1;
//         for (int kx = 0; kx < nb_dof_loc_x; kx++) {
//           if (nx_g[kx] == -1) {
//             testx = 0;
//           }
//         }
//         for (int ky = 0; ky < nb_dof_loc_y; ky++) {
//           if (ny_g[ky] == -1) {
//             testy = 0;
//           }
//         }
//         if (testx && testy) {
//           elt_mat = biop(jx, jy);

//           for (int kx = 0; kx < I.size(); kx++) {
//             if (nx[kx] != -1) {
//               for (int ky = 0; ky < I.size(); ky++) {
//                 if (ny[ky] != -1) {
//                   mat(nx[kx], ny[ky]) += elt_mat(kx, ky);
//                 }
//               }
//             }
//           }
//         }
//       }
//     }
//   }
};

} // namespace bemtool
#endif
