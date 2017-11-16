#ifndef BEMTOOL_MISC_OVERLAP_HPP
#define BEMTOOL_MISC_OVERLAP_HPP

#include <set>
#include <mpi.h>
#include <algorithm>
#include <numeric>

template<class Dof>
void Partition(const std::vector<std::pair<int,int>>& MasterOffset, const std::vector<int>& perm, const Dof& dof, std::vector<int>& cluster_to_ovr_subdomain, std::vector<int>& ovr_subdomain_to_global, std::vector<int>& neighbors, std::vector<std::vector<int> >& intersections){

  // Get the number of processes
  int sizeWorld;
  MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);

  // Get the rank of the process
  int rankWorld;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankWorld);

  //
  int nbdof = NbDof(dof);
  int nbelt = NbElt(dof);
  std::vector<bool> part_overlap(nbdof);

  // Partitionnement des dofs sans overlap
  std::vector<int> part(nbdof);
  for (int i=0; i<sizeWorld; i++)
	for (int j=0; j< MasterOffset[i].second; j++)
		part[perm[MasterOffset[i].first+j]] = i;
    // part[MasterClusters[i][j]] = i;

  // Tag des elements qui contiennent mes dofs: P1 -> P0
  std::vector<int> elts(nbelt,0);
  for(int j=0; j<nbelt; j++){
    for(int k=0; k<Dof::Trait::dim+1; k++){
      elts[j]=elts[j]||(part[dof[j][k]]==rankWorld);
    }
  }

  // tag de mes dofs avec part_overlap
  for (int i =0;i<nbdof;i++){
    part_overlap[i]=(part[i]==rankWorld);
  }

	// Tag de mes dofs qui sont contenus dans les elements: P0 -> P1
  for(int j=0; j<nbelt; j++){
    for(int k=0; k<Dof::Trait::dim+1; k++){
      part_overlap[dof[j][k]]= part_overlap[dof[j][k]] || elts[j];
    }
  }

	// On va deux fois plus loin pour trouver les voisins à partir de part
	std::vector<bool> part_find_neighbors=part_overlap;
	for(int j=0; j<nbelt; j++){
		for(int k=0; k<Dof::Trait::dim+1; k++){
			elts[j]=elts[j]||(part_find_neighbors[dof[j][k]]);
		}
	}

	for(int j=0; j<nbelt; j++){
    for(int k=0; k<Dof::Trait::dim+1; k++){
      part_find_neighbors[dof[j][k]]= part_find_neighbors[dof[j][k]] || elts[j];
    }
  }

  std::set<int> neighbors_set;
  neighbors.clear();
  for (int i=0;i<nbdof;i++){
    if (part_find_neighbors[i]){
      if (part[i]!=rankWorld){
        neighbors_set.insert(part[i]);
      }
    }
  }
  for(auto iter=neighbors_set.begin(); iter!=neighbors_set.end();++iter) {
    neighbors.push_back(*iter);
  }


  // Get the global to overlapping subdomain numbering
  std::vector<int> temp(nbdof,-1);
  int count =0;
  for (int i =0;i<nbdof;i++){
    if (part_overlap[i]){
      temp[i]=count++;
    }
  }

  // Get the cluster to overlapping subdomain numbering
  cluster_to_ovr_subdomain.resize(MasterOffset[rankWorld].second);
  for (int j=0; j< MasterOffset[rankWorld].second; j++)
		cluster_to_ovr_subdomain[j] = temp[perm[MasterOffset[rankWorld].first+j]];

  // Get the overlapping subdomain to global numbering
  int size_ovr_subdomain = std::accumulate(part_overlap.begin(),part_overlap.end(),0,[](int a, bool b){ return a+(std::size_t)(b);});

  ovr_subdomain_to_global.resize(size_ovr_subdomain);
  for (int i=0; i< nbdof; i++){
    if (temp[i]!=-1){
      ovr_subdomain_to_global[temp[i]]=i;
    }
  }

	// Calcul de l'overlap des voisins
	intersections.clear();
  for (int i_n=0;i_n<neighbors.size();i_n++){
		std::vector<bool> part_overlap_neighbors(nbdof);

		// Tag des elements qui contiennent les dofs du voisin: P1 -> P0
	  fill(elts.begin(),elts.end(),0);
	  for(int j=0; j<nbelt; j++){
	    for(int k=0; k<Dof::Trait::dim+1; k++){
	      elts[j]=elts[j]||(part[dof[j][k]] == neighbors[i_n]);
	    }
	  }

		// tag des dofs du voisin avec part_overlap_neighbors
	  for (int i =0;i<part.size();i++){
	    part_overlap_neighbors[i]=(part[i] == neighbors[i_n]);
	  }

		// Tag des dofs du voisin qui sont contenus dans les elements: P0 -> P1
	  for(int j=0; j<nbelt; j++){
	    for(int k=0; k<Dof::Trait::dim+1; k++){
	      part_overlap_neighbors[dof[j][k]]= part_overlap_neighbors[dof[j][k]] || elts[j];
	    }
	  }

		std::transform(part_overlap_neighbors.begin(),part_overlap_neighbors.end(),part_overlap.begin(),part_overlap_neighbors.begin(),[](bool a,bool b){return a&&b;});

		std::size_t size_intersection = std::accumulate(part_overlap_neighbors.begin(),part_overlap_neighbors.end(),0,[](int a, bool b){return a+(std::size_t)(b);});

		std::vector<int> intersection;
		intersection.reserve(size_intersection);
		for (int i = 0; i < nbdof; i++) {
			if (part_overlap_neighbors[i]){
				intersection.push_back(temp[i]);
			}
		}
		intersections.push_back(intersection);
  }
  // MPI_Barrier(MPI_COMM_WORLD);
  // if (rankWorld==0){
  //   std::cout << size_ovr_subdomain << std::endl;
  //   std::cout << intersections[0].size() << std::endl;
  //   for (int i=0 ; i<intersections[0].size();i++){
  //     std::cout << intersections[0][i] << " ";
  //   }
  //   for (int i =0;i<neighbors.size();i++){
  //     std::cout << neighbors[i] << std::endl;
  //   }
  //   std::cout<<std::endl;
  //   std::cout<<std::endl;
  // }
  // MPI_Barrier(MPI_COMM_WORLD);
  // if (rankWorld==1){
  //   std::cout << size_ovr_subdomain << std::endl;
  //   std::cout << intersections[0].size() << std::endl;
  //   for (int i=0 ; i<intersections[0].size();i++){
  //     std::cout << intersections[0][i] << " ";
  //   }
  //   for (int i =0;i<neighbors.size();i++){
  //     std::cout << neighbors[i] << std::endl;
  //   }
  //   std::cout<<std::endl;
  //   std::cout<<std::endl;
  // }
  // MPI_Barrier(MPI_COMM_WORLD);
}




#endif
