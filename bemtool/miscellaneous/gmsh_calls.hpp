#ifndef GMSH_CALL_H
#define GMSH_CALL_H

#include <string>
#include <fstream>
#include "../calculus/calculus.hpp"
#include "misc.hpp"

namespace bemtool{



////=============================================================////
////===========================  Circle =========================////
////=============================================================////
inline void gmsh_circle(std::string mesh_name, Real R, Real lc, int verbose=0){
	std::ofstream circle((mesh_name+".geo").c_str());
	if (circle.is_open()) {
		circle << "lc = "+NbrToStr(lc)+";\n";
		circle << "Point(0) = { 0 , 0 , 0 , lc};\n";

		//// Droite
		circle << "Point(1) = { "+NbrToStr(R)+" , "+NbrToStr(0)+" , "+NbrToStr(0)+" , lc}; \n";
		//// Haut
		circle << "Point(2) = { "+NbrToStr(0)+" , "+NbrToStr(R)+" , "+NbrToStr(0)+" , lc}; \n";
		//// Gauche
		circle << "Point(3) = { "+NbrToStr(-R)+" , "+NbrToStr(0)+" , "+NbrToStr(0)+" , lc}; \n";
		//// Bas
		circle << "Point(4) = { "+NbrToStr(0)+" , "+NbrToStr(-R)+" , "+NbrToStr(0)+" , lc}; \n";

		//// Droite
		circle << "Circle(1) = { "+NbrToStr(1)+" , "+NbrToStr(0)+" , "+NbrToStr(2)+"}; \n";
		//// Haut
		circle << "Circle(2) = { "+NbrToStr(2)+" , "+NbrToStr(0)+" , "+NbrToStr(3)+"}; \n";
		//// Gauche
		circle << "Circle(3) = { "+NbrToStr(3)+" , "+NbrToStr(0)+" , "+NbrToStr(4)+"}; \n";
		//// Bas
		circle << "Circle(4) = { "+NbrToStr(4)+" , "+NbrToStr(0)+" , "+NbrToStr(1)+"}; \n";


		circle << "Physical Line(0) ={ "+NbrToStr(1)+" , "+NbrToStr(2)+" , "+NbrToStr(3)+" , "+NbrToStr(4)+"};\n";


		circle.close();
	}
	else std::cout << "Unable to open file \n" <<std::endl;

	system(("gmsh -2 -v "+NbrToStr(verbose)+" "+mesh_name+".geo").c_str());
}

////=============================================================////
////===========================  Disc =========================////
////=============================================================////
inline void gmsh_disc(std::string mesh_name, Real R, Real lc, int verbose=0){
	std::ofstream circle((mesh_name+".geo").c_str());
	if (circle.is_open()) {
		circle << "lc = "+NbrToStr(lc)+";\n";
		circle << "Point(0) = { 0 , 0 , 0 , lc};\n";

		//// Droite
		circle << "Point(1) = { "+NbrToStr(R)+" , "+NbrToStr(0)+" , "+NbrToStr(0)+" , lc}; \n";
		//// Haut
		circle << "Point(2) = { "+NbrToStr(0)+" , "+NbrToStr(R)+" , "+NbrToStr(0)+" , lc}; \n";
		//// Gauche
		circle << "Point(3) = { "+NbrToStr(-R)+" , "+NbrToStr(0)+" , "+NbrToStr(0)+" , lc}; \n";
		//// Bas
		circle << "Point(4) = { "+NbrToStr(0)+" , "+NbrToStr(-R)+" , "+NbrToStr(0)+" , lc}; \n";

		//// Droite
		circle << "Circle(1) = { "+NbrToStr(1)+" , "+NbrToStr(0)+" , "+NbrToStr(2)+"}; \n";
		//// Haut
		circle << "Circle(2) = { "+NbrToStr(2)+" , "+NbrToStr(0)+" , "+NbrToStr(3)+"}; \n";
		//// Gauche
		circle << "Circle(3) = { "+NbrToStr(3)+" , "+NbrToStr(0)+" , "+NbrToStr(4)+"}; \n";
		//// Bas
		circle << "Circle(4) = { "+NbrToStr(4)+" , "+NbrToStr(0)+" , "+NbrToStr(1)+"}; \n";


		circle << "Line Loop(0) ={ "+NbrToStr(1)+" , "+NbrToStr(2)+" , "+NbrToStr(3)+" , "+NbrToStr(4)+"};\n";
		circle << "Plane Surface(0) = {"+NbrToStr(0)+"};\n";
		circle << "Physical Surface(0) = {"+NbrToStr(0)+"};\n";

		circle.close();
	}
	else std::cout << "Unable to open file \n" <<std::endl;

	system(("gmsh -2 -v "+NbrToStr(verbose)+" "+mesh_name+".geo").c_str());
}

////=============================================================////
////===========================  Sphere =========================////
////=============================================================////
inline void gmsh_sphere(std::string mesh_name, Real R, Real lc, int verbose=0){
	std::ofstream sphere((mesh_name+".geo").c_str());
	if (sphere.is_open()) {
		sphere << "lc = "+NbrToStr(lc)+";\n";

		sphere << "Point(1) = { 0 , 0 , 0 , lc};\n";

		sphere << "Point(2) = { "+NbrToStr(R)+" , "+NbrToStr(0)+" , "+NbrToStr(0)+" , lc}; \n";
		sphere << "Point(3) = { "+NbrToStr(0)+" , "+NbrToStr(R)+" , "+NbrToStr(0)+" , lc}; \n";
		sphere << "Point(4) = { "+NbrToStr(0)+" , "+NbrToStr(0)+" , "+NbrToStr(R)+" , lc}; \n";
		sphere << "Point(5) = { "+NbrToStr(-R)+" , "+NbrToStr(0)+" , "+NbrToStr(0)+" , lc}; \n";
		sphere << "Point(6) = { "+NbrToStr(0)+" , "+NbrToStr(-R)+" , "+NbrToStr(0)+" , lc}; \n";
		sphere << "Point(7) = { "+NbrToStr(0)+" , "+NbrToStr(0)+" , "+NbrToStr(-R)+" , lc}; \n";

		sphere << "Circle(1) = {2,1,3}; \n";
		sphere << "Circle(2) = {3,1,5}; \n";
		sphere << "Circle(3) = {5,1,6}; \n";
		sphere << "Circle(4) = {6,1,2}; \n";
		sphere << "Circle(5) = {2,1,7}; \n";
		sphere << "Circle(6) = {7,1,5}; \n";
		sphere << "Circle(7) = {5,1,4}; \n";
		sphere << "Circle(8) = {4,1,2}; \n";
		sphere << "Circle(9) = {6,1,7}; \n";
		sphere << "Circle(10) = {7,1,3}; \n";
		sphere << "Circle(11) = {3,1,4}; \n";
		sphere << "Circle(12) ={4,1,6}; \n";


		sphere << "Line Loop(1) = {1,11,8}; \n";
		sphere << "Line Loop(2) = {2,7,-11};  \n";
		sphere << "Line Loop(3) = {3,-12,-7};  \n";
		sphere << "Line Loop(4) = {4,-8,12};  \n";
		sphere << "Line Loop(5) = {5,10,-1};  \n";
		sphere << "Line Loop(6) = {-2,-10,6}; \n";
		sphere << "Line Loop(7) = {-3,-6,-9};  \n";
		sphere << "Line Loop(8) = {-4,9,-5};  \n";


		sphere << "Ruled Surface(1) = {1};\n";
		sphere << "Ruled Surface(2) = {2};\n";
		sphere << "Ruled Surface(3) = {3};\n";
		sphere << "Ruled Surface(4) = {4};\n";
		sphere << "Ruled Surface(5) = {5};\n";
		sphere << "Ruled Surface(6) = {6};\n";
		sphere << "Ruled Surface(7) = {7};\n";
		sphere << "Ruled Surface(8) = {8};\n";


		sphere << "Surface Loop (1) = {1,2,3,4,5,6,7,8};\n";
		sphere << "Volume (1) = {1};\n";
		sphere.close();
	}
	else std::cout << "Unable to open file \n" <<std::endl;

	system(("gmsh -2 -v "+NbrToStr(verbose)+" "+mesh_name+".geo").c_str());
}

////=============================================================////
////===========================  Ball ===========================////
////=============================================================////
inline void gmsh_ball(std::string mesh_name, Real R, Real lc, int verbose=0){
	std::ofstream ball((mesh_name+".geo").c_str());
	if (ball.is_open()) {
		ball << "lc = "+NbrToStr(lc)+";\n";

		ball << "Point(1) = { 0 , 0 , 0 , lc};\n";

		ball << "Point(2) = { "+NbrToStr(R)+" , "+NbrToStr(0)+" , "+NbrToStr(0)+" , lc}; \n";
		ball << "Point(3) = { "+NbrToStr(0)+" , "+NbrToStr(R)+" , "+NbrToStr(0)+" , lc}; \n";
		ball << "Point(4) = { "+NbrToStr(0)+" , "+NbrToStr(0)+" , "+NbrToStr(R)+" , lc}; \n";
		ball << "Point(5) = { "+NbrToStr(-R)+" , "+NbrToStr(0)+" , "+NbrToStr(0)+" , lc}; \n";
		ball << "Point(6) = { "+NbrToStr(0)+" , "+NbrToStr(-R)+" , "+NbrToStr(0)+" , lc}; \n";
		ball << "Point(7) = { "+NbrToStr(0)+" , "+NbrToStr(0)+" , "+NbrToStr(-R)+" , lc}; \n";

		ball << "Circle(1) = {2,1,3}; \n";
		ball << "Circle(2) = {3,1,5}; \n";
		ball << "Circle(3) = {5,1,6}; \n";
		ball << "Circle(4) = {6,1,2}; \n";
		ball << "Circle(5) = {2,1,7}; \n";
		ball << "Circle(6) = {7,1,5}; \n";
		ball << "Circle(7) = {5,1,4}; \n";
		ball << "Circle(8) = {4,1,2}; \n";
		ball << "Circle(9) = {6,1,7}; \n";
		ball << "Circle(10) = {7,1,3}; \n";
		ball << "Circle(11) = {3,1,4}; \n";
		ball << "Circle(12) ={4,1,6}; \n";


		ball << "Line Loop(1) = {1,11,8}; \n";
		ball << "Line Loop(2) = {2,7,-11};  \n";
		ball << "Line Loop(3) = {3,-12,-7};  \n";
		ball << "Line Loop(4) = {4,-8,12};  \n";
		ball << "Line Loop(5) = {5,10,-1};  \n";
		ball << "Line Loop(6) = {-2,-10,6}; \n";
		ball << "Line Loop(7) = {-3,-6,-9};  \n";
		ball << "Line Loop(8) = {-4,9,-5};  \n";


		ball << "Ruled Surface(1) = {1};\n";
		ball << "Ruled Surface(2) = {2};\n";
		ball << "Ruled Surface(3) = {3};\n";
		ball << "Ruled Surface(4) = {4};\n";
		ball << "Ruled Surface(5) = {5};\n";
		ball << "Ruled Surface(6) = {6};\n";
		ball << "Ruled Surface(7) = {7};\n";
		ball << "Ruled Surface(8) = {8};\n";


		ball << "Surface Loop (1) = {1,2,3,4,5,6,7,8};\n";
		ball << "Volume (1) = {1};\n";
		ball.close();
	}
	else std::cout << "Unable to open file \n" <<std::endl;

	system(("gmsh -3 -v "+NbrToStr(verbose)+" "+mesh_name+".geo").c_str());
}

////=============================================================////
////==========================  Segment =========================////
////=============================================================////
inline void gmsh_seg(std::string mesh_name, R3 A, R3 B, int verbose=0){
	std::ofstream segment((mesh_name+".geo").c_str());
	if (segment.is_open()) {
		Real lc =0.1;
		segment << "lc = "+NbrToStr(lc)+";\n";

		//// Droite
		segment << "Point(0) = { "+NbrToStr(A[0])+" , "+NbrToStr(A[1])+" , "+NbrToStr(A[2])+" , lc}; \n";
		//// Haut
		segment << "Point(1) = { "+NbrToStr(B[0])+" , "+NbrToStr(B[1])+" , "+NbrToStr(B[2])+" , lc}; \n";

		segment << "Line(2) =  {0,1}; \n";

		segment << "Physical Line(3) ={3};\n";


		segment.close();
	}
	else std::cout << "Unable to open file \n" <<std::endl;

	system(("gmsh -2 -v "+NbrToStr(verbose)+" "+mesh_name+".geo").c_str());
}

////=============================================================////
////===========================  Clean ==========================////
////=============================================================////
inline void gmsh_clean(std::string mesh_name){
	system(("rm "+mesh_name+".geo").c_str());
	system(("rm "+mesh_name+".msh").c_str());
}

////=============================================================////
////==========================  Plot  ===========================////
////=============================================================////
template <class m_t, class f_t> void write_gmsh(const m_t & m, const f_t & f, std::string filename, std::string name="\"\""){

	const vect<R3>& node = get_node(get_geometry(m));
	int nbnode = size(node);
	int nbelt = nb_elt(m);

	std::ofstream file; file.open((filename+".msh").c_str());
	file << "$MeshFormat\n";
	file << "2.2 0 8\n";
	file << "$EndMeshFormat\n";

	//== Nodes
	file << "$Nodes\n";
	file << nbnode << std::endl;
	for(int j=0; j<nbnode; j++){
		file << j+1 <<" "<< node[j][0]<<" "<<node[j][1]<<" "<<node[j][2] << "\n";
	}
	file << "$EndNodes\n";

	//== Elements
	file << "$Elements\n";
	file << nbelt << std::endl;
	for (int j=0;j<nbelt;j++){
		if (get_dim(m[j])==1){
			file << j+1 <<" " <<1 <<" "<<2<<" "<<0<<" "<<0<<" "; // number type_elt (ligne) nbr_tags
			file << &m[j][0]-&node[0] +1<<" ";
			file << &m[j][1]-&node[0] +1<<" "<<std::endl;
		}
		if (get_dim(m[j])==2){
			file << j+1 <<" " <<2 <<" "<<2<<" "<<0<<" "<<0<<" "; // number type_elt (triangle) nbr_tags
			file << &m[j][0]-&node[0] +1<<" ";
			file << &m[j][1]-&node[0] +1<<" ";
			file << &m[j][2]-&node[0] +1<<" "<<std::endl;
		}
		if (get_dim(m[j])==3){
			file << j+1 <<" " <<4 <<" "<<2<<" "<<0<<" "<<0<<" "; // number type_elt (tetrahedron) nbr_tags
			file << &m[j][0]-&node[0] +1<<" ";
			file << &m[j][1]-&node[0] +1<<" ";
			file << &m[j][2]-&node[0] +1<<" ";
			file << &m[j][3]-&node[0] +1<<" "<<std::endl;
		}



	}
	file << "$EndElements\n";

	//== Nodedata
	file << "$NodeData\n";
	file << 1<<std::endl; // number-of-string-tags
	file <<name<<std::endl;
	file <<"1"<<std::endl; // number-of-real-tags
	file <<"0.0"<<std::endl;
	file <<"3"<<std::endl; // number-of-int-tags
	file <<"0"<<std::endl;
	file <<"1"<<std::endl;
	file <<nbnode<<std::endl;
	for (int j=0;j<nbnode;j++){
		file << j+1 <<" " <<f[j]<<std::endl; // number type_elt (triangle) nbr_tags
	}
	file << "$EndNodeData\n";

	file.close();
}

template <class m_t, class space, class f_t> void write_gmsh_2(const m_t& Omega, const space & dof, const f_t & f, std::string filename, std::string name="\"\""){
	int nbdof = nb_dof(dof);
  int nbelt = nb_elt(dof);

	std::vector<R3> node(nbdof);
	for (int i =0 ; i<nbelt;i++){
		for (int j =0;j<dof.dim_loc;j++){
			node[dof[i][j]]=Omega[i][j];
		}
	}



	std::ofstream file; file.open((filename+".msh").c_str());
	file << "$MeshFormat\n";
	file << "2.2 0 8\n";
	file << "$EndMeshFormat\n";

	//== Nodes
	file << "$Nodes\n";
	file << nbdof << std::endl;
	for(int j=0; j<nbdof; j++){
		file << j+1 <<" "<< node[j][0]<<" "<<node[j][1]<<" "<<node[j][2] << "\n";
	}
	file << "$EndNodes\n";

	//== Elements
	file << "$Elements\n";
	file << nbelt << std::endl;
	for (int j=0;j<nbelt;j++){
		if (get_dim(Omega[j])==1){
			file << j+1 <<" " <<1 <<" "<<2<<" "<<0<<" "<<0<<" "; // number type_elt (ligne) nbr_tags
			file << dof[j][0] +1<<" ";
			file << dof[j][1] +1<<" "<<std::endl;
		}
		if (get_dim(Omega[j])==2){
			file << j+1 <<" " <<2 <<" "<<2<<" "<<0<<" "<<0<<" "; // number type_elt (triangle) nbr_tags
			file << dof[j][0] +1<<" ";
			file << dof[j][1] +1<<" ";
			file << dof[j][2] +1<<" "<<std::endl;
		}
		if (get_dim(Omega[j])==3){
			file << j+1 <<" " <<4 <<" "<<2<<" "<<0<<" "<<0<<" "; // number type_elt (tetrahedron) nbr_tags
			file << dof[j][0] +1<<" ";
			file << dof[j][1] +1<<" ";
			file << dof[j][2] +1<<" ";
			file << dof[j][3] +1<<" "<<std::endl;
		}



	}
	file << "$EndElements\n";

	//== Nodedata
	file << "$NodeData\n";
	file << 1<<std::endl; // number-of-string-tags
	file <<name<<std::endl;
	file <<"1"<<std::endl; // number-of-real-tags
	file <<"0.0"<<std::endl;
	file <<"3"<<std::endl; // number-of-int-tags
	file <<"0"<<std::endl;
	file <<"1"<<std::endl;
	file <<nbdof<<std::endl;
	for (int j=0;j<nbdof;j++){
		file << j+1 <<" " <<f[j]<<std::endl; // number type_elt (triangle) nbr_tags
	}
	file << "$EndNodeData\n";

	file.close();
}

}
#endif
