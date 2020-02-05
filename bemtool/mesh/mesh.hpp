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
#ifndef BEMTOOL_MESH_MESH_HPP
#define BEMTOOL_MESH_MESH_HPP

#include <vector>
#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include "element.hpp"
#include "../calculus/calculus.hpp"


namespace bemtool{
    
    //==========================//
    //   Numerotation locale    //
    //   d'un element           //
    //==========================//
    
    class EltData{
        
    private:
        typedef std::map<const void*,int>       map_t;
        typedef map_t::const_iterator      it_t;
        
        std::map<const void*,int> num;
        
    public:
        EltData(const EltData& ed): num(ed.num) {};
        
        template <typename m_t>
        EltData(const m_t& m, const int& j) {num[&m]=j;}
        
        template <typename m_t>
        void Push(const m_t& m, const int& j){num[&m]=j;}
        
        template <typename m_t>
        int operator[](const m_t& m) const {
            it_t it = num.find(&m);
            if(it==num.end()){return -1;}
            return it->second;}
        
        friend const int Size(const EltData& data){return data.num.size();}
        
        friend std::ostream& operator<<(std::ostream& os, const EltData& ed){
            for(it_t it = ed.num.begin(); it != ed.num.end(); it++){
                os << "mesh:\t" << it->first  << "\t";
                os << "num:\t"  << it->second << std::endl; }
            return os;}
        
    };
    
    
    //==========================//
    //    Registre d'elements   //
    //==========================//
    
    template <typename EltType>
    class EltList{
        
    private:
        typedef EltList<EltType>   this_t;
        typedef EltType            e_t;
        typedef EltData            ed_t;
        static const int           end = -1;
        
        std::vector<e_t >   elt;
        std::vector<ed_t>   data;
        std::vector<int >   first;
        std::vector<int >   next;
        
    public:
        EltList(){};
        
        void Init(const int& nb_node){first.resize(nb_node,-1);}
        
        const e_t&                operator[](const int&    j) const {return elt[j];}
        friend const int          NbElt     (const this_t& l)       {return l.elt.size();}
        friend const std::vector<e_t>& EltOf     (const this_t& l)       {return l.elt;}
        
        const ed_t& operator[](const e_t& e) const{
            const int& begin = first[Key(e)];
            for(int p=begin; p!=end; p=next[p]){
                if(e==elt[p]){return data[p];} }
            std::cout << "error: missing element" << std::endl;
            exit(EXIT_FAILURE);
            return data[0];
        }
        
        template <typename m_t> int Push(e_t e, const m_t& m, const int& j){
            int& begin = first[Key(e)];
            for(int p=begin; p!=end; p=next[p]){
                if(e==elt[p]){data[p].Push(m,j); return p;}}
            
            elt.push_back(e);
            next.push_back(begin);
            begin = elt.size()-1;
            data.push_back(EltData(m,j));
            return begin;
        }
        
        
    };
    
    
    typedef EltList<Elt0D> ListElt0D;
    typedef EltList<Elt1D> ListElt1D;
    typedef EltList<Elt2D> ListElt2D;
    typedef EltList<Elt3D> ListElt3D;
    
    
    //==============//
    //    Geometry  //
    //==============//
    
    class Geometry{
        
    public:
        typedef R3 v_t;
        
    private:
        std::string      mesh_file;
        std::vector<R3>  node;
        ListElt1D   elt1D;
        ListElt2D   elt2D;
        ListElt3D   elt3D;
        
        inline void LoadNodes();
        
    public:
        
        template <typename f_t> Geometry(const f_t& f): mesh_file(f) {(*this).LoadNodes();}
        Geometry(const std::vector<v_t>&);
        Geometry(){}; // add default constructor AF 07/19
        
        inline void initEltData();
        inline void setnodes(const R3& p);
        
        const R3&       operator[](const int&   j) const {return node[j]; }
        const EltData&  operator[](const Elt1D& e) const {return elt1D[e];}
        const EltData&  operator[](const Elt2D& e) const {return elt2D[e];}
        const EltData&  operator[](const Elt3D& e) const {return elt3D[e];}
        
        template <class i_t> const subarray<const Geometry,i_t> operator[](const i_t& i_) const {
            return subarray<const Geometry,i_t>(*this,i_);}
        
        template <typename m_t> int Push(Elt1D e, const m_t& m, const int& j){return elt1D.Push(e,m,j);}
        template <typename m_t> int Push(Elt2D e, const m_t& m, const int& j){return elt2D.Push(e,m,j);}
        template <typename m_t> int Push(Elt3D e, const m_t& m, const int& j){return elt3D.Push(e,m,j);}
        
        inline friend const int             NbNode    (const Geometry& g){return g.node.size();  }
        inline friend const int             NbElt1D   (const Geometry& g){return NbElt(g.elt1D); }
        inline friend const int             NbElt2D   (const Geometry& g){return NbElt(g.elt2D); }
        inline friend const int             NbElt3D   (const Geometry& g){return NbElt(g.elt3D); }
        inline friend const std::vector<Elt1D>&  Elt1DOf   (const Geometry& g){return EltOf(g.elt1D); }
        inline friend const std::vector<Elt2D>&  Elt2DOf   (const Geometry& g){return EltOf(g.elt2D); }
        inline friend const std::vector<Elt3D>&  Elt3DOf   (const Geometry& g){return EltOf(g.elt3D); }
        inline friend const std::string&         MeshFileOf(const Geometry& g){return g.mesh_file;    }
        
    };
    
    
    template<int D> struct GetElt{static inline const std::vector< Elt<D> >& Of(const Geometry& g);};
    template<> inline const std::vector<Elt1D>& GetElt<1>::Of(const Geometry& g){return Elt1DOf(g);}
    template<> inline const std::vector<Elt2D>& GetElt<2>::Of(const Geometry& g){return Elt2DOf(g);}
    template<> inline const std::vector<Elt3D>& GetElt<3>::Of(const Geometry& g){return Elt3DOf(g);}
    
    
    Geometry::Geometry(const std::vector<v_t>& x){
        for(int j=0; j<x.size(); j++){node.push_back(x[j]);}
        initEltData();
    }
    
    inline void Geometry::LoadNodes(){
        
        assert(!node.size());
        std::ifstream file(mesh_file.c_str());
        
        if( file.fail() ){
            std::cout << "error: reading nodes" << std::endl;
            exit(EXIT_FAILURE);}
        
        std::string line;
        std::istringstream iss;
        int poubelle;
        R3 p;
        
        while( line != "$Nodes" ){
            getline(file,line);}
        
        // Chargt nb de noeuds
        file >> poubelle;
        getline(file,line);
        
        // Chargt donnees noeuds
        getline(file,line);
        while( line != "$EndNodes" ){
            iss.str(line);
            iss >> poubelle;
            iss>>p;
            node.push_back(p);
            iss.clear();
            getline(file,line);
        }
        
        file.close();
        
        initEltData();
        
    }
    
    // factorization for the initialization of the element vectors   AF 07/2019
    inline void Geometry::initEltData(){
        elt1D.Init(node.size());
        elt2D.Init(node.size());
        elt3D.Init(node.size());
    }
    // nodes setter AF 07/2019
    inline void Geometry::setnodes(const R3& p){node.push_back(p);}
    
    
    
    
    //==========================//
    //     Mesh Class           //
    //==========================//
    
    enum Boundedness{bounded, unbounded};
    
    template <int D> class Mesh{
        
    public:
        typedef Mesh<D>   this_t;
        typedef Elt<D>       e_t;
        static const int dim = D;
        
    private:
        Geometry*               node;
        const std::vector<e_t>*      elt;
        std::vector<int>             num_elt;
        std::vector<R3>              normal;
        Boundedness             boundedness;
        
    public:
        Mesh<D>(): boundedness(bounded){};
        Mesh<D>(Geometry& g): node(&g), boundedness(bounded), elt(&GetElt<D>::Of(g)) {};
        Mesh<D>(const Mesh<D>& m): node(m.node), boundedness(bounded) {assert(!NbElt(m));}
        
        friend const Geometry&    GeometryOf(const this_t& m)       {return *(m.node);}
        friend const std::vector<e_t>& EltOf     (const this_t& m)       {return *(m.elt);}
        friend int                NbElt     (const this_t& m)       {return m.num_elt.size();}
        const e_t&                operator[](const int& j   ) const {return (*elt)[num_elt[j]];}
        bool IsBounded() const {return boundedness==bounded;}
        
        void operator=(const Boundedness& b){
            if(boundedness!=b){
                boundedness=b;
                for(int j=0; j<normal.size(); j++){normal[j] = (-1.)*normal[j];}
            }
        }
        
        void operator+=(const this_t& m){
            node = m.node;
            elt  = m.elt;
            for(int j=0; j<NbElt(m); j++){
                int J = node->Push( m[j], *this, num_elt.size() );
                num_elt.push_back(J);} }
        
        template <class r_t> void operator<<(const r_t& r_){
            int J = node->Push( e_t(r_,&(*node)[0]), *this, num_elt.size() );
            num_elt.push_back(J); }
        
        friend void WriteMedit(const this_t& m, char const * const name){
            Geometry& node = *(m.node);
            int nb_node = NbNode(node);
            std::ofstream file; file.open(name);
            file << "MeshVersionFormatted 1\nDimension 3\nVertices\n";
            file << nb_node << std::endl;
            for(int j=0; j<nb_node; j++){file<<node[j]<<"\t 0 \n";}
            file << "\nTriangles\n" << NbElt(m) << std::endl;
            for(int j=0; j<NbElt(m); j++){ file << Num(m[j])+1 << "\t0 \n";}
            file.close();
        }
        
        friend const std::vector<R3>& NormalTo(const this_t& m){return m.normal;}
        template <typename v_t> void Orienting(const v_t& v){
            normal.resize(num_elt.size());
            for(int j=0; j<normal.size(); j++){normal[j]=v[j];}}
        
        // void Load(Geometry& g){node = &g;elt  = &GetElt<dim>::Of(g);}
        
        void Load(Geometry& g, int ref=-1){ // -1 : everything is loaded
            
            node = &g;
            elt  = &GetElt<dim>::Of(g);
            int nb_loaded_elt = 0;
            
            //  Geometry& node = *(m.node);
            array<dim+1,int> I;
            std::string filename = MeshFileOf(g);
            
            // Variables  auxiliaires
            int poubelle, elt_type;
            int tag, nb_tags;
            
            // Traitement chaines caractere
            std::string line;
            std::istringstream iss;
            std::ifstream file(filename.c_str(), std::ifstream::in);
            if(file.fail()){std::cout << "error: loading mesh\n"; exit(EXIT_FAILURE);}
            
            // Deplacement rubrique elements
            while(line != "$Elements"){getline(file,line);}
            getline(file,line);
            getline(file,line);
            
            // Lecture rubrique elements
            while( line != "$EndElements" ){
                iss.str(line);
                iss >> poubelle;
                iss >> elt_type;
                iss >> nb_tags;
                iss >> tag;
                
                if((elt_type==dim || (dim==3 && elt_type==4)) && (tag == ref || ref==-1)  ){
                    for(int j=0; j<nb_tags-1; j++){iss>>poubelle;}
                    iss >> I; *this << g[I-1];
                    nb_loaded_elt++;
                }
                
                iss.clear();
                getline(file,line);
            }
            
            file.close();
            
            if(nb_loaded_elt==0){
                std::cout << "error: no element loaded" << std::endl;
                exit(EXIT_FAILURE);}
            
        }
        
        void setNodesGeomety(Geometry& myg) {
            this->node=myg;
        }
        
        // geometry and element vectors setter AF 07/2019
        inline void set_elt(Geometry& g){node = &g; elt=&GetElt<D>::Of(g);}
        // setter for an element AF 07/2019
        inline void setOneElt(Geometry& g,array<dim+1,int> I){ *this << g[I]; }
       
        
    };
    
    typedef Mesh<1> Mesh1D;
    typedef Mesh<2> Mesh2D;
    typedef Mesh<3> Mesh3D;
    
} //namespace

#endif
