#ifndef MISC_H
#define MISC_H

//===================================================================
//
//  Copyright 2015 Xavier Claeys
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


#include <ctime>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

namespace bemtool{



////========================================================////
////////////===== Barre de progression ======///////////////////

class progress{
    
private:
    const char* title;
    const int length;
    int       it;
    int       prg;
    clock_t   t0;
    int       verbose;
    
public:
    progress(const char* aff,const int& l, int v=1): title(aff), length(l) {
        t0 = clock();
        it=0; prg=0;
        verbose=v;
        
        if (verbose>0){
            std::cout << "\r";
            std::cout << title << ": \t";
            std::cout << 0 << "%";
            std::cout.flush();
        }
    }
    
    void operator++(int n){
        it++;

        if( int(it*100./length)>prg & verbose>0){
            prg=int(it*100./length);
            std::cout << "\r";
            std::cout << title << ": \t";
            std::cout << prg << "%";
            std::cout.flush();
            
        }
    }
    
    void end(){
        t0 = clock()-t0;
        time_t now; time(&now);
        if (verbose>0){
            std::cout << "\r";
            std::cout << title << ": \t";
            std::cout << ((float)t0)/CLOCKS_PER_SEC << " sec." << std::endl;
        }
    }
    
};


////========================================================////
/////////////////===== Conversions ======///////////////////////

template <typename nbr>
std::string NbrToStr(nbr N){
	std::ostringstream strs;
	strs << N;
	std::string str = strs.str();
	return str;
}

template <typename T>
T StrToNbr ( const std::string &Text )
{
	std::istringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
}


////========================================================////
///////////////////////===== Input ======///////////////////////
inline std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
inline std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

}


#endif
