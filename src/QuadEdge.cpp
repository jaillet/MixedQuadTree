/*
 <Mix-mesher: region type. This program generates a mixed-elements 2D mesh>

 Copyright (C) <2013,2018>  <Claudio Lobos> All rights reserved.

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/lgpl.txt>
 */
/**
* @file QuadEdge.cpp
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "QuadEdge.h"

namespace Clobscode
{
	
// vector<unsigned int> info; //info[2] midpoint

	QuadEdge::QuadEdge(){
		info.assign(2,0);
	}
	
    QuadEdge::QuadEdge(unsigned int point1, unsigned int point2){
        info.resize(2);
        assign(point1, point2);
    }

	QuadEdge::~QuadEdge(){
		
	}

    // REM: possible name confusion with std::assign
//    void QuadEdge::assign(const unsigned int &point1, const unsigned int &point2){
    void QuadEdge::assign(unsigned int point1, unsigned int point2){
        if (point1<point2) {
			info[0]=point1;
			info[1]=point2;
		}
		else {
			info[1]=point1;
			info[0]=point2;
		}
	}
	
	/*bool QuadEdge::split(set<QuadEdge> &allQuadEdges, unsigned int maxp){
		
		pair<set<QuadEdge>::iterator , bool> result;
		//create possible new QuadEdge
		QuadEdge ne1(info[0],maxp);
		
		cout << "trying to insert " << ne1 << "\n";
		
		result = allQuadEdges.insert(ne1);
		
		
		if (!result.second) {
			cout << "cannot split this QuadEdge " << *this << "\n";
			return false;
		}
		else {
			
			QuadEdge ne2(info[1],maxp);
			allQuadEdges.insert(ne2);
			cout << "QuadEdge " << *this << " was split\n";
		}
		return true;

	}*/
	
	ostream& operator<<(ostream& o, const QuadEdge &e){
		o << e[0] << " ";
        o << e[1];
		return o;
	}
	
    bool operator==(const QuadEdge &e1, const QuadEdge &e2) {
		//this possible as QuadEdges are sorted by min index
		if (e1[0]==e2[0] && e1[1]==e2[1])
			return true;
		return false;
	}
	
    bool operator!=(const QuadEdge &e1, const QuadEdge &e2) {
		return !(e1==e2);
	}

    bool operator<(const QuadEdge &e1, const QuadEdge &e2) {
		if (e1[0]!=e2[0]){
			return e1[0]<e2[0];
		}
		return e1[1] < e2[1];
	}
	
    //FJA using default operator=
//    QuadEdge& QuadEdge::operator=(const QuadEdge &e){
//        if (this != &e) { // self-assignment check expected
//            info[0]=e[0];
//            info[1]=e[1];
//            info[2]=e[2];
//        }
//        return *this;
//    }

/* FJA just here for recall on op=, to check bug in the one right above ;o(
// assume the object holds reusable storage, such as a heap-allocated buffer mArray
T& operator=(const T& other) // copy assignment
{
    if (this != &other) { // self-assignment check expected
        if (other.size != size) {         // storage cannot be reused
            delete[] mArray;              // destroy storage in this
            size = 0;
            mArray = nullptr;             // preserve invariants in case next line throws
            mArray = new int[other.size]; // create storage in this
            size = other.size;
        } 
        std::copy(other.mArray, other.mArray + other.size, mArray);
    }
    return *this;
}
T& operator=(T&& other) noexcept // move assignment
{
    if(this != &other) { // no-op on self-move-assignment (delete[]/size=0 also ok)
        delete[] mArray;                               // delete this storage
        mArray = std::exchange(other.mArray, nullptr); // leave moved-from in valid state
        size = std::exchange(other.size, 0);
    }
    return *this;
}*/

}
