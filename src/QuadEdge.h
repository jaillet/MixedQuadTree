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
* @file QuadEdge.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#ifndef QuadEdge_h
#define QuadEdge_h 1

#include <iostream>
#include <vector>
#include <set>

using std::vector;
using std::ostream;
using std::set;
using std::cout;
using std::pair;

namespace Clobscode
{
	
	//class QuadEdge;
	
	class QuadEdge{
		
	public:
        
		QuadEdge();
		
        QuadEdge(unsigned int idx1, unsigned int idx2);

		virtual ~QuadEdge();
		
        virtual void assign(unsigned int idx1, unsigned int idx2);
		
//        virtual void updateMidPoint(unsigned int idx);
		
//		virtual bool split(set<QuadEdge> &allQuadEdges, unsigned int maxp);
		
//        virtual void setMidPoint(unsigned int mid);
		
        virtual unsigned int operator[](unsigned int pos) const;
		
//        virtual QuadEdge& operator=(const QuadEdge &e);

		friend ostream& operator<<(ostream& o, const QuadEdge &e);
		
        friend bool operator==(const QuadEdge &e1, const QuadEdge &e2);
		
        friend bool operator!=(const QuadEdge &e1, const QuadEdge &e2);
		
        friend bool operator<(const QuadEdge &e1, const QuadEdge &e2);
		
		
	protected:
		
        vector<unsigned int> info; //info[2] midpoint

	};
	
    /*inline void QuadEdge::updateMidPoint(unsigned int idx){
		info[2] = idx;
	}*/

    inline unsigned int QuadEdge::operator[](unsigned int pos) const{
//		return info.at(pos);
        return info[pos]; //no position validity check
    }
	
    //FJA TODO: this is a duplicate of updateMidPoint
    /*inline void QuadEdge::setMidPoint(unsigned int mid){
		info[2] = mid;
	}*/
	
}
#endif
