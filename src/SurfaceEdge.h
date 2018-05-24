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
* @file SurfaceEdge.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#ifndef SurfaceEdge_h
#define SurfaceEdge_h 1

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
	
    // Describe connectivity
    // For each point, vector of indices of edges
	class SurfaceEdge{
		
	public:
        
		SurfaceEdge();
		
		SurfaceEdge(const unsigned int &idx1, const unsigned int &idx2);
		
		virtual ~SurfaceEdge();
		
		virtual void assign(const unsigned int &idx1, const unsigned int &idx2);
		
		virtual unsigned int operator[](unsigned const int &pos) const;
		
		virtual void operator=(const SurfaceEdge &e);
		
		friend ostream& operator<<(ostream& o, const SurfaceEdge &e);
		
		friend bool operator==(const SurfaceEdge &e1, const SurfaceEdge &e2);
		
		friend bool operator!=(const SurfaceEdge &e1, const SurfaceEdge &e2);
		
		friend bool operator<(const SurfaceEdge &e1, const SurfaceEdge &e2);
		
		
	protected:
		
        //
		vector<unsigned int> info;

	};
		
	inline unsigned int SurfaceEdge::operator[](unsigned const int &pos) const{
		return info.at(pos);
	}
	
}
#endif
