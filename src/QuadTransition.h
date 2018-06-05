/*
 <Mix-mesher: region type. This program generates a mixed-elements mesh>
 
 Copyright (C) <2013,2017>  <Claudio Lobos>
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/gpl.txt>
 */

#ifndef QuadTransition_h
#define QuadTransition_h 1

#include <iostream>
#include <vector>
#include <list>
#include "MeshPoint.h"

using std::list;
using std::vector;
using std::ostream;
using Clobscode::MeshPoint;

namespace patterns {

	class QuadTransition{
		
	public:
		
		QuadTransition(vector<unsigned int> &nodes);
		
        virtual bool getNewElements(vector<vector<unsigned int> > &sub_elements,
                                    const unsigned int &nedges);
        
    protected:
        
        virtual bool oneEdge(vector<vector<unsigned int> > &sub_elements);
        
        virtual bool twoEdges(vector<vector<unsigned int> > &sub_elements);
        
        virtual bool threeEdges(vector<vector<unsigned int> > &sub_elements);
        
        virtual bool fourEdges(vector<vector<unsigned int> > &sub_elements);
        
        virtual vector<unsigned int> rotate(const vector<unsigned int> &nodes,
                                             const unsigned int &times);
        
        virtual vector<unsigned int> inverseRotate(const vector<unsigned int> &nodes,
                                                    const unsigned int &times);
		
	protected:
		
		vector<unsigned int> nodes;
		
	};
	
}
#endif
