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
 * @file OneIrregularVisitor.cpp
 * @author Claudio Lobos, Fabrice Jaillet
 * @version 0.1
 * @brief
 **/

#include "OneIrregularVisitor.h"

#include "../Quadrant.h"

namespace Clobscode
{
    
    
    //Note: this class is now used for debugging proposes only. For this reason
    //it will print a lot of information of non-balanced Octants.
    //The actual balancing is performed in the Mesher class.
    
    
    OneIrregularVisitor::OneIrregularVisitor() :edges(NULL)
    { }
    
    OneIrregularVisitor::OneIrregularVisitor(const map<QuadEdge, EdgeInfo> *MapEdges)
    :edges(MapEdges)
    { }
    
    void OneIrregularVisitor::setEdges(const map<QuadEdge, EdgeInfo> &MapEdges) {
        this->edges = &MapEdges;
    }
    
    bool OneIrregularVisitor::visit(Quadrant *o) {
        
        unsigned int oidx = o->getIndex();
        
        //EdgeVisitor ev;
        for (unsigned int i=0; i<4; i++) {
            
            QuadEdge ee;
            EdgeVisitor::getEdge(o,i,ee);
            auto my_edge = edges->find(ee);
            
            if (my_edge==edges->end()) {
                cerr << "Error at Quadrant::isOneIrregular!!!\n";
                cerr << "Edge not found\n";
                std::abort();
            }
            
            unsigned int mid_idx = (my_edge->second)[0];
            
            //if the edge is not split, check the others
            if (mid_idx==0) {
                continue;
            }
            //At this point, the edge is split so both
            //"sub-edges" must be checked. If at least one of
            //them is also split, then this element is not
            //one-irregular
            
            auto sub_edge = edges->find(QuadEdge (ee[0],mid_idx));
            if ((sub_edge->second)[0]!=0) {
                
                cout << "Quad " << oidx << " not balanced\n";
                
                return false;
            }

            sub_edge = edges->find(QuadEdge (ee[1],mid_idx));
            if ((sub_edge->second)[0]!=0) {
                
                cout << "Quad " << oidx << " not balanced\n";
                
                return false;
            }
        }
        return true;
    }
}