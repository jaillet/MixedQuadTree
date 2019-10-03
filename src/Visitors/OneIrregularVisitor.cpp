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

//set<QuadEdge> *edges;
//const unsigned short *max_ref_level;


    OneIrregularVisitor::OneIrregularVisitor() :MapEdges(NULL), max_ref_level(NULL)
    { }

    OneIrregularVisitor::OneIrregularVisitor(const map<QuadEdge, EdgeInfo> *MapEdges,
                                             const unsigned short *max_ref_level)
        :MapEdges(MapEdges), max_ref_level(max_ref_level)
    { }

    
    void OneIrregularVisitor::setMapEdges(const map<QuadEdge, EdgeInfo> &MapEdges) {
        this->MapEdges = &MapEdges;
    }

    void OneIrregularVisitor::setMaxRefLevel(const unsigned short &max_ref_level) {
        this->max_ref_level = &max_ref_level;
    }

    bool OneIrregularVisitor::visit(Quadrant *o) {

        if (*max_ref_level==o->ref_level) {
            return true;
        }

        const vector<unsigned int> &pi = o->pointindex;
        for (unsigned int i=0; i<4; i++) {

            QuadEdge ee(pi[i],pi[(i+1)%4]);
            
            auto search = MapEdges->find(ee);
            if (search==MapEdges->end()) {
                cerr << ee << " ";
                cerr << "Edge not found at OneIrregularVisitor::visit\n";
            }

            unsigned int mid_idx = (search->second)[0];

            //if the edge is not split, check the others
            if (mid_idx==0) {
                continue;
            }
            //At this point, the edge is split so both
            //"sub-edges" must be checked. If at least one of
            //them is also split, then this element is not
            //one-irregular
            
            search = MapEdges->find(QuadEdge (ee[0],mid_idx));
            if (search==MapEdges->end()) {
                cerr << "Edge not found at OneIrregularVisitor::visit\n";
            }
            
            if ((search->second)[0]!=0) {
                return false;
            }
            
            search = MapEdges->find(QuadEdge (ee[1],mid_idx));
            if (search==MapEdges->end()) {
                cerr << "Edge not found at OneIrregularVisitor::visit\n";
            }
            
            if ((search->second)[0]!=0) {
                return false;
            }
        }
        return true;
    }
}
