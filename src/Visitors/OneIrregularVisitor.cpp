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


    OneIrregularVisitor::OneIrregularVisitor() :mapedges(NULL), max_ref_level(NULL)
    { }

    OneIrregularVisitor::OneIrregularVisitor(const map<QuadEdge, unsigned int> *mapedges,
                                             const unsigned short *max_ref_level)
        :mapedges(mapedges), max_ref_level(max_ref_level)
    { }

    
    void OneIrregularVisitor::setMapEdges(const map<QuadEdge, unsigned int> &mapedges) {
        this->mapedges = &mapedges;
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
            
            auto search = mapedges->find(ee);
            if (search==mapedges->end()) {
                cerr << "Edge not found at OneIrregularVisitor::visit\n";
            }

            unsigned int mid_idx = mapedges->at(ee);

            //if the edge is not split, check the others
            if (mid_idx==0) {
                continue;
            }
            //At this point, the edge is split so both
            //"sub-edges" must be checked. If at least one of
            //them is also split, then this element is not
            //one-irregular

            QuadEdge sub1(ee[0],mid_idx);
            
            search = mapedges->find(sub1);
            if (search==mapedges->end()) {
                cerr << "Edge not found at OneIrregularVisitor::visit\n";
            }
            
            if (mapedges->at(sub1)!=0) {
                return false;
            }
            QuadEdge sub2(ee[1],mid_idx);
            
            search = mapedges->find(sub2);
            if (search==mapedges->end()) {
                cerr << "Edge not found at OneIrregularVisitor::visit\n";
            }
            
            if (mapedges->at(sub2)!=0) {
                return false;
            }
        }
        return true;
    }
}
