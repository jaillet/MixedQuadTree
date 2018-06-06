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
* @file TransitionPatternVisitor.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "TransitionPatternVisitor.h"
#include "../Quadrant.h"

namespace Clobscode
{

    TransitionPatternVisitor::TransitionPatternVisitor()
        :points(NULL),edges(NULL),max_ref_level(NULL)
    {    }

    void TransitionPatternVisitor::setPoints(vector<MeshPoint> &points) {
        this->points = &points;
    }

    void TransitionPatternVisitor::setEdges(const set<QuadEdge> &edges) {
        this->edges = &edges;
    }

    void TransitionPatternVisitor::setMaxRefLevel(const unsigned short &max_ref_level) {
        this->max_ref_level = &max_ref_level;
    }

    bool TransitionPatternVisitor::visit(Quadrant *o) {

        //if this Quadrant is refined to the maximum level, return it immediately
        if (*max_ref_level == o->ref_level) {
            return true;
        }
        
        const vector<unsigned int> &pointindex = o->pointindex;
        
        //number of mid-edges.
        unsigned int quantity = 0;
        // 4 vertex + one node per Edge => 8
        vector<unsigned int> nodes (8,0);

        //update the 4 nodes of the Quadrant.
        for (unsigned int i=0; i<4; i++) {
            nodes[i] = pointindex[i];
        }
        //search for nodes inserted in edges
        for (unsigned int i=0; i<4; i++) {
            QuadEdge ee;
            EdgeVisitor::getEdge(o,i,ee);
            set<QuadEdge>::const_iterator my_edge = edges->find(ee);
            if (my_edge==edges->end()) {
                cout << "  edge " << ee << " not found at applyTransitionPattern\n";
            }
            else {
                if ((*my_edge)[2]!=0) {
                    nodes[i+4] = (*my_edge)[2];
                    quantity++;
                }
            }
        }
        //if this elements do not present nodes inserted in its edges
        //then return true (meaning this case is already considered in
        //the transition patterns) and add this element to the vector
        //of "new elements"
        if (quantity==0) {
            return true;
        }
        
        //The middle node of the Quadrant can never be inserted
        //otherwise this Quadrant was already removed from the list
        //and replaced with 4 new Quadrants.
        
        //------------------------------------------------------
        //Finally, apply the transition pattern
        //------------------------------------------------------
        
        //creat the pattern
        patterns::QuadTransition qt(nodes);
        
        //the subelements of this Quadrant will no longer be a Quad.
        //It will now contain mixed-elements.
        
        cout << "\n\n\n hello!!! \n";
        
        bool done = qt.getNewElements(o->sub_elements,quantity);
        
        cout << "new number of subele: " << (o->sub_elements).size() << endl;
        return done;
    }
}
