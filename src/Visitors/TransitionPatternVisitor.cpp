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

#include "TransitionPatternVisitor.h"
#include "../Quadrant.h"

namespace Clobscode
{

    TransitionPatternVisitor::TransitionPatternVisitor() {
        points = NULL;
        edges = NULL;
        max_ref_level = NULL;
    }

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
        
        vector<unsigned int> &pointindex = o->pointindex;
        EdgeVisitor ev;
        
        vector<unsigned int> nodes (8,0);
        bool splitted = false;
        //update the 4 nodes of the Quadrant.
        for (unsigned int i=0; i<4; i++) {
            nodes[i] = pointindex[i];
        }
        //search for nodes inserted in edges
        for (unsigned int i=0; i<4; i++) {
            QuadEdge e;
            ev.getEdge(o,i,e);
            set<QuadEdge>::iterator my_edge = edges->find(e);
            if (my_edge==edges->end()) {
                cout << "  edge " << e << " not found at applyTransitionPattern\n";
            }
            else {
                if ((*my_edge)[2]!=0) {
                    nodes[i+4] = (*my_edge)[2];
                    splitted = true;
                }
            }
        }
        //if this elements do not present nodes inserted in its edges
        //then return true (meaning this case is already considered in
        //the transition patterns) and add this element to the vector
        //of "new elements"
        if (!splitted) {
            return true;
        }
        
        //The middle node of the Quadrant can never be inserted
        //otherwise this Quadrant was already removed from the list
        //and replaced with 4 new Quadrants.
        
        //------------------------------------------------------
        //Finally, apply the transition pattern
        //------------------------------------------------------
        
        //creat the pattern
        patterns::QuadTransition qt (nodes);
        
        //the subelements of this Quadrant will no longer be a Quad.
        //It will now contain mixed-elements.
        vector<vector<unsigned int>> &sub_elements = o->sub_elements;
        sub_elements.clear();
        
        return qt.getNewElements(*points,sub_elements);
    }
}