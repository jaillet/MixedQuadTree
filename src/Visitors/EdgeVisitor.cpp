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
* @file EdgeVisitor.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "EdgeVisitor.h"
#include "../Quadrant.h"


namespace Clobscode
{
void EdgeVisitor::insertEdges(Quadrant *q, set<QuadEdge> &edges) {

    QuadEdge ee;
    // assume 4 edges in Quadrant
    for (unsigned int i=0; i<4; i++) {
        getEdge(q,i,ee);
        edges.insert(ee);
    }
    //FJA for debugging purpose
    // => OK, pb with qtcreator/gdb not updating after std::set reallocation
//    std::cerr << "edges contains:";
//    for (auto it=edges.begin(); it!=edges.end(); ++it)
//        std::cerr << " - " << *it << "-" << &(*it);
//    std::cerr << '\n'<< std::flush;

    }

    void EdgeVisitor::getEdge(Quadrant *q, unsigned int idx, QuadEdge &e) {
        const vector<unsigned int> &pointindex = q->pointindex;
        unsigned int e0,e1;
        switch (idx) {
            case 0:
                e0 = pointindex[0];
                e1 = pointindex[1];
                break;
            case 1:
                e0 = pointindex[1];
                e1 = pointindex[2];
                break;
            case 2:
                e0 = pointindex[2];
                e1 = pointindex[3];
                break;
            case 3:
                e0 = pointindex[3];
                e1 = pointindex[0];
                break;
            default:
                std::cerr << "In EdgeVisitor::getEdge(Quadrant *q, unsigned int idx, QuadEdge &e),"
                          << "bad idx value=" << idx << std::endl;
                break;
        }
        // REM: possible name confusion with std::assign
        e.assign(e0,e1); // set QuadEdge indices of extremities
    }
}
