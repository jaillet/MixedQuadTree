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


namespace Clobscode
{
    void EdgeVisitor::insertEdges(Quadrant *q, map<QuadEdge, EdgeInfo> &MapEdges) {
        
        vector<unsigned int> qpts = q->getPointIndex();
        unsigned int i = q->getIndex();
        
        auto found = MapEdges.find(QuadEdge (qpts[0],qpts[1]));
        if (found!=MapEdges.end()) {
            (found->second)[1] = i;
        }
        else {
            MapEdges.emplace(QuadEdge (qpts[0],qpts[1]), EdgeInfo (1,i));
        }
        
        found = MapEdges.find(QuadEdge (qpts[1],qpts[2]));
        if (found!=MapEdges.end()) {
            (found->second)[1] = i;
        }
        else {
            MapEdges.emplace(QuadEdge (qpts[1],qpts[2]), EdgeInfo (1,i));
        }
        
        found = MapEdges.find(QuadEdge (qpts[2],qpts[3]));
        if (found!=MapEdges.end()) {
            (found->second)[2] = i;
        }
        else {
            MapEdges.emplace(QuadEdge (qpts[2],qpts[3]), EdgeInfo (2,i));
        }
        
        found = MapEdges.find(QuadEdge (qpts[3],qpts[0]));
        if (found!=MapEdges.end()) {
            (found->second)[2] = i;
        }
        else {
            MapEdges.emplace(QuadEdge (qpts[3],qpts[0]), EdgeInfo (2,i));
        }
        
    }
    
    void EdgeVisitor::removeEdges(Quadrant *q, map<QuadEdge, EdgeInfo> &MapEdges) {
        
        vector<unsigned int> qpts = q->getPointIndex();
        
        auto e1 = MapEdges.find(QuadEdge (qpts[0], qpts[1]));
        e1->second[1] = std::numeric_limits<unsigned int>::max();
        auto e2 = MapEdges.find(QuadEdge (qpts[1], qpts[2]));
        e2->second[1] = std::numeric_limits<unsigned int>::max();
        auto e3 = MapEdges.find(QuadEdge (qpts[2], qpts[3]));
        e3->second[2] = std::numeric_limits<unsigned int>::max();
        auto e4 = MapEdges.find(QuadEdge (qpts[3], qpts[0]));
        e4->second[1] = std::numeric_limits<unsigned int>::max();
    }
    
    void EdgeVisitor::getEdge(const Quadrant *q, unsigned int idx, QuadEdge &e) {
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
                exit (EXIT_FAILURE);
        }
        // REM: possible name confusion with std::assign
        e.assign(e0,e1); // set QuadEdge indices of extremities
    }
}