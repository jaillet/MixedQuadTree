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

#ifndef TransitionPatternVisitor_h
#define TransitionPatternVisitor_h 1

#include "../QuadEdge.h"
#include "../MeshPoint.h"
#include "../QuadTransition.h"

#include "Visitor.h"
#include "EdgeVisitor.h"

#include <list>
#include <set>
#include <vector>

using Clobscode::MeshPoint;
using Clobscode::QuadEdge;
using Clobscode::EdgeVisitor;
using std::list;
using std::set;
using std::vector;



namespace Clobscode
{
    class TransitionPatternVisitor : public Visitor{
    public:
        TransitionPatternVisitor();

        bool visit(Quadrant *o) override;

        void setPoints(vector<MeshPoint> &points);
        void setEdges(const set<QuadEdge> &edges);
        void setMaxRefLevel(const unsigned short &max_ref_level);

    private:
        vector<MeshPoint> *points;
        const set<QuadEdge> *edges;
        const unsigned short *max_ref_level;
    };
}


#endif //MESHER_ROI_TRANSITIONPATTERNVISITOR_H
