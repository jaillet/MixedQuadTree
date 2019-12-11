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

#ifndef OneIrregularVisitor_h
#define OneIrregularVisitor_h 1

#include "../MeshPoint.h"
#include "../QuadEdge.h"
#include "../Point3D.h"
#include "../EdgeInfo.h"

#include <map>
#include <vector>

#include "Visitor.h"
#include "EdgeVisitor.h"

namespace Clobscode
{
    class OneIrregularVisitor : public Visitor {
    public:
        OneIrregularVisitor();
        
        OneIrregularVisitor(const map<QuadEdge, EdgeInfo> *MapEdges);
        
        bool visit(Quadrant *o) override; // Quad is const, but no const for override...
        
        void setEdges(const map<QuadEdge, EdgeInfo> &MapEdges);
        
    protected:
        const map<QuadEdge, EdgeInfo> *edges;

    };
    
}


#endif //MESHER_ROI_ONEIRREGULARVISITOR_H