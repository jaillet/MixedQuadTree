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
* @file RemoveSubElementsVisitor.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#ifndef RemoveSubElementsVisitor_h
#define RemoveSubElementsVisitor_h 1

#include "../MeshPoint.h"
#include "Visitor.h"

#include <list>
#include <set>
#include <vector>

using Clobscode::MeshPoint;
using std::list;
using std::vector;


namespace Clobscode
{
    class RemoveSubElementsVisitor : public Visitor{
    public:
        RemoveSubElementsVisitor();

        bool visit(Quadrant *o) override;

        void setPoints(const vector<MeshPoint> &points);

    private:
        const vector<MeshPoint> *points;
    };
}



#endif //MESHER_ROI_REMOVESUBELEMENTSVISITOR_H
