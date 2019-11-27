/*
 <Mix-mesher: region type. This program generates a mixed-elements 2D mesh>

 Copyright (C) <2013,2019>  <Claudio Lobos> All rights reserved.

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
* @file RefinementDrawingRegion.cpp
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "RefinementDrawingRegion.h"
#define _USE_MATH_DEFINES
#include <math.h>

namespace Clobscode
{
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    RefinementDrawingRegion::RefinementDrawingRegion(const Polyline &aPly, const unsigned short &aLevel)
        :RefinementRegion(aLevel,false,false),ply(aPly)
    {

    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    RefinementDrawingRegion::~RefinementDrawingRegion()
    {

    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    GeometricTransform RefinementDrawingRegion::rotateWithinYou(Polyline &input) {
        //rotate input and return transformation
        GeometricTransform gt;
        gt.rotatePolyline(input);
        return gt;
        //gt defined in parent class
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void RefinementDrawingRegion::rotate(GeometricTransform &gt){

        //do nothing

    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    bool RefinementDrawingRegion::intersectsQuadrant(const vector<MeshPoint> &points, Quadrant &q) const
    {

        IntersectionsDrawingVisitor idv;
        idv.setPolyline(ply);
        idv.setPoints(points);

        return (idv.visit(&q));
    }
}
