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
* @file RefinementDrawingRegion.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#ifndef RefinementDrawingRegion_h
#define RefinementDrawingRegion_h 1

#include <vector>
#include <iostream>
#include "RefinementRegion.h"
#include "Visitors/IntersectionsVisitor.h"

using Clobscode::Point3D;
using std::vector;
using std::cout;
using Clobscode::RefinementRegion;
using Clobscode::IntersectionsDrawingVisitor;
using Clobscode::Quadrant;

namespace Clobscode
{
    class RefinementDrawingRegion : public RefinementRegion
	{
	public:
				
	RefinementDrawingRegion(const Polyline&, const unsigned short &level);
				
	virtual ~RefinementDrawingRegion();
        
        virtual GeometricTransform rotateWithinYou(Polyline &input) override;
		
        virtual void rotate(GeometricTransform &gt) override;
		
        virtual bool intersectsQuadrant(const vector<MeshPoint> &points, Quadrant &q) const override;

	protected:
	    Polyline ply;
	};
}

#endif
