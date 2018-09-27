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
* @file RefinementAllRegion.cpp
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "RefinementAllRegion.h"

namespace Clobscode
{	
	//--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
	RefinementAllRegion::RefinementAllRegion(const unsigned short &level)
        :RefinementRegion(level)
	{
//		refine_level = level;
//        local_rot = false;
//        input_rot = false;
        
        //Use force input rotation (Defined in parent class) if the input domain
        //should be re-aligned.
	}

	//--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
	RefinementAllRegion::~RefinementAllRegion()
	{
		
	}

    
    //--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
    GeometricTransform RefinementAllRegion::rotateWithinYou(Polyline &input) {
        //rotate input and return transformation
        GeometricTransform gt;
        gt.rotatePolyline(input);
        return gt;
        //gt defined in parent class
    }
    
	//--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
    void RefinementAllRegion::rotate(GeometricTransform &gt) {
    
        //do nothing
    }
    	
	//--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
    bool RefinementAllRegion::intersectsQuadrant(const vector<MeshPoint> &points, Quadrant &oct) const
    {
        return true;
	}
}
