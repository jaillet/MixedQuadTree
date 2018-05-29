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
* @file RefinementRegion.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "RefinementRegion.h"

namespace Clobscode
{
//unsigned short refine_level;
//bool local_rot, input_rot;
//vector<Point3D> pts;
//GeometricTransform gt;

	//--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
    RefinementRegion::RefinementRegion(unsigned short level, bool local_rot, bool input_rot )
        :refine_level(level),local_rot(local_rot),input_rot(local_rot)
	{	}
	//--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
	RefinementRegion::~RefinementRegion()
	{
		
	}
	
}
