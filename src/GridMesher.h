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
* @file GridMesher.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#ifndef GridMesher_h
#define GridMesher_h 1

#include "MeshPoint.h"
#include <math.h>

using std::vector;
using std::list;
using Clobscode::MeshPoint;

namespace Clobscode
{
	class GridMesher{
		
	public:
		
		GridMesher();
		
		virtual ~GridMesher();
		
        virtual void generatePoints(vector<double> &bounds,
                                    vector<double> &all_x,
                                    vector<double> &all_y);
		
		virtual void generateMesh(vector<double> &all_x,
						vector<double> &all_y,
						vector<MeshPoint> &points,
						vector<vector<unsigned int> > &elements);
		
	protected:
		
		virtual void generateVector(vector<double> &coords, 
									double min, double max, 
									double step);
		
	};
	
}
#endif
