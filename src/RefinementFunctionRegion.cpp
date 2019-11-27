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
* @file RefinementAllRegion.cpp
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "RefinementFunctionRegion.h"

namespace Clobscode
{
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    RefinementFunctionRegion::RefinementFunctionRegion(const unsigned short &level)
        :RefinementRegion(level,false,false)
    {

    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    RefinementFunctionRegion::~RefinementFunctionRegion()
    {

    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    GeometricTransform RefinementFunctionRegion::rotateWithinYou(Polyline &input) {
        //rotate input and return transformation
        GeometricTransform gt;
        gt.rotatePolyline(input);
        return gt;
        //gt defined in parent class
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void RefinementFunctionRegion::rotate(GeometricTransform &gt){

        //do nothing

    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    bool RefinementFunctionRegion::intersectsQuadrant(const vector<MeshPoint> &points, Quadrant &q) const
    {
        uint nbIn=0;
        //Ellipse parameters
        Point3D center(0.5,0.75,0.0);
        double xradius=.25,yradius=.5;

        // Check whether some points are inside the ellipse
        for (auto quaNoIdx:q.getPointIndex()){
            if (pow((points[quaNoIdx].getPoint()[0]-center[0])/xradius,2) +
                    pow((points[quaNoIdx].getPoint()[1]-center[1])/yradius,2 )<= 1)
                ++nbIn;
        }

        // Some points (but not all) are in = Quad intersects the ellipse
        if ( nbIn>0 && nbIn<4 )
            return true;
        else if (nbIn==0) { // all points are outside
            // test if ellipse center is inside Quad
            const Point3D& p1 = points[q.getPointIndex()[0]].getPoint();
            const Point3D& p2 = points[q.getPointIndex()[2]].getPoint();
            if ( p1[0]<=center[0] && center[0]<=p2[0] &&
                 p1[1]<=center[1] && center[1]<=p2[1]) {
                return true;
            }
        }
        return false;
    }
}
