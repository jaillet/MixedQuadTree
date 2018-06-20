/*
 <Mix-mesher: region type. This program generates a mixed-elements 2D mesh>

 Copyright (C) <2013,2018>  <Claudio Lobos> All rights reserved.

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU  Lesser General Public License as published by
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
* @file Point3D.cpp
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "Point3D.h"

namespace Clobscode
{

    string Point3D::print() const {
        ostringstream o;
        string out;

        o << x << " ";
        o << y << " ";
        o << z;

        out = o.str();
        return out;
    }

    std::ostream& operator<<(std::ostream& o, const Point3D &p){
        o << p.print().c_str();
        return o;
    }

    void Point3D::translateRotate(const Point3D &t, double xangle,
                                  double yangle, double zangle) {
        x-=t[0];
        y-=t[1];
        z-=t[2];

        if (xangle!=0) {
            xAxisRotation(xangle);
        }
        if (yangle!=0) {
            yAxisRotation(yangle);
        }
        if (zangle!=0) {
            zAxisRotation(zangle);
        }

    }

    void Point3D::rotateTranslate(const Point3D &t, double xangle,
                                  double yangle, double zangle) {

        if (zangle!=0) {
            zAxisRotation(-zangle);
        }
        if (yangle!=0) {
            yAxisRotation(-yangle);
        }
        if (xangle!=0) {
            xAxisRotation(-xangle);
        }

        x+=t[0];
        y+=t[1];
        z+=t[2];
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    // Compute Angle at P0P1-P1P2 in radians
    // P0 previous, this=P1 mid, P2 next point
    double Point3D::angle3Points(const Point3D &P0, const Point3D &P2) const {

        Point3D V1 = P0-(*this), V2 = P2-(*this);
        V1.normalize();
        V2.normalize();
        return (acos(V1.dot(V2)));

    }


}
