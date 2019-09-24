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
* @file Quadrant.cpp
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "Quadrant.h"
#include "GeometricTransform.h"

namespace Clobscode
{
	
	
	//--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
	Quadrant::Quadrant(vector<unsigned int> &epts, 
                   const unsigned short &ref_level)
        :pointindex(epts),ref_level(ref_level),
          surface(false),max_dis(numeric_limits<double>::infinity()) {
        
        /***** BEGIN Debugging variables *******/
              debugging = false;
        /***** END Debugging variables *******/
		sub_elements.assign(1,pointindex);
	}

    Quadrant::Quadrant(vector<unsigned int> &epts, const Quadrant &quad) : pointindex(epts),
                                                                           ref_level(quad.ref_level), surface(quad.surface), debugging(quad.debugging),
                                                                           max_dis(quad.max_dis) {
        sub_elements.assign(1,pointindex);
        intersected_edges = quad.intersected_edges;
        intersected_features = quad.intersected_features;
    }
	
	//--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
	Quadrant::~Quadrant(){
		
	}

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    bool Quadrant::badAngle(unsigned int &nIdx, const vector<MeshPoint> &mp) const {

        if (pointindex.size()!=4) {
            cerr << "Warning in Quadrant::badAngle";
            cerr << " Splitted quadrant not managed yet\n";
            return false;
        }

        const Point3D &P0 = mp[pointindex[(nIdx+3)%4]].getPoint(); //previous point
        const Point3D &P1 = mp[pointindex[nIdx]].getPoint();       //mid point
        const Point3D &P2 = mp[pointindex[(nIdx+1)%4]].getPoint(); //next point
        
        double angle = P1.angle3Points(P0,P2);

        //cerr << "Angle " << angle << "\n";
        
        if (isnan(angle)) {
            return true;
        }
        
        //if the angle is smoth (close to 180) or if it is going to produce inverted
        //elements in a quad (>180), then it must be managed with transition patterns
        //if (angle>2.61799 && angle<3.66519) {
        if (angle>150.0) { // && angle<210.) {
            return true;
        }
        return false;
        
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    double Quadrant::getAngle(unsigned int &nIdx, const vector<MeshPoint> &mp) const {
        
        if (pointindex.size()!=4) {
            cerr << "Warning in Quadrant::badAngle";
            cerr << " Splitted quadant not managed yet\n";
            return false;
        }
        
//        unsigned int i0 = (4+nIdx-1)%4, i2=(nIdx+1)%4;
        
        const Point3D &P0 = mp[pointindex[(nIdx+3)%4]].getPoint(); //previous point
        const Point3D &P1 = mp[pointindex[nIdx]].getPoint();       //mid point
        const Point3D &P2 = mp[pointindex[(nIdx+1)%4]].getPoint(); //next point        
        
        //if two nodes are the same, we return NaN
        if ((P0-P1).Norm()<0.001) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if ((P2-P1).Norm()<0.001) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        
        //double angle = toDegrees( P1.angle3Points(P0,P2));
        //function angle returns degrees in 360.
        double angle = P1.angle3Points(P0,P2);
        
        if (isnan(angle)) {
            return 180;
        }
        
        return angle;
        
    }
    
	//--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    bool Quadrant::pointInside(const vector<MeshPoint> &mp, const Point3D &p) const {
        const Point3D &minp = mp[pointindex[0]].getPoint();
        const Point3D &maxp = mp[pointindex[2]].getPoint();
        if (minp[0]>p[0] || maxp[0]<p[0]) {
            return false;
        }
        if (minp[1]>p[1] || maxp[1]<p[1]) {
            return false;
        }
        if (minp[2]>p[2] || maxp[2]<p[2]) {
            return false;
        }
        
        //cout << "Quadrant with feature\n" << minp << " " << maxp << " -> " << p << "\n";
        return true;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    bool Quadrant::accept(Visitor *v)
    {
        return v->visit(this);
    }

	//--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    std::ostream& operator<<(std::ostream& o, const Quadrant &q){
        std::vector<unsigned int> points = q.getPointIndex();
		for (unsigned int i=0; i<points.size(); i++)
			o << " " << points[i];
		return o;
	}

}
