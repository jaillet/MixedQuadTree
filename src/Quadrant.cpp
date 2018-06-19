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


namespace Clobscode
{
	
	
	//--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
	Quadrant::Quadrant(vector<unsigned int> &epts, 
                   const unsigned short &ref_level)
        :pointindex(epts),ref_level(ref_level),
          n_influences(0),influence_commit(false),surface(false),max_dis(numeric_limits<double>::infinity()),feature(false) {
        
		sub_elements.assign(1,pointindex);
	}
	
	//--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
	Quadrant::~Quadrant(){
		
	}
	
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    bool Quadrant::badAngle(const unsigned int &nIdx, vector<MeshPoint> &mp) const {
        
        
        unsigned int i0, i2;
        Point3D P1 = mp[nIdx].getPoint(), P0, P2;

        P0 = mp[(4+nIdx-1)%4].getPoint();
        P2 = mp[(nIdx+1)%4].getPoint();
        
        Point3D V1 = P0-P1, V2 = P2-P1;
        V1.normalize();
        V2.normalize();
        double angle = acos(V1.dot(V2));
        /*if (angle<2.61799 || angle>3.66519) {
            return true;
        }*/
        return false;
        
    }
    
	//--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    bool Quadrant::pointInside(vector<MeshPoint> &mp, const Point3D &p) const {
        Point3D minp = mp[pointindex[0]].getPoint();
        Point3D maxp = mp[pointindex[2]].getPoint();
        if (minp[0]>p[0] || maxp[0]<p[0]) {
            return false;
        }
        if (minp[1]>p[1] || maxp[1]<p[1]) {
            return false;
        }
        if (minp[2]>p[2] || maxp[2]<p[2]) {
            return false;
        }
        
        cout << "Quadrant with feature\n" << minp << " " << maxp << " -> " << p << "\n";
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
	std::ostream& operator<<(std::ostream& o,Quadrant &e){
        std::vector<unsigned int> points = e.getPointIndex();
		for (unsigned int i=0; i<points.size(); i++)
			o << " " << points[i];
		return o;
	}

}
