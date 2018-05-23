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
* @file PolyEdge.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#ifndef PolyEdge_h
#define PolyEdge_h 1

#include "Point3D.h"
#include <vector>
#include <iostream>

using Clobscode::Point3D;
using std::vector;
using std::cout;
using std::ostream;

namespace PolyMesh
{
    class PolyEdge
	{
	public:		
        // Construction / destruction
        PolyEdge() {}

        PolyEdge(unsigned int point1, unsigned int point2){
            nodes.reserve(2);
            if(point1 < point2){
                nodes.push_back(point1);
                nodes.push_back(point2);
            }
            else{
                nodes.push_back(point2);
                nodes.push_back(point1);
            }
        }
        PolyEdge(vector<unsigned int> &fpts) : PolyEdge(fpts[0],fpts[1]) {}

        virtual ~PolyEdge(){

        }
        //--------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------

        //returns first node
        virtual unsigned int getKey() const;

        //returns second node
        virtual unsigned int getVal() const;

        virtual void computeNormal(vector<Point3D> &pts);

        //FJA useful?? besides, should be get?
        // as computeNormal set the normal, here, calculates a normal and returns
        //        virtual Point3D computeNormalAtNode(unsigned int nidx,vector<Point3D> &pts) const;

        //FJA usefull? why 2 points as parameters
        //        virtual Point3D projection(const Point3D &p1, const Point3D &p2) const;
		
        virtual vector<unsigned int> &getPoints();
        virtual const vector<unsigned int> &getPoints() const;

        virtual const Point3D &getNormal() const;
		
        virtual Point3D getNormalizedNormal() const;
		
        virtual Point3D getSemiNormalizedNormal() const;
		
		virtual unsigned int &operator[](unsigned int pos);
		
		virtual unsigned int operator[](unsigned int pos) const;
		
        //FJA transform triangle/edge intersection into what?
        //        virtual bool segmentIntersection(vector<Point3D> &pts,
        //                                         const Point3D &ep1,
        //                                         const Point3D &ep2);

        //computes the distant of pPoint to the line generated
        //by this segment.
        virtual double distance(const Point3D &pPoint) const;

        friend ostream& operator<<(ostream& o,const PolyEdge &p);

        friend bool operator==(const PolyEdge &p1,const PolyEdge &p2);

        friend bool operator!=(const PolyEdge &p1,const PolyEdge &p2);

	protected:
		// Data
        vector<unsigned int> nodes;
        Point3D mEdgeNormal, onepoint; //FJA onepoint will be used only in distance
		
	};
	
    inline vector<unsigned int> &PolyEdge::getPoints(){
        return nodes;
    }
    inline const vector<unsigned int> &PolyEdge::getPoints() const{
        return nodes;
    }

    inline const Point3D &PolyEdge::getNormal() const {
        return mEdgeNormal;
	}	
	
    inline Point3D PolyEdge::getNormalizedNormal() const {
        Point3D normal = mEdgeNormal/(mEdgeNormal.Norm());
		return normal;
	}
	
    inline Point3D PolyEdge::getSemiNormalizedNormal() const {
        std::cerr << "WARNING FJA: returns Normal divided by 1000.!\n";
        Point3D normal = mEdgeNormal/10000.;
		return normal;
    }

    inline unsigned int PolyEdge::getKey() const{
        return nodes[0];
    }

    inline unsigned int PolyEdge::getVal() const{
        return nodes[1];
    }

}

#endif
