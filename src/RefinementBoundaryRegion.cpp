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
* @file RefinementBoundaryRegion.cpp
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "RefinementBoundaryRegion.h"

namespace Clobscode
{	
	//--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
    RefinementBoundaryRegion::RefinementBoundaryRegion(const Polyline &input, const unsigned short &level)
        :RefinementRegion(level,true,true),ply(input)
	{
//		refine_level = level;
//        local_rot = false;
//        input_rot = false;
        
        //Use force input rotation (Defined in parent class) if the input domain
        //should be re-aligned.
	}

	//--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
    RefinementBoundaryRegion::~RefinementBoundaryRegion()
	{
		
	}

    
    //--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
    GeometricTransform RefinementBoundaryRegion::rotateWithinYou(Polyline &input) {
        //rotate input and return transformation
        GeometricTransform gt;
        gt.rotatePolyline(input);
        return gt;
        //gt defined in parent class
    }
    
	//--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
    void RefinementBoundaryRegion::rotate(GeometricTransform &gt) {
    
        //do nothing
    }


    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    unsigned int RefinementBoundaryRegion::nbEdgeIntersection(const Point3D &oct_p1, const Point3D &oct_p2,
                                            const Point3D &seg_p1, const Point3D &seg_p2) const
    {
        //compute the parametric equations of the given segment
        double dx = seg_p2[0]-seg_p1[0];
        double dy = seg_p2[1]-seg_p1[1];
        double dz = seg_p2[2]-seg_p1[2];
        double X,Y,Z,t;
        unsigned int nbInter=0;

        //Each 't' value is computed regarding a particular
        //plane. For instance, X_min (pmin[0]) is YZ plane
        //defined at X_min. When the 't' value is between 0
        //and 1, it means that the point where the segment
        //cuts the plane is between the two points that define
        //the segment. In any other case, the intersection is
        //produced at some other point outside the limits of
        //this segment (p1,p2), therefore the clipping can be
        //directly discard.

        ////////////////////////////////////////////////////
        // Important note: if tangent segment should be   //
        // considered as intersecting the Quadrant, replace //
        // each < by a <= from this point to the end of   //
        // this method.                                   //
        ////////////////////////////////////////////////////

        //test X axis over X_min
        t = (oct_p1[0]-seg_p1[0])/dx;

        if (0<=t && t<=1) {
            //test Y axis for this 't'
            Y = seg_p1[1] + t*dy;
            if (oct_p1[1]<=Y && Y<=oct_p2[1]) {
                ++nbInter;
//                //test Z axis for this 't'
//                Z = seg_p1[2] + t*dz;
//                if (oct_p1[2]<=Z && Z<=oct_p2[2]) {
//                    ++nbInter;
//                }
            }
        }

        //test X axis over X_max
        t = (oct_p2[0]-seg_p1[0])/dx;

        if (0<=t && t<=1) {
            //test Y axis for this 't'
            Y = seg_p1[1] + t*dy;
            if (oct_p1[1]<=Y && Y<=oct_p2[1]) {
                ++nbInter;
//                //test Z axis for this 't'
//                Z = seg_p1[2] + t*dz;
//                if (oct_p1[2]<=Z && Z<=oct_p2[2]) {
//                    ++nbInter;
//                }
            }
        }

        //test Y axis over Y_min
        t = (oct_p1[1]-seg_p1[1])/dy;

        if (0<=t && t<=1) {
            //test Y axis for this 't'
            X = seg_p1[0] + t*dx;
            if (oct_p1[0]<=X && X<=oct_p2[0]) {
                ++nbInter;
//                //test Z axis for this 't'
//                Z = seg_p1[2] + t*dz;
//                if (oct_p1[2]<=Z && Z<=oct_p2[2]) {
//                    ++nbInter;
//                }
            }
        }

        //test Y axis over Y_max
        t = (oct_p2[1]-seg_p1[1])/dy;
        if (0<=t && t<=1) {
            //test Y axis for this 't'
            X = seg_p1[0] + t*dx;
            if (oct_p1[0]<=X && X<=oct_p2[0]) {
                    ++nbInter;
//                //test Z axis for this 't'
//                Z = seg_p1[2] + t*dz;
//                if (oct_p1[2]<=Z && Z<=oct_p2[2]) {
//                    ++nbInter;
//                }
            }
        }

        //        //test Z axis over Z_min
        //        t = (oct_p1[2]-seg_p1[2])/dz;
        //        if (0<=t && t<=1) {
        //            //test Y axis for this 't'
        //            X = seg_p1[0] + t*dx;
        //            if (oct_p1[0]<=X && X<=oct_p2[0]) {
        //                //test Z axis for this 't'
        //                Y = seg_p1[1] + t*dy;
        //                if (oct_p1[1]<=Y && Y<=oct_p2[1]) {
        //                    ++nbInter;
        //                }
        //            }
        //        }

        //        //test Z axis over Z_max
        //        t = (oct_p2[2]-seg_p1[2])/dz;
        //        if (0<=t && t<=1) {
        //            //test Y axis for this 't'
        //            X = seg_p1[0] + t*dx;
        //            if (oct_p1[0]<=X && X<=oct_p2[0]) {
        //                //test Z axis for this 't'
        //                Y = seg_p1[1] + t*dy;
        //                if (oct_p1[1]<=Y && Y<=oct_p2[1]) {
        //                    ++nbInter;
        //                }
        //            }
        //        }

        //If the above test didn't succeed, the segment
        //is not intersected by this Quadrant
        return nbInter;
    }

	//--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
    bool RefinementBoundaryRegion::intersectsQuadrant(const vector<MeshPoint> &points, Quadrant &q) const
    {
        if (q.getRefinementLevel() >= getRefinementLevel() )
            return false;

        // test1: if nb Features >= 2
        // first condition is optimization if NbFeatures has been already computed
        if (q.getIntersectedFeatures().size()>=2 || ply.getNbFeatures(q,points)>=2)
            return true;

        // test2: if distFProjectionFeature > distMax for each node
        // done: compute distMax before...
        for (const uint idFeat:q.getIntersectedFeatures()) { // quad has feature?
            double minFeatDist=std::numeric_limits<double>::max();
            const Point3D &ptFeat=ply.getPoints()[idFeat];
            for (auto id:q.getPointIndex()) {
                Point3D pt = points[id].getPoint();
//                if (ply.pointIsInMesh(pt,q.getIntersectedEdges()))
                    minFeatDist=min(minFeatDist, (pt - ptFeat).Norm() );
                //TODO: should break here if minFeatDist < q.getMaxDistance()
            }
            if ( minFeatDist*1. > q.getMaxDistance() )
                return true;
        }

        // test3: if the number of intersection of the Polyline and Quadrant edge >= 3
        if (q.getIntersectedEdges().size()>=2) {
            unsigned int nbInter=0;
            //(Point3D &p1, Point3D &p2), diagonal extremities
            const Point3D &p1 = points[q.getPointIndex()[0]].getPoint();
            const Point3D &p2 = points[q.getPointIndex()[2]].getPoint();

            for (unsigned int idEdge:q.getIntersectedEdges()) {
                const PolyEdge &edge= ply.getEdges().operator [](idEdge);
                const Point3D &p3 = ply.getPoints()[edge[0]];
                const Point3D &p4 = ply.getPoints()[edge[1]];
                nbInter += nbEdgeIntersection(p1,p2,p3,p4);

            }
//            if (nbInter>=3)
//                return true;
            if (nbInter>4 || (nbInter>2 && q.hasIntersectedFeatures()) )
                return true;
        }

        // test4: if nb IntersectedEdges >= 3
        //TODO: remove this condition when others have been implemented...
//        if (q.getIntersectedEdges().size()>=3)
//            return true;

        return false;
	}
}
