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
* @file PolyEdge.cpp
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "PolyEdge.h"

namespace PolyMesh
{

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    // pts contains all the vertices of the polyline
    void PolyEdge::computeNormal(vector<Point3D> &pts)
    {
        // edge normal computation
        // vector orthogonal to A and B
        mEdgeNormal[0] = -(pts[nodes[1]][1]-pts[nodes[0]][1]);
        mEdgeNormal[1] =   pts[nodes[1]][0]-pts[nodes[0]][0] ;

        //FJA will be used only in distance
        onepoint = pts[nodes[0]];
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    double PolyEdge::distance(const Point3D &p) const{
        Point3D vd = p - onepoint;
        return vd*mEdgeNormal;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

//    Point3D SurfTriangle::projection(const Point3D &p1, const Point3D &p2) const{


//        double d1 = distance(p1), d2 = distance(p2), s;
//        s = d1/(d1-d2);

//        Point3D dif = p2 - p1;
//        dif *= s;

//        return p1 - dif;
//    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    unsigned int &PolyEdge::operator[](unsigned int pos){
        return nodes.at(pos);
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    unsigned int PolyEdge::operator[](unsigned int pos) const{
        return nodes.at(pos);
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

//    bool SurfTriangle::segmentIntersection(vector<Point3D> &pts,
//                                           const Point3D &ep1,
//                                           const Point3D &ep2){
//        Point3D p0 = pts[mIdxV[0]];
//        Point3D p1 = pts[mIdxV[1]];
//        Point3D p2 = pts[mIdxV[2]];

//        //compute the distance to each point of the edge
//        double dis_min = (ep1 - p0).dot(mTriangleNormal);
//        double dis_max = (ep2 - p0).dot(mTriangleNormal);

//        //if both, ep1 and ep2 are at the same side of the
//        //plane, there is no intersection.
//        if ((dis_max>0 && dis_min>0) || (dis_max<0 && dis_min<0)) {
//            return false;
//        }

//        /*if (dis_min==0) {

//        }*/

//        //Compute the intersection between the segment
//        //and triangle's plane.
//        double s = dis_min/(dis_min - dis_max);

//        Point3D point_in_plane = ep1 + (ep2-ep1)*s;

//        Point3D u = p1 - p0;
//        Point3D v = p2 - p0;
//        Point3D w = point_in_plane - p0;

//        double uu=u.dot(u);
//        double vv=v.dot(v);
//        double uv=u.dot(v);
//        double wv=w.dot(v);
//        double wu=w.dot(u);
//        double d =uv*uv - uu*vv;

//        if(d==0){
//            cout << "error !!! malformed triangle\n";
//            cout << "at HexTriIntersection::";
//            cout << "edgeTriangleIntersection";
//            return false;
//        }

//        double si= (uv*wv - vv*wu)/d;
//        double ti= (uv*wu - uu*wv)/d;
//        if(si>=0 && ti>=0 && (si+ti)<=1){
//            return true;
//        }

//        return false;
//    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    Point3D PolyEdge::getNormalAtNode(unsigned int nidx,
                                      vector<Point3D> &pts ) const
    {
        //return mEdgeNormal;
        unsigned int lidx = 2; // trick to avoid true/false variable

        if (nidx == nodes[0]) {
            lidx = 0;
        }
        else {
            if (nidx == nodes[1]) {
                lidx = 1;
            }
        }

        if (lidx==2) { //should be 0 or 1 if found
            Point3D v;
            std::cerr << "warning at PolyMesh::PolyEdge::getNormalAtNode";
            std::cerr << " node index not found in edge\n";
            return v;
        }

        // edge lenght computation
        double  length = (pts[nodes[(lidx+1)%2]] - pts[nodes[lidx]]).Norm();
        return length * mEdgeNormal;
    }


    ostream& operator<<(ostream& o, const PolyEdge &e){
        o << "e(" << e.getKey() << ",";
        o << e.getVal() << ")";
        return o;
    }

    bool operator==(const PolyEdge &e1,const PolyEdge &e2){
        if(e1.getKey()==e2.getKey() && e1.getVal()==e2.getVal())
            return true;
        return false;
    }

    bool operator!=(const PolyEdge &p1,const PolyEdge &p2){
        return !(p1==p2);
    }

}
