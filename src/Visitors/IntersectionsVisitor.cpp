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
* @file IntersectionsVisitor.cpp
* @author Created by nanairo on 08-03-16.
* @authors Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "IntersectionsVisitor.h"

#include <stdexcept>
#include "../Quadrant.h"

namespace Clobscode

//Polyline *mesh;
//vector<MeshPoint> *points;
//list<unsigned int> *edges;
//vector<Point3D> *coords;
//bool select_edges;

{
    IntersectionsVisitor::IntersectionsVisitor()
        :ply(NULL),points(NULL),edges(NULL),coords(NULL),select_edges(false)
    {    }


    IntersectionsVisitor::IntersectionsVisitor(bool select_edges)
        :ply(NULL),points(NULL),edges(NULL),coords(NULL),select_edges(select_edges)
    {    }

    void IntersectionsVisitor::setPolyline( Polyline &ply) {
        this->ply = &ply;
    }
    void IntersectionsVisitor::setPoints(vector<MeshPoint> &points){
        this->points = &points;
    }
    void IntersectionsVisitor::setEdges(list<unsigned int> &edges){
        this->edges = &edges;
    }
    void IntersectionsVisitor::setCoords(vector<Point3D> &coords){
        this->coords = &coords;
    }


    bool IntersectionsVisitor::visit(Quadrant *q) {
        //cout << "IntersectionsVisitor" << endl;
        //check intersections with selected input edges
        if (select_edges) //checkIntersections(PolyEdge &ply,list<unsigned int> &edges,vector<Point3D> &coords);
        {
            if (edges == NULL || coords == NULL)
                throw std::runtime_error(std::string("Calling wrong IntersectionsVisitor!"));
            list<unsigned int> &intersected_edges = q->intersected_edges;

            list<unsigned int>::const_iterator e_iter;
            const vector<Point3D> &input_pts = ply->getPoints();

            for (e_iter = edges->begin(); e_iter!=edges->end(); e_iter++) {
                const PolyEdge &edge = ply->getEdges()[*e_iter];
                if (intersectsEdge(edge,input_pts,coords->at(0),coords->at(1))) {
                    intersected_edges.push_back(*e_iter);
                }
            }

            return !intersected_edges.empty();
        }
        //check intersections with all input edges
        else //checkIntersections(PolyEdge &ply, vector<MeshPoint> &pts)
        {
            if (points == NULL)
                throw std::runtime_error(std::string("Calling wrong IntersectionsVisitor!"));
            const vector<unsigned int> &pointindex = q->pointindex;
            list<unsigned int> &intersected_edges = q->intersected_edges;

            const vector<PolyEdge> &edges = ply->getEdges();
            const vector<Point3D> &input_pts = ply->getPoints();

            for (unsigned int j=0; j<edges.size(); j++) {
                if (intersectsEdge(edges[j],input_pts,
                                       points->at(pointindex[0]).getPoint(),  //bbox vertices
                                       points->at(pointindex[2]).getPoint()))
                {
                    intersected_edges.push_back(j);
                }
            }
            return !intersected_edges.empty();
        }

    }

    //auxiliary functions
    bool IntersectionsVisitor::intersectsEdge(const PolyMesh::PolyEdge &pEdge, const vector<Point3D> &input_pts, const Point3D &pmin,
                                             const Point3D &pmax) const {

        //pmin and pmax describe the bounding box of the
        //Quadrant. They are used to detect intersections
        //between input edges and the Quadrant with a clipping
        //method.

        //This code is based on the Cohen-Sutherland algorithm
        //for image clipping over a given window. In this case
        //the "window" is the Quadrant (square), and the "line"
        //to be clipped is the edge.

        const vector<unsigned int> &iEdge = pEdge.getPoints(); //indices of edge extremities

        unsigned int n_pts = iEdge.size(); //FJA always 2 ?????

        vector<unsigned int> sides(n_pts,0);

        for (unsigned int i=0; i<n_pts; i++) {
            sides[i] = computePosition(input_pts[iEdge[i]],
                                       pmin,pmax);
            if (sides[i] == 0) {
                return true;
            }
        }

        //Compute logic and between bits associated to each
        //triangle node. If all of them are at same side (e.g.
        //all at left), the result will be different from zero,
        //meaning that there is no intersection.
        unsigned int all_and = sides[0] & sides[1];
        //If all are at the same side, return false.
        if (all_and != 0) {
            return false;
        }

        //Patch: all other cases, we cannot be sure.
        //return true;

        //TEST EDGES//

        //Find out if the edge cross the quad. Two cases: both are
        //"centered", but one is on the left and the other on
        //the right (same for Top/bottom)
        // FJA@discussion (and Front/Back).
        //Tested with logic or. The other case is that their
        //are in exact opposite "2D corners". Tested with XoR.
        // FJA@discussion1 I think this is not true for 00 11 11
        // if side[p1]=9 up-left, and p[2]=6 bottom-right for example
        // and p1 high enough and p2 right enough, there is no intersection
        // are we supposing they are close enougth from the quadrant???
        // I guess no, so this test (15 <> 00 11 11) is wrong !
        // and the 3D one as well (63 <> 11 11 11) is wrong !
        // FJA@discussion2 and exact opposite "3D corners"
        unsigned int result_or = sides[0]|sides[1];
//        unsigned int result_xor = sides[0]^sides[1];
        if (result_or == 48) {
            //48 <> 11 00 00 z
            std::cerr << "Warning in IntersectionsVisitor::intersectsEdge()\n";
            std::cerr << "  case 48 should never happen, as z=0.0\n";
            return true;
        }
        if (result_or == 12) {
            //12 <> 00 11 00 y
            return true;
        }
        if (result_or == 3) {
            //3 <> 00 00 11 x
            return true;
        }
//        if (result_xor == 63) {
//            //63 <> 11 11 11
//            return true;
//        }

        //General case of clipping this edge against
        //the Quadrant.
        const Point3D &p1 = input_pts[iEdge[0]];
        const Point3D &p2 = input_pts[iEdge[1]];
        /* FJA: for me, if we are here it's because clipping
        * returned false for all sides of the quadrant
        * I think we should trust it, and return false,
        * the edge is not intersecting this quadrant */
//        if (clipGeneralCase(p1,p2,pmin,pmax)) {
//            return true;
//        }
        return (clipGeneralCase(p1,p2,pmin,pmax));

        //Last "hard" test: Apply Ray Tracing to detect
        //intersection between every edge of the Quadrant
        //(square) and the polyline. If at least one
        //of them intersects, then this triangle
        //intersects the Quadrant.
        //If the test didn't succeed, the triangle
        //is not intersected by this Quadrant

        //return edgeTriangleIntersection(pEdge,input_pts,pmin,pmax);
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    bool IntersectionsVisitor::clipGeneralCase(const Point3D &p1, const Point3D &p2, const Point3D &pmin,
                                          const Point3D &pmax) const {

        //compute the parametric equations of the given segment
        double dx = p2[0]-p1[0];
        double dy = p2[1]-p1[1];
        double dz = p2[2]-p1[2];
        double X,Y,Z,t;

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
        t = (pmin[0]-p1[0])/dx;

        //if (0<t && t<1) {
        if (0<=t && t<=1) {
            //test Y axis for this 't'
            Y = p1[1] + t*dy;
            //if (pmin[1]<Y && Y<pmax[1]) {
            if (pmin[1]<=Y && Y<=pmax[1]) {
                //test Z axis for this 't'
                Z = p1[2] + t*dz;
                //if (pmin[2]<Z && Z<pmax[2]) {
                if (pmin[2]<=Z && Z<=pmax[2]) {
                    return true;
                }
            }
        }

        //test X axis over X_max
        t = (pmax[0]-p1[0])/dx;

        //if (0<t && t<1) {
        if (0<=t && t<=1) {
            //test Y axis for this 't'
            Y = p1[1] + t*dy;
            //if (pmin[1]<Y && Y<pmax[1]) {
            if (pmin[1]<=Y && Y<=pmax[1]) {
                //test Z axis for this 't'
                Z = p1[2] + t*dz;
                //if (pmin[2]<Z && Z<pmax[2]) {
                if (pmin[2]<=Z && Z<=pmax[2]) {
                    return true;
                }
            }
        }


        //test Y axis over Y_min
        t = (pmin[1]-p1[1])/dy;

        //if (0<t && t<1) {
        if (0<=t && t<=1) {
            //test Y axis for this 't'
            X = p1[0] + t*dx;
            //if (pmin[0]<X && X<pmax[0]) {
            if (pmin[0]<=X && X<=pmax[0]) {
                //test Z axis for this 't'
                Z = p1[2] + t*dz;
                //if (pmin[2]<Z && Z<pmax[2]) {
                if (pmin[2]<=Z && Z<=pmax[2]) {
                    return true;
                }
            }
        }

        //test Y axis over Y_max
        t = (pmax[1]-p1[1])/dy;

        //if (0<t && t<1) {
        if (0<=t && t<=1) {
            //test Y axis for this 't'
            X = p1[0] + t*dx;
            //if (pmin[0]<X && X<pmax[0]) {
            if (pmin[0]<=X && X<=pmax[0]) {
                //test Z axis for this 't'
                Z = p1[2] + t*dz;
                //if (pmin[2]<Z && Z<pmax[2]) {
                if (pmin[2]<=Z && Z<=pmax[2]) {
                    return true;
                }
            }
        }

        std::cerr << "Warning IntersectionsVisitor::clipGeneralCase()\n";
        std::cerr << "  FJA: usefull to test z in 2D???? as dz=0 and t=NaN\n";

        //test Z axis over Z_min
        t = (pmin[2]-p1[2])/dz;

        //if (0<t && t<1) {
        if (0<=t && t<=1) {
            //test Y axis for this 't'
            X = p1[0] + t*dx;
            //if (pmin[0]<X && X<pmax[0]) {
            if (pmin[0]<=X && X<=pmax[0]) {
                //test Z axis for this 't'
                Y = p1[1] + t*dy;
                //if (pmin[1]<Y && Y<pmax[1]) {
                if (pmin[1]<=Y && Y<=pmax[1]) {
                    return true;
                }
            }
        }

        //test Z axis over Z_max
        t = (pmax[2]-p1[2])/dz;

        //if (0<t && t<1) {
        if (0<=t && t<=1) {
            //test Y axis for this 't'
            X = p1[0] + t*dx;
            //if (pmin[0]<X && X<pmax[0]) {
            if (pmin[0]<=X && X<=pmax[0]) {
                //test Z axis for this 't'
                Y = p1[1] + t*dy;
                //if (pmin[1]<Y && Y<pmax[1]) {
                if (pmin[1]<=Y && Y<=pmax[1]) {
                    return true;
                }
            }
        }

        return false;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    bool IntersectionsVisitor::edgeTriangleIntersection(const PolyMesh::PolyEdge &st,
                                                   const vector<Point3D> &input_pts,
                                                   const Point3D &pmin,
                                                   const Point3D &pmax) const {
        std::cerr << "Warning IntersectionsVisitor::edgeTriangleIntersection()\n";
        std::cerr << "  FJA: I didn't checked this function in 2D yet! And you, Claudio?\n";

        //compute the coords of all Quadrant edges
        const vector<vector<Point3D> > &oct_edges = getEdges(pmin,pmax);

        //test each edge against the triangle
        for (unsigned int i=0; i<oct_edges.size(); i++) {

            if (st.segmentIntersection(input_pts,oct_edges[i][0],oct_edges[i][1])) {
                return true;
            }
        }

        return false;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    vector<vector<Point3D>> IntersectionsVisitor::getEdges(const Point3D &pmin,
                                                      const Point3D &pmax) const {
        vector<vector<Point3D> > edges;
        edges.reserve(4);

        //create the 4 nodes of the Quadrant
        Point3D p0 (pmin[0],pmin[1],0.0); //FJA assume pmin[2]=pmax[2]=0.0;
        Point3D p1 (pmax[0],pmin[1],0.0);
        Point3D p2 (pmax[0],pmax[1],0.0);
        Point3D p3 (pmin[0],pmax[1],0.0);

        vector<Point3D> edge(2, Point3D ());
        //CCW
        edge[0] = p0;
        edge[1] = p1;
        edges.push_back(edge);

        edge[0] = p1;
        edge[1] = p2;
        edges.push_back(edge);

        edge[0] = p2;
        edge[1] = p3;
        edges.push_back(edge);

        edge[0] = p3;
        edge[1] = p0;
        edges.push_back(edge);

        return edges;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    unsigned int IntersectionsVisitor::computePosition(const Point3D &p, const Point3D &pmin, const Point3D &pmax) const {

        unsigned int sides = 0;

        //decide left/right X axis.
        if (p[0]<pmin[0]) {
            sides += 1;
        }
        else if (pmax[0] < p[0]) {
            sides += 2;
        }

        //decide bottom/top Y axis.
        if (p[1]<pmin[1]) {
            sides += 4;
        }
        else if (pmax[1] < p[1]) {
            sides += 8;
        }

        //decide back/front Z axis.
        if (p[2]<pmin[2]) {
            sides += 16;
            std::cerr << "FJA: Warning in IntersectionsVisitor::computePosition()\n";
            std::cerr << "  case 16 should never happen, as z=0.0\n";
        }
        else if (pmax[2]<p[2]) {
            sides += 32;
            std::cerr << "FJA: Warning in IntersectionsVisitor::computePosition()\n";
            std::cerr << "  case 32 should never happen, as z=0.0\n";
        }

        return sides;
    }
}


