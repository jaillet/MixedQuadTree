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
* @file Polyline.cpp
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "Polyline.h"

namespace Clobscode
{
	
	//--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    Polyline::Polyline() {
    
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    Polyline::Polyline(vector<Point3D> &pts,
                       vector<vector<unsigned int> > &edg_ind) : mVertices(pts) {

		if (pts.empty()) {
			std::cout << "Error in Mesh::Init input mesh without points\n";
			return;
		}
		
        // init bouding box with first Point3D
        computeBounds();
        checkNormalToPlane();

        /*for (unsigned int i=0; i<bounds.size(); ++i) {
            cout << bounds[i] << endl;
        }*/
        
        mEdges.reserve(edg_ind.size());
        // initialising the faces
        for (unsigned int iedg = 0; iedg < edg_ind.size(); ++iedg) {
            PolyEdge edg(edg_ind[iedg]);
            edg.computeNormal(mVertices);
            mEdges.push_back(edg);
        }

		// computing the pseudo normal at each surface node
        computeNodesPseudoNormal();
		
	}

    Polyline::~Polyline(){
		
	}
	
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void Polyline::computeBounds(){

        bounds.resize(6,0.0);

        // .at(0) throw exception if empty polyline
        bounds[0]=bounds[3]=mVertices.at(0)[0];
        bounds[1]=bounds[4]=mVertices.at(0)[1];
        bounds[2]=bounds[5]=mVertices.at(0)[2];

        // search for max and min coodinates in X, Y and Z.
        double x,y,z;
        for (unsigned int i=1; i< mVertices.size(); ++i) {

            x = mVertices[i][0];
            y = mVertices[i][1];
            z = mVertices[i][2];

            if(bounds[0]>x)
                bounds[0]=x;
            else if(bounds[3]<x)
                bounds[3]=x;
            if(bounds[1]>y)
                bounds[1]=y;
            else if(bounds[4]<y)
                bounds[4]=y;
            if(bounds[2]>z)
                bounds[2]=z;
            else if(bounds[5]<z)
                bounds[5]=z;
        }
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void Polyline::computeEdgesNormal(){
        for (auto &e:mEdges) {
            e.computeNormal(mVertices);
        }
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void Polyline::computeNodesPseudoNormal(){

        unsigned int npts = mVertices.size();
        mVerticePseudoNormals.assign(npts,Point3D());
        vector<list<unsigned int> > seg_per_node(npts);

        //save a reference to the segment for each node
        for (unsigned int i=0; i<mEdges.size(); i++) {
            seg_per_node[mEdges[i][0]].push_back(i);
            seg_per_node[mEdges[i][1]].push_back(i);
        }
		
		//compute normal of each node
        for (unsigned int i=0; i<npts; i++) {
            for (auto iseg:seg_per_node[i]) {
                mVerticePseudoNormals[i] += mEdges[iseg].getNormalAtNode(i,mVertices);
			}
			
            if (seg_per_node[i].empty()) {
				continue;
			}
			//normalize
			mVerticePseudoNormals[i].normalize();
		}
	}
	
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    // cn_PnPoly(): crossing number test for a point in a polygon
    //      Input:   P = a point,
    //               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
    //      Return:  0 = outside, 1 = inside
    // This code is patterned after http://geomalgorithms.com/a03-_inclusion.html
    int Polyline::crossingNumber(const Point3D &P) const
    {
        int    cn = 0;    // the  crossing number counter

        // loop through all edges of the polygon
        for (unsigned int i=0; i<mEdges.size(); i++) {    // edge from V[i]  to V[i+1]
            const Point3D &P0=mVertices[mEdges[i][0]];
            const Point3D &P1=mVertices[mEdges[i][1]];
            if (((P0.Y() <= P.Y()) && (P1.Y() > P.Y()))     // an upward crossing
            || ((P0.Y() > P.Y()) && (P1.Y() <=  P.Y()))) { // a downward crossing
                // compute  the actual edge-ray intersect x-coordinate
                double vt = (P.Y()  - P0.Y()) / (P1.Y() - P0.Y());
                if (P.X() <  P0.X() + vt * (P1.X() - P0.X())) // P.x < intersect
                     ++cn;   // a valid crossing of y=P.y right of P.x
            }
        }
        return (cn&1);    // 0 if even (out), and 1 if  odd (in)
    }

    // wn_PnPoly(): winding number test for a point in a polygon
    //      Input:   P = a point,
    //               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
    //      Return:  wn = the winding number (=0 only when P is outside)
    // This code is patterned after http://geomalgorithms.com/a03-_inclusion.html
    int Polyline::windingNumber(const Point3D &P) const
    {
        int    wn = 0;    // the  winding number counter

        // loop through all edges of the polygon
        for (unsigned int i=0; i<mEdges.size(); i++) {   // edge from V[i] to  V[i+1]
            const Point3D &P0=mVertices[mEdges[i][0]];
            const Point3D &P1=mVertices[mEdges[i][1]];
            if (P0.Y() <= P.Y()) {          // start y <= P.y
                if (P1.Y()  > P.Y()) {     // an upward crossing
                     if (P.isLeft( P0, P1) > 0)  // P left of  edge
                         ++wn;            // have  a valid up intersect
                }
            }
            else {                        // start y > P.y (no test needed)
                if (P1.Y()  <= P.Y())     // a downward crossing
                     if (P.isLeft(P0, P1) < 0)  // P right of  edge
                         --wn;            // have  a valid down intersect
            }
        }
        return wn;
    }

	//--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    // returns if pPoint lies on edge of indice #iEdg
    // computes non-signed distance pDist, and projected Point pProjP
    bool Polyline::closestPointToEdge(const Point3D & pPoint, unsigned int iEdg, double & pDist,
                                       Point3D & pProjP) const
	
    {
        const Point3D &P0=mVertices[mEdges[iEdg][0]];
        const Point3D &P1=mVertices[mEdges[iEdg][1]];

        Point3D V = P1 - P0;
        Point3D W = pPoint - P0;

        double c1 = W.dot(V);
        if ( c1 <= 0 ) {
            pProjP=P0;
            pDist=P0.distance(pPoint);
        }
        else {
            double c2 = V.dot(V);
            if ( c2 <= c1 ){
                pProjP=P1;
                pDist=P1.distance(pPoint);
            } else {
                double b = c1 / c2;
                Point3D Pb= P0 + b * V;
                pProjP=Pb;
                pDist=Pb.distance(pPoint);
            }
        }
        // sign of distance (negative if pP inside)
//        if (W.dot(mEdges[iEdg].getNormal())<0)
//            pDist=-pDist;

        //True if the given node is co-linear to this edge
        //return (c1<1E-8);

        //True if the given node lies on this edge
        return (pDist<1E-8);

    }

	//--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    bool Polyline::pointIsInMesh(const Point3D & pPoint ) const{
		// define if a point is inside a mesh or not

        if (mEdges.empty()) {
            return false;
        }
        //the winding number (=0 only when P is outside)
        return( windingNumber(pPoint)>0);
	}
	
	//--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    bool Polyline::pointIsInMesh(const Point3D & pPoint, const list<unsigned int> &lEdges) const{

        // index of the closest edge
        unsigned int closestEdge = 0;
        // closest point on the list of Edges (on edge, or vertice)
        Point3D pProjP,pProjP_tmp;
        // distance to this closest point (always positive)
        double pDist;
        //current closest distance: positive infinity
        double closestDist = numeric_limits<double>::infinity();
//		// true if this node is inside the surface
        bool bIsIn = false;
		
        if (mEdges.empty() || lEdges.empty()) {
            std::cerr << "Error at Polyline::pointIsInMesh, nowhere to project a point\n";
            exit(1); //or: return false;
        }
		
        bool found=false;
        // loop through all edges of the list
        for (unsigned int iEdg=0; iEdg<lEdges.size(); iEdg++) {

            bool isOn=closestPointToEdge(pPoint, iEdg, pDist, pProjP_tmp);

            if (pDist < closestDist) {
                pProjP = pProjP_tmp;
                closestEdge = iEdg;
                closestDist = pDist;
                found = true;
                if (isOn) // lies on edge, no need to continue
                    break;
            }

            if (found) {
                const Point3D &P0=mVertices[mEdges[closestEdge][0]];
                const Point3D &P1=mVertices[mEdges[closestEdge][1]];
                bIsIn = pProjP.isLeft(P0,P1);
            }
            else {
                std::cerr << "Error in Polyline::getProjection";
                std::cerr << " couldn't project node\n";
            }
        }
        return bIsIn;
    }

    //projection of a Point to Edge iedg
    Point3D Polyline::getProjection(const Point3D &pPoint, int iedg) const {
        {
            const Point3D &P0=mVertices[mEdges[iedg][0]];
            const Point3D &P1=mVertices[mEdges[iedg][1]];
            Point3D V = P1 - P0;
            Point3D W = pPoint - P0;

            double c1 = W.dot(V);
            if ( c1 <= 0 )
                return (P0);

            double c2 = V.dot(V);
            if ( c2 <= c1 )
                return (P1);

            double b = c1 / c2;
            Point3D Pb= P0 + b * V;
            return (Pb);
        }
    }

	//--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    Point3D Polyline::getProjection(const Point3D & pPoint) const{

        // FJA TODO merge with
        // Point3D Polyline::getProjection(const Point3D & pPoint, list<unsigned int> &lEdges) const;
        // set lEdges (list<>) with mEdges (vector<>) and call
        // return (pPoint, lEdges);

        // closest point on the edge (on edge, or vertice)
        Point3D pProjP_tmp,pProjP;
        // distance to this closest point (always positive)
        double pDist;
        //current closest distance: positive infinity
        double closestDist = numeric_limits<double>::infinity();

        if (mEdges.empty()) {
            cerr << "Error at Polyline::getProjection nowhere to project a point\n";
            exit(1); //or: return pProjP;
        }

        bool found = false;
        // browsing all the surface faces for min distance.
        for (unsigned int iEdge = 0; iEdge < mEdges.size(); iEdge++) {
            // computing the distance for this edge (segment)
            bool isOn=closestPointToEdge(pPoint,iEdge,pDist,pProjP_tmp);

            if (pDist < closestDist) {
                pProjP = pProjP_tmp;
                closestDist = pDist;
                found = true;
                if (isOn) // lies on edge, no need to continue
                    break;
            }
        }
        if (!found) {
            std::cerr << "Error in Polyline::getProjection";
            std::cerr << " couldn't project node\n";
        }

        return pProjP;
    }
		
	
	//--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    Point3D Polyline::getProjection(const Point3D & pPoint, list<unsigned int> &lEdges) const{

        // closest point on the edge (on edge, or vertice)
        Point3D pProjP_tmp,pProjP;
        // distance to this closest point (always positive)
        double pDist;
        //current closest distance: positive infinity
        double closestDist = numeric_limits<double>::infinity();

        if (mEdges.empty() || lEdges.empty()) {
            std::cerr << "Error at Polyline::getProjection nowhere to project a point\n";
            exit(1); //or: return pProjP;
        }

        bool found = false;
        // browsing all the surface faces for min distance.
        for (unsigned int iEdge = 0; iEdge < lEdges.size(); iEdge++) {
            // computing the distance for this edge (segment)
            bool isOn=closestPointToEdge(pPoint,iEdge,pDist,pProjP_tmp);

            if (pDist < closestDist) {
                pProjP = pProjP_tmp;
                closestDist = pDist;
                found = true;
                if (isOn) // lies on edge, no need to continue
                    break;
            }
        }
        if (!found) {
            std::cerr << "Error in Polyline::getProjection";
            std::cerr << " couldn't project node\n";
        }
        return pProjP;
    }

    Point3D Polyline::checkNormalToPlane() {
        Point3D W;
        if (mVertices.size()>=3) {
            Point3D U=mVertices[1]-mVertices[0], V;
            uint i=2;
            do { // suppose not all points are colinear....
                V=mVertices[i]-mVertices[0];
                ++i;
            } while (U.is2DCollinear(V));
//            W=(U^V).normalize();
        }

//        if (W.distance(Point3D(0.0,0.0,1.0)) >1E-8) {
        if ( ! W.is3DCollinear(Point3D(0.0,0.0,1.0)) ) {
            std::cerr << "warning Polyline::computeNormalToPlane()\n";
            std::cerr << "Polyline not aligned to XY plane\n";
            exit(1);
        }
        return (W);
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    vector<Point3D> Polyline::getNormals() const{
		vector<Point3D> normals;
        for (unsigned int i=0; i<getEdges().size(); i++){
            //FJA?? semiNormalized => /1000. ........... sure we want that ?
            normals.push_back(getEdges()[i].getSemiNormalizedNormal());
        }
        return normals;
	} 

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    ostream& operator<<(ostream& os, const Polyline &ply){
        for (auto &v:ply.getPoints()) {
            os << v << "; ";
        }
        os << std::endl;
        for (auto &n:ply.getNormals()) {
            os << n << "; ";
        }
        os << std::endl;
        for (auto &e:ply.getEdges()) {
            os << e << " ";
        }
        os << std::endl;
        return os;
    }
}
