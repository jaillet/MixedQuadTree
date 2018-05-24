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
                       vector<vector<unsigned int> > &edg_ind){
		
		if (pts.empty()) {
			std::cout << "Error in Mesh::Init input mesh without points\n";
			return;
		}
		
        // init bouding box with first Point3D
        bounds.resize(6);
		bounds[0]=bounds[3]=pts[0][0];
		bounds[1]=bounds[4]=pts[0][1];
		bounds[2]=bounds[5]=pts[0][2];
		
		
        // initialising the vertices
        // and also, search for max and min coodinates in X, Y and Z.
        mVertices.reserve(pts.size());
        mVertices.push_back(pts[0]);
        for (unsigned int i=1; i< pts.size(); ++i) {
			
			double x = pts[i][0];
			double y = pts[i][1];
			double z = pts[i][2];
			
            mVertices.push_back(pts[i]);

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
        for (int i=0; i<mEdges.size(); i++) {    // edge from V[i]  to V[i+1]
            const Point3D &P0=mVertices[mEdges[i][0]];
            const Point3D &P1=mVertices[mEdges[i][1]];
            if (((P0.Y() <= P.Y()) && (P1.Y() > P.Y()))     // an upward crossing
            || ((P0.Y() > P.Y()) && (P1.Y() <=  P.Y()))) { // a downward crossing
                // compute  the actual edge-ray intersect x-coordinate
                float vt = (float)(P.Y()  - P1.Y()) / (P1.Y() - P0.Y());
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
        for (int i=0; i<mEdges.size(); i++) {   // edge from V[i] to  V[i+1]
            const Point3D &P0=mVertices[mEdges[i][0]];
            const Point3D &P1=mVertices[mEdges[i][1]];
            if (P0.Y() <= P.Y()) {          // start y <= P.y
                if (P1.Y()  > P.Y())      // an upward crossing
                     if (P.isLeft( P0, P1) > 0)  // P left of  edge
                         ++wn;            // have  a valid up intersect
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
    bool Polyline::SignedDistToSegment(const Point3D & pP, const unsigned int &iEdg,
                                       const double &current_min_dist, double & pDist,
                                       Point3D & pProjP, bool & pIsIn, int & edgeNode) const
	
    {
        std::cerr << "bool Polyline::SignedDistToSegment(...) const\n" ;
        std::cerr << "not implemented yet\n" ;

        // compute the shortest distance between pP and the Segment (index iT).
		
        // inputs:
        // Point3D pP, the point to be projected
        // uInt iEdg, the index of the edge/segment to test
		
        // outputs
        // double & pDist, the signed distance between the point and the edge
        // VPoint3D & pProjP, the projected point
        // bool pIsIn, boolean expressing if the node P is inside or outside the polyline
		
		
        // declare points and normals of this edge.
        // p* the vertices of the edge
        // pN*, the pseudo normals at each vertices of the edge
        // pN**, the pseudo normals at each edges of the polyline
        Point3D pA = mVertices[mEdges[iEdg][0]];
        Point3D pB = mVertices[mEdges[iEdg][1]];
		
        Point3D pNA = mVerticePseudoNormals[mEdges[iEdg][0]];
        Point3D pNB = mVerticePseudoNormals[mEdges[iEdg][1]];
		
        // Edge vectors
        Point3D AB = pB - pA;
		
//        // Normal vector pointing towards positive half-space, assuming (ABC) is CCW
//        Point3D N  = AB^AC;
//        double N2  = N * N;
		
        // normal used to compute the sign
        Point3D pProjN;
		
        //True if the given node is co-linear to this edge
        bool coplanar = false;
		
//        //------------------------------
//        // Projection on plane (ABC)
//        //------------------------------
		
//        // P1 = P projected on plane
//        Point3D P1;
//        Point3D AP1;
//        {
//            // Projection of P on ABC plane
//            Point3D AP = pP - pA;
//            double k = ( AP * N ) / N2;
			
//			// do not trust this face when the node is co-planar
//			// to it and the result is "outside", unless all
//			// faces be co-planar to the node.
//			if (fabs(k)<1E-8) {
//				coplanar = true;
//			}
			
//            P1 = pP - k*N;
						
//            // Check if P1 is in triangle
//            AP1 = P1 - pA;
//            double x = (AP1 ^  AC) * N / N2;
//            double y = (AB  ^ AP1) * N / N2;
			
//            // AP1 = x*AB + y*AC
//            if( x>=0 && y>=0 && x+y<=1 )
//            {
//                // Normal interpolation
//				// real normal of the triangle
//				pProjN = N;
				
//                // Signed distance
//                double lAbsDist = pP.distance(P1);
//                double lSgn = pProjN * ( pP - P1 );
				
//                // Projection params
//                pProjP = P1;
//                pDist = lSgn<0 ? -lAbsDist : lAbsDist ;
				
//				// computing if the node is inside the surface
//				pIsIn = lSgn < 0;
				
//				faceEdgeNode = 0;
                
//				// return false as the node is inside the triangle
//				return false;
//            }
//        }
		
//        //-------------------------------------------------------------------
//        // P1 is out of the triangle : compute projections on triangle edges
//        //-------------------------------------------------------------------
		
//        //-------------------------------
//        // Distance to segment AB
//        //-------------------------------
//        {
//            Point3D AP = P - pA;
//            double AB2 = AB * AB;
//            double t = ( AP * AB ) / AB2;
			
//            if( t < 0 )
//            {
//                pProjP = pA;
//                pProjN = pNA;
				
//                faceEdgeNode = 2;
//            }
//            else
//                if( t > 1 )
//                {
//                    pProjP = pB;
//                    pProjN = pNB;
					
//                    faceEdgeNode = 2;
//                }
//                else
//                {
//                    pProjP = pA + t * AB;
//                    // pseudo normal of the edge AB
//                    pProjN = pNab;
//                    faceEdgeNode = 1;
//                }
			
//            pDist = pProjP.distance( pP );
//        }

//        //-------------------------------
//        // Orientation
//        //-------------------------------
		
//        // Change distance sign
//        if( pProjN * ( pP - pProjP ) < 0)
//		{
//            pDist = -pDist;
//		}
		
//		pIsIn = pDist < 1E-8;
		
        return coplanar;
		
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
    bool Polyline::pointIsInMesh(const Point3D & pPoint, list<unsigned int> &lFaces) const{
        std::cerr << "bool Polyline::pointIsInMesh(const Point3D & pPoint, list<unsigned int> &lFaces) const\n" ;
        std::cerr << "not implemented yet\n" ;

//		// define if a point is inside a mesh or not
		
//		// index of the closest triangle
//		unsigned int closestTriangle = 0;
//		// closest point on the triangle (on triangle face, on edge, or vertice)
//		Point3D pProjP;
//		// distance to this closest point (always positive)
//		double pDist;
//		//current closest distance: positive infinity
//		double closestDist = numeric_limits<double>::infinity();
//		// true if this node is inside the surface
//		bool pIsIn;
        bool bIsIn = false;
//		// 0 if close to a face, 1 if close to an edge, 2 if close to a vertice
//		int faceEdgeNode;
//		int iFaceEdgeNode;
		
//		if (mTriangles.empty() || lFaces.empty()) {
//			return false;
//		}
		
//		list<unsigned int>::iterator iSurfF;
		
//		//one_good recalls if a non co-planar face conserve the
//		//current_min_dis
//		bool one_good = false;
		
//		// browsing all the surface faces
//		for (iSurfF = lFaces.begin(); iSurfF!=lFaces.end(); iSurfF++)
//		{
			
//			bool coplanar = SignedDistToTriangle(pPoint,*iSurfF,closestDist,
//												 pDist,pProjP,pIsIn,faceEdgeNode);
			
//			if (coplanar) {
//				if (!one_good) {
//					pDist = fabs(pDist);
					
//					if (pDist < closestDist) {
//						closestTriangle = *iSurfF;
//						closestDist = pDist;
//						bIsIn = pIsIn;
//						iFaceEdgeNode = faceEdgeNode;
//					}
//				}
//			}
//			else {
//				pDist = fabs(pDist);
				
//				if (!one_good || pDist < closestDist) {
//					closestTriangle = *iSurfF;
//					closestDist = pDist;
//					bIsIn = pIsIn;
//					iFaceEdgeNode = faceEdgeNode;
//				}
//				one_good = true;
//			}
//		}
		
        return bIsIn;
    }
	
	//--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    Point3D Polyline::getProjection(const Point3D & pPoint) const{
        std::cerr << "Point3D Polyline::getProjection(const Point3D & pPoint) const\n" ;
        std::cerr << "not implemented yet\n" ;

//		// define if a point is inside a mesh or not
		
//		// index of the closest triangle
//		unsigned int closestTriangle = 0;
//		// closest point on the triangle (on triangle face, on edge, or vertice)
        Point3D pProjP_tmp,pProjP;
//		// distance to this closest point (always positive)
//		double pDist;
//		//current closest distance: positive infinity
//		double closestDist = numeric_limits<double>::infinity();
//		// true if this node is inside the surface
//		bool pIsIn;
//		bool bIsIn = false;
//		// 0 if close to a face, 1 if close to an edge, 2 if close to a vertice
//		int faceEdgeNode;
//		int iFaceEdgeNode;
		
//		if (mTriangles.empty()) {
//			cout << "Error at Polyline::getProjection nowhere to project a point\n";
//			return pProjP;
//		}
		
//		// browsing all the surface faces for min distance.
//		for (unsigned int iSurfF = 0; iSurfF < mTriangles.size(); iSurfF++) {
//			// computing the distance for this face (triangle)
//            SignedDistToEdge(pPoint,iSurfF,closestDist,pDist,pProjP_tmp,pIsIn,faceEdgeNode);
			
//			pDist = fabs(pDist);
			
//			if (pDist < closestDist) {
//				pProjP = pProjP_tmp;
//				closestTriangle = iSurfF;
//				closestDist = pDist;
//				bIsIn = pIsIn;
//				iFaceEdgeNode = faceEdgeNode;
//			}
//		}
		
        return pProjP;
    }
		
	
	//--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    Point3D Polyline::getProjection(const Point3D & pPoint, list<unsigned int> &lFaces) const{
        std::cerr << "Point3D Polyline::getProjection(const Point3D & pPoint, list<unsigned int> &lFaces)\n" ;
        std::cerr << "not implemented yet\n" ;

//        // define if a point is inside a mesh or not
		
//		// index of the closest triangle
//		unsigned int closestTriangle = 0;
//		// closest point on the triangle (on triangle face, on edge, or vertice)
//		Point3D pProjP_tmp,pProjP;
//		// distance to this closest point (always positive)
//		double pDist;
//		//current closest distance: positive infinity
//		double closestDist = numeric_limits<double>::infinity();
//		// true if this node is inside the surface
//		bool pIsIn;
//		bool bIsIn = false;
//		// 0 if close to a face, 1 if close to an edge, 2 if close to a vertice
//		int faceEdgeNode;
//		int iFaceEdgeNode;
		
//		if (mTriangles.empty() || lFaces.empty()) {
//			cout << "Error at Polyline::getProjection nowhere to project a point\n";
//			return pProjP;
//		}
		
//		list<unsigned int>::iterator iSurfF;
		
//		bool found = false;
		
//		// browsing all the surface faces
//		for (iSurfF = lFaces.begin(); iSurfF!=lFaces.end(); iSurfF++)
//		{
			
//			// computing the distance for this face (triangle)
//			SignedDistToTriangle(pPoint,*iSurfF,closestDist,pDist,pProjP_tmp,pIsIn,faceEdgeNode);
			
//			pDist = fabs(pDist);
			
//			if (pDist < closestDist) {
//				pProjP = pProjP_tmp;
//				closestTriangle = *iSurfF;
//				closestDist = pDist;
//				bIsIn = pIsIn;
//				iFaceEdgeNode = faceEdgeNode;
//				found = true;
//			}
//		}
		
//		if (!found) {
//			cout << "Error in Polyline::getProjection";
//			cout << " couldn't project node\n";
//		}
		
//		return pProjP;
	}
	
    vector<Point3D> Polyline::getNormals() const{
		vector<Point3D> normals;
        for (unsigned int i=0; i<getEdges().size(); i++){
            //FJA?? semiNormalized => /1000. ........... sure we want that ?
            normals.push_back(getEdges()[i].getSemiNormalizedNormal());
        }
        return normals;
	} 
}
