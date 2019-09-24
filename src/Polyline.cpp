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
#include "GeometricTransform.h"

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
        
        computeNodesAngle();
		
	}

    Polyline::~Polyline(){
		
	}
	
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void Polyline::computeBounds(){

        bounds.resize(6,0.0);

        // .at(0) throw exception if empty polyline
        bounds[0]=bounds[3]=mVertices[0].X();
        bounds[1]=bounds[4]=mVertices[0].Y();
        bounds[2]=bounds[5]=mVertices[0].Z();

        // search for max and min coodinates in X, Y and Z.
        double x,y,z;
        for (unsigned int i=1; i< mVertices.size(); ++i) {

            x = mVertices[i].X();
            y = mVertices[i].Y();
            z = mVertices[i].Z();

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
    void Polyline::computeNodesAngle(){
        
        unsigned int npts = mVertices.size();
        mVerticesAngles.assign(npts,0);
        vector<vector<unsigned int> > seg_per_node(npts);
        
        //save a reference to the segment for each node
        for (unsigned int i=0; i<mEdges.size(); i++) {
            seg_per_node[mEdges[i][0]].push_back(i);
            seg_per_node[mEdges[i][1]].push_back(i);
        }
        
        //compute normal of each node
        for (unsigned int i=0; i<npts; i++) {
            if (seg_per_node[i].size()!=2) {
                continue;
            }
            unsigned int i0, i2;
            i0 = mEdges[seg_per_node[i][0]][1];
            i2 = mEdges[seg_per_node[i][1]][0];
            if (i2==i0) {
                i0 = mEdges[seg_per_node[i][1]][1];
                i2 = mEdges[seg_per_node[i][0]][0];
            }
            const Point3D &P1 = mVertices[i],
                    &P0 = mVertices[i0],
                    &P2 = mVertices[i2];
            
            // take care, polyline is oriented reversely to quadrant ;o(
            mVerticesAngles[i]= P1.angle3Points(P2,P0);

//            cerr << "Poly[" << i << "] angle = " << i0  << " " << i  << " " << i2 << " " << mVerticesAngles[i] << "\n";
            //normalize
            //mVerticePseudoNormals[i].normalize();
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

//            cerr << "Poly[" << i << "] normal = " << mVerticePseudoNormals[i] << "\n";
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
                                       Point3D & pProjP, int &projWhere) const
	
    {
        const Point3D &P0=mVertices[mEdges[iEdg][0]];
        const Point3D &P1=mVertices[mEdges[iEdg][1]];
        projWhere=-1; // default, if projection on the edge

        Point3D V = P1 - P0;
        Point3D W = pPoint - P0;

        double c1 = W.dot(V);
        if ( c1 <= 0 ) {
            pProjP=P0;
            pDist=P0.distance(pPoint);
            projWhere=mEdges[iEdg][0]; // if projection on P0
        }
        else {
            double c2 = V.dot(V);
            if ( c2 <= c1 ){
                pProjP=P1;
                pDist=P1.distance(pPoint);
                projWhere=mEdges[iEdg][1]; // if projection on P1
            } else {
                double b = c1 / c2;
                pProjP=P0 + b * V;
                pDist=pProjP.distance(pPoint); // proj on edge
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

        // index of the closest edges (if many)
        vector <unsigned int> closestEdges;
        // closest point on the list of Edges (on edge, or vertice)
        Point3D pProjP,pProjP_tmp;
        // -1 if projection on edge, index of the point if projection on an extremity
        int projWhere=-1, projWhere_tmp;
        // distance to this closest point (always positive)
        double pDist;
        //current closest distance: positive infinity
        double closestDist = numeric_limits<double>::infinity();
//		// true if this node is inside the surface
        bool bIsIn = false;
		
        if (mEdges.empty() || lEdges.empty()) {
            std::cerr << "Error at Polyline::pointIsInMesh, nowhere to project a point\n";
            exit(EXIT_FAILURE); //or: return false;
        }
		
        bool found=false;
        // loop through all edges of the list
        for (auto iEdg:lEdges) {

            bool isOn=closestPointToEdge(pPoint, iEdg, pDist, pProjP_tmp, projWhere_tmp);

            // if projection on the same point for 2 edges
            if (found && projWhere_tmp>=0 && projWhere_tmp==projWhere ) {
                closestEdges.push_back(iEdg); //just update list
            } else {
                if (pDist < closestDist ) {
                    pProjP = pProjP_tmp;
                    closestEdges.clear();
                    closestEdges.push_back(iEdg);
                    closestDist = pDist;
                    projWhere=projWhere_tmp;
                    found = true;
                    if (isOn) // lies on edge, no need to continue
                        return false; //break;
                }
            }
        }
        if (found) {

            uint iClosest=0;

            if (closestEdges.size()>1) { // if closest point lies on an edge extremity
                Point3D v=(pPoint-pProjP).normalize();
                double dot0= v.dot(mEdges[closestEdges[0]].getNormalizedNormal());
                double dot1= v.dot(mEdges[closestEdges[1]].getNormalizedNormal());
                // check with the edge that normal is most colinear to the (pPoint-pProjP) line
                // in other terms, the direction of the projection is most orthogonal to that edge
                if (fabs(dot0)>=fabs(dot1)) {
                    iClosest=0; // first edge
                } else {
                    iClosest=1; // second edge
                }
            }

            const Point3D &P0=mVertices[mEdges[closestEdges[iClosest]][0]];
            const Point3D &P1=mVertices[mEdges[closestEdges[iClosest]][1]];
            //bIsIn = (pPoint.isLeft(P0,P1)>=0); // pPoint left of the closest edge, "in" if "on"
            bIsIn = (pPoint.isLeft(P0,P1)>0); // pPoint left of the closest edge, "out" if "on"
            //                std::cerr << pPoint << "  e:" << mEdges[closestEdge[0]] << " dot: " << dot0 << "  -  " << mEdges[closestEdge[1]] << " dot: " << dot1 << "->" << (fabs(dot0)>=fabs(dot1))
            //                          << " -- " << bIsIn << "=" << windingNumber(pPoint) << "\n" << std::flush;

            //FJA temporary test to get removed once one is sure this is working
            /*if (bIsIn != pointIsInMesh(pPoint)) {
                std::cerr << "big bug here\n";
                //FJA in case this is not working, back to the computation for the whole polyline ;o(
                //return(pointIsInMesh(pPoint));
                exit(EXIT_FAILURE);
            }*/

            return bIsIn;
        }
        else {
            std::cerr << "Error in Polyline::getProjection";
            std::cerr << " couldn't project node\n";
        }

        return bIsIn;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    inline bool Polyline::pointIsFeature(unsigned int i) const {
        return (mVerticesAngles[i]<MinFeatureAngle || mVerticesAngles[i]>MaxFeatureAngle);
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    list<unsigned int> Polyline::getFeatureProjection(const Quadrant &q,
                                                 vector<MeshPoint> &mp) const {
        
        list<unsigned int> iEdges = q.getIntersectedEdges();
        list<unsigned int> fes;
        if (iEdges.size()<2) {
            return fes;
        }
        
        vector<unsigned int> iNodes;
        for (unsigned int iEdg:iEdges) {
            iNodes.push_back(mEdges[iEdg][0]);
            iNodes.push_back(mEdges[iEdg][1]);
        }
        std::sort(iNodes.begin(),iNodes.end());
        
        //FJA: we want to test only the intersection, not the extremities of the edges
        // the feature should lies IN the quadrant
        list<unsigned int> iCommonNodes;
        for (unsigned int i=1; i<iNodes.size(); ++i) {
            if (iNodes.at(i-1)==iNodes.at(i)){
                iCommonNodes.push_back(iNodes[i]); //insert duplicate
                ++i; //skip next duplicate item
            }
        }
        
        for (unsigned int iNd:iCommonNodes) {
            //if angle is between 175 and 185 grades, then is not a sharp feature:
            if ( pointIsFeature(iNd) ) {
                //std::cout << " " << iNd << std::flush;
                fes.push_back(iNd);
            }
        }
        return fes;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    unsigned int Polyline::getNbFeatures(Quadrant &q, const vector<MeshPoint> &mp) const {
        unsigned int nbFeat=0;

        const list<unsigned int> &iEdges = q.getIntersectedEdges();
        list<unsigned int> iFeatures;

        if (iEdges.size()<2) {
            return false;
        }

        vector<unsigned int> iNodes;
        for (const auto iEdg:iEdges) {
            iNodes.push_back(mEdges[iEdg][0]);
            iNodes.push_back(mEdges[iEdg][1]);
        }
        std::sort(iNodes.begin(),iNodes.end());

        //FJA: we want to test only the intersection, not the extremities of the edges
        // the feature should lie IN the quadrant
        list<unsigned int> iCommonNodes;
        for (unsigned int i=1; i<iNodes.size(); ++i) {
            if (iNodes.at(i-1)==iNodes.at(i)){
                iCommonNodes.push_back(iNodes[i]); //insert duplicate
                ++i; //skip next duplicate item
           }
        }
        
        for (auto iNd:iCommonNodes) {
            //if angle is between 175 and 185 grades, then is not a sharp feature:
            if ( pointIsFeature(iNd) ) {
                //std::cout << " " << iNd << std::flush;
                if (q.pointInside(mp,mVertices[iNd])) {
                    iFeatures.push_back(iNd);
                    ++nbFeat;
                }
            }
        }
        q.setIntersectedFeatures(iFeatures);
        return nbFeat;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
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
        // -1 if projection on edge, index of the point if projection on an extremity
        int projWhere_tmp;
        // distance to this closest point (always positive)
        double pDist;
        //current closest distance: positive infinity
        double closestDist = numeric_limits<double>::infinity();

        if (mEdges.empty()) {
            cerr << "Error at Polyline::getProjection nowhere to project a point\n";
            exit(EXIT_FAILURE); //or: return pProjP;
        }

        bool found = false;
        // browsing all the surface faces for min distance.
        for (unsigned int iEdge = 0; iEdge < mEdges.size(); iEdge++) {
            // computing the distance for this edge (segment)
            bool isOn=closestPointToEdge(pPoint,iEdge,pDist,pProjP_tmp,projWhere_tmp);

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
    Point3D Polyline::getProjection(const Point3D & pPoint, const list<unsigned int> &lEdges) const{

        // closest point on the edge (on edge, or vertice)
        Point3D pProjP_tmp,pProjP;
        // -1 if projection on edge, index of the point if projection on an extremity
        int projWhere, projWhere_tmp;
        // distance to this closest point (always positive)
        double pDist;
        //current closest distance: positive infinity
        double closestDist = numeric_limits<double>::infinity();

        if (mEdges.empty() || lEdges.empty()) {
            std::cerr << "Error at Polyline::getProjection nowhere to project a point\n";
            exit(EXIT_FAILURE); //or: return pProjP;
        }

        bool found = false;
        // browsing all the surface faces for min distance.
        for (unsigned int iEdge:lEdges) {
            // computing the distance for this edge (segment)
            bool isOn=closestPointToEdge(pPoint,iEdge,pDist,pProjP_tmp,projWhere_tmp);

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

    Point3D Polyline::checkNormalToPlane() const {
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
            exit(EXIT_FAILURE);
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
        for (const auto &v:ply.getPoints()) {
            os << v << "; ";
        }
        os << std::endl;
        for (const auto &n:ply.getNormals()) {
            os << n << "; ";
        }
        os << std::endl;
        for (const auto &e:ply.getEdges()) {
            os << e << " ";
        }
        os << std::endl;
        return os;
    }
}
