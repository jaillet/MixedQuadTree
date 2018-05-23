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
* @file Polyline.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#ifndef Polyline_h
#define Polyline_h 1

#include "Point3D.h"
#include "PolyEdge.h"
//FJA no more used #include "SurfEdgeContainer.h"
#include "SurfaceEdge.h"
#include <limits>
#include <set>

using std::vector;
using std::set;
using Clobscode::Point3D;
using PolyMesh::PolyEdge;
//FJA?? using SurfMesh::SurfEdgeContainer;

namespace Clobscode
{
    class Polyline {

    public:
        Polyline();
        
        Polyline(vector<Point3D> &pts,
				vector<vector<unsigned int> > &fv);
		
        virtual ~Polyline();
		
        virtual vector<Point3D> &getPoints();
        virtual const vector<Point3D> &getPoints() const;

        virtual vector<Point3D> &getVerticePseudoNormals();
        virtual const vector<Point3D> &getVerticePseudoNormals() const;

        virtual set<SurfaceEdge> &getEdges();
        virtual const set<SurfaceEdge> &getEdges() const;

        virtual vector<double> &getBounds();
        virtual const vector<double> &getBounds() const;

		virtual bool pointIsInMesh(const Point3D &pPoint);
		
		virtual bool pointIsInMesh(const Point3D & pPoint, 
								   list<unsigned int> &lFaces);
		
        virtual Point3D getProjection(const Point3D &pPoint) const;
		
		virtual Point3D getProjection(const Point3D & pPoint, 
                                   list<unsigned int> &lFaces) const;
								 
        virtual Point3D getCentroid() const;
		
        virtual vector<Point3D> getNormals() const;
		
	protected:
		
		// compute the pseudo normal at each surface node
		virtual void computeNodePseudoNormal();
		
		// compute the pseudo normal at each surface edge
		virtual void computeEdgePseudoNormal();
		
		virtual bool SignedDistToTriangle(const Point3D & pP, 
										  const unsigned int &iT, 
										  const double &current_min_dist, 
										  double & pDist, 
										  Point3D & pProjP, bool & pIsIn, 
										  int & faceEdgeNode);
		
		
	protected:
		
		vector<Point3D> mVertices;
        vector<PolyEdge> mSegments;
		
		//one pseudo normal per vertices (same size)
		vector<Point3D> mVerticePseudoNormals; 
		
		//one pseudo normal per edge, so 3 pseudo normals per 
		//triangle (3*size)
		vector<Point3D> mEdgePseudoNormals;
		
		//bounding box of this surface mesh
		vector<double> bounds;	
		
		//Point used by Ray Tracing Method (unstable)
		Point3D outside;
        
        //the set of edges
        //FJA to be decided: the best dataset for this...
        set<SurfaceEdge> edges;
        
        
        
        
	};
	
    inline vector<Point3D> &Polyline::getPoints() {
        return mVertices;
    }
    inline const vector<Point3D> &Polyline::getPoints() const {
        return mVertices;
    }

    inline vector<double> &Polyline::getBounds() {

        //cout << "tri bounds: " << bounds.size() << endl;
        return bounds;
    }
    inline const vector<double> &Polyline::getBounds() const {

        //cout << "tri bounds: " << bounds.size() << endl;
        return bounds;
    }

    inline set<SurfaceEdge> &Polyline::getEdges() {
        return edges;
    }
    inline const set<SurfaceEdge> &Polyline::getEdges() const {
        return edges;
    }

    inline vector<Point3D> &Polyline::getVerticePseudoNormals() {
        return mVerticePseudoNormals;
    }
    inline const vector<Point3D> &Polyline::getVerticePseudoNormals() const {
        return mVerticePseudoNormals;
    }

    inline Point3D Polyline::getCentroid() const {
		Point3D mCentroid;
		mCentroid.X() = (bounds[0]+bounds[3])/2;
        mCentroid.Y() = (bounds[1]+bounds[4])/2;
        mCentroid.Z() = (bounds[2]+bounds[5])/2;
		return mCentroid;
	}
	   
}
#endif
