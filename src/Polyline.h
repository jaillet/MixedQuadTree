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
#include <limits>
#include <set>
#include <math.h>

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
                vector<vector<unsigned int> > &edg_ind);

        virtual ~Polyline();
		
        virtual vector<Point3D> &getPoints();
        virtual const vector<Point3D> &getPoints() const;

        virtual vector<Point3D> &getVerticePseudoNormals();
        virtual const vector<Point3D> &getVerticePseudoNormals() const;

        virtual vector<double> &getBounds();
        virtual const vector<double> &getBounds() const;

        virtual vector<PolyEdge> &getEdges();
        virtual const vector<PolyEdge> &getEdges() const;

        // 0 = outside, 1 = inside
        virtual int crossingNumber(const Point3D &P) const;
        // =0 only when P is outside
        virtual int windingNumber(const Point3D &P) const;

        virtual bool pointIsInMesh(const Point3D &pPoint) const;
		virtual bool pointIsInMesh(const Point3D & pPoint, 
                                   const list<unsigned int> &lEdges) const;
		
        //returns true if the list of edges containts at least one feature.
        virtual bool hasFeature(const list<unsigned int> &iEdges) const;

        //projection of a Point to Edge iedg
        virtual Point3D getProjection(const Point3D &pPoint, int iedg) const;

        //projection of a Point to the Polyline
        virtual Point3D getProjection(const Point3D &pPoint) const;

        //projection of a Point to a list of Edges, identified by index
		virtual Point3D getProjection(const Point3D & pPoint, 
                                   list<unsigned int> &lEdges) const;
								 
        virtual void update();

        virtual Point3D checkNormalToPlane();

        virtual Point3D getCentroid() const;

        // compute the pseudo normal at each surface node
        virtual void computeBounds();

        // compute the pseudo normal at each surface node
        virtual void computeEdgesNormal();
        
        // compute the angle of incident edges for each node
        virtual void computeNodesAngle();

        virtual vector<Point3D> getNormals() const;

        // compute the pseudo normal at each surface node
        virtual void computeNodesPseudoNormal();

        friend ostream& operator<<(ostream& os, const Polyline &ply);
		
	protected:
				
        // returns if pPoint co-linear to edge of indice #iEdg
        // computes non-signed distance pDist, and projected Point pProjP
        virtual bool closestPointToEdge(const Point3D & pPoint,
                                          unsigned int iEdg,
                                          double & pDist,
                                          Point3D & pProjP) const;
		
		
	protected:
		
		vector<Point3D> mVertices;
        vector<PolyEdge> mEdges;
		
		//one pseudo normal per vertices (same size)
		vector<Point3D> mVerticePseudoNormals;
        vector<double> mVerticesAngles;
		
        //bounding box of this surface mesh
        vector<double> bounds;
    };
	
    inline void Polyline::update() {
        computeBounds();
        checkNormalToPlane();
        computeEdgesNormal();
        computeNodesPseudoNormal();
    }

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

    inline vector<PolyMesh::PolyEdge> &Polyline::getEdges() {
        return mEdges;
    }
    inline const vector<PolyMesh::PolyEdge> &Polyline::getEdges() const {
        return mEdges;
    }

    inline vector<Point3D> &Polyline::getVerticePseudoNormals() {
        return mVerticePseudoNormals;
    }
    inline const vector<Point3D> &Polyline::getVerticePseudoNormals() const {
        return mVerticePseudoNormals;
    }

    inline Point3D Polyline::getCentroid() const {
        Point3D mCentroid ( (bounds[0]+bounds[3])/2.,
                            (bounds[1]+bounds[4])/2.,
                            (bounds[2]+bounds[5])/2.);;
		return mCentroid;
	}

}
#endif
