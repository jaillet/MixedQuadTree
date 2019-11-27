/*
 <Mix-mesher: region type. This program generates a mixed-elements 2D mesh>

 Copyright (C) <2013,2019>  <Claudio Lobos> All rights reserved.

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
* @file Quadrant.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#ifndef Quadrant_h
#define Quadrant_h 1

#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <limits>
#include <algorithm>    // std::sort

#include "Visitors/Visitor.h"

#include "MeshPoint.h"
#include "Point3D.h"
#include "QuadEdge.h"

using Clobscode::MeshPoint;
using Clobscode::QuadEdge;
using Clobscode::Point3D;
using std::vector;
using std::list;
using std::set;
using std::pair;

namespace Clobscode
{
	
	class Quadrant{
        friend class IntersectionsVisitor;
        friend class IntersectionsDrawingVisitor;
        friend class OneIrregularVisitor;
        friend class PointMovedVisitor;
        friend class SplitVisitor;
        friend class EdgeVisitor;
        friend class TransitionPatternVisitor;
        friend class SurfaceTemplatesVisitor;
        friend class RemoveSubElementsVisitor;

	public:
		
		Quadrant(vector<unsigned int> &epts, const unsigned short &ref_level,
                 const unsigned int &q_id);
		
		virtual ~Quadrant();

        bool accept(Visitor *v);
		
		//access methods
        virtual const vector<unsigned int> &getPointIndex() const; // read only
        virtual unsigned int getPointIndex(unsigned int i) const;

        virtual const vector<unsigned int> getSubPointIndex() const; // read only
        virtual const vector<unsigned int> getSubPointIndexOnlyOnEdges() const; // read only

        virtual bool isInside() const;
		
        virtual bool intersectsSurface() const;
        
        virtual bool pointInside(const vector<MeshPoint> &mp, const Point3D &p) const;
		
        virtual unsigned short getRefinementLevel() const;
		
        virtual const list<unsigned int> &getIntersectedEdges() const; //read only
        virtual list<unsigned int> &getIntersectedEdges() ; //modification

        virtual const vector<vector<unsigned int> > &getSubElements() const; //read only
        virtual const vector<unsigned int> &getSubElement(unsigned int) const; //read only
        virtual       vector<unsigned int> &getSubElement(unsigned int); //modif

        virtual void updateSubElements(vector<vector <unsigned int> > &nsubs);
        virtual void removeSubElement(unsigned int i);

        virtual void computeMaxDistance(vector<MeshPoint> &mp);
		
        virtual double getMaxDistance() const;
        
        virtual bool badAngle(unsigned int &nIdx, const vector<MeshPoint> &mp) const;
		
		//flag for inside Quadrants that due to "inside node" moved
		//to the input domain, it must be treated as a surface
		//element by the surfacePatterns
		virtual void setSurface();
		
        virtual bool isSurface() const;
		
        virtual void setIntersectedEdges(list<unsigned int> &iedges);
        
        virtual void setIntersectedFeatures(const list<unsigned int> &iFeatures);
        virtual const list<unsigned int>& getIntersectedFeatures() const;
        virtual list<unsigned int>& getIntersectedFeatures();
        virtual bool hasIntersectedFeatures() const;
        
        //virtual void setInRegionState(const bool &value);
        //virtual const bool &isInRegion();
        
        virtual const unsigned int&getIndex();
        
        
        /***** BEGIN Debugging methods *******/
        virtual void setDebugging();
        
        virtual bool isDebugging() const;
        
        virtual double getAngle(unsigned int &nIdx, const vector<MeshPoint> &mp) const;
        /***** END Debugging methods *******/

    protected:
        
		//protected:
		vector<unsigned int> pointindex;
        vector<vector<unsigned int> > sub_elements; //, possibles, continuity;
        list<unsigned int> intersected_edges;
        list<unsigned int> intersected_features;
        //the level at which this Quadrant is found in the
        //the tree structure (Quadtree). Used for optimization
		unsigned short ref_level;
        
        //the quad unique identifier
        unsigned int q_id;
		
        //inregion is set to true when this quadrant, and therfore their childs,
        //intersects a RefinementRegion, allowing to avoid asking again.
        bool surface;//, inregion;
		
        /***** BEGIN Debugging variables *******/
        bool debugging;
        /***** END Debugging variables *******/
        
		double max_dis;
	};
	
    /***** BEGIN Debugging methods *******/
    inline void Quadrant::setDebugging() {
        debugging = true;
    }
    
    inline bool Quadrant::isDebugging() const {
        return debugging;
    }
    /***** END Debugging methods *******/
    
    /*inline void Quadrant::setInRegionState(const bool &value) {
        inregion = value;
    }
    
    inline const bool &Quadrant::isInRegion() {
        return inregion;
    }*/
    
    
    inline void Quadrant::setIntersectedFeatures(const list<unsigned int> &iFeatures) {
        intersected_features = iFeatures;
    }
    
    inline const list<unsigned int>& Quadrant::getIntersectedFeatures() const {
        return intersected_features;
    }
    inline list<unsigned int>& Quadrant::getIntersectedFeatures() {
        return intersected_features;
    }

    inline bool Quadrant::hasIntersectedFeatures() const {
        return (intersected_features.size()>0);
    }

    inline const vector<unsigned int> &Quadrant::getPointIndex() const{
		return pointindex;
	}
    inline unsigned int Quadrant::getPointIndex(unsigned int i) const {
        return pointindex[i];
    }

    inline const vector<unsigned int> Quadrant::getSubPointIndex() const{
        vector<unsigned int> subpointindex;

        if (sub_elements.size()==1) { //not subdivided, return corners
            return getPointIndex();
        } else {
            //go through each subelement, construct inner node list
            for (const auto & subelem:sub_elements ) {
                subpointindex.insert(subpointindex.end(),subelem.begin(),subelem.end());
            }
            //remove duplicate
            sort(subpointindex.begin(),subpointindex.end());
//            auto it=unique(subpointindex.begin(),subpointindex.end());
//            subpointindex.resize( std::distance(subpointindex.begin(),it) );
            //std::unique doesn't remove duplicates, just throws them to the end
            //that's the reason for std::erase...
            subpointindex.erase(unique(subpointindex.begin(),subpointindex.end()),subpointindex.end());
            subpointindex.shrink_to_fit();

            return subpointindex;
        }
    }

    inline const vector<unsigned int> Quadrant::getSubPointIndexOnlyOnEdges() const {

        //first, get all index from subelems
        vector<unsigned int> subpointindex(getSubPointIndex());

        // next circulate the corner's index, and remove from
        for (auto i:pointindex) {
            //surprisingly, std::remove doesn't remove duplicates, just throws them to the end...
            //that's the reason for std::erase...
            subpointindex.erase(std::remove(subpointindex.begin(), subpointindex.end(), i),
                                subpointindex.end());
        }

        return subpointindex;
    }


    inline bool Quadrant::isInside() const {
        return intersected_edges.empty();
	}
	
    inline bool Quadrant::intersectsSurface() const {
        return !intersected_edges.empty();
	}
	
    inline const list<unsigned int> &Quadrant::getIntersectedEdges() const {
        return intersected_edges;
    }
    inline list<unsigned int> &Quadrant::getIntersectedEdges() {
        return intersected_edges;
    }

    inline unsigned short Quadrant::getRefinementLevel() const {
		return ref_level;
	}
	
    inline const vector<vector<unsigned int> > &Quadrant::getSubElements() const {
        return sub_elements;
    }
    inline const vector<unsigned int> &Quadrant::getSubElement(unsigned int i) const {
        return sub_elements[i];
    }
    inline vector<unsigned int> &Quadrant::getSubElement(unsigned int i) {
        return sub_elements[i];
    }

    inline void Quadrant::computeMaxDistance(vector<MeshPoint> &mp){
        const Point3D &p0 = mp[pointindex[0]].getPoint();
        const Point3D &p1 = mp[pointindex[2]].getPoint();
        //max_dis = 0.3 * (p0 - p1).Norm();
        max_dis = 0.3 * (p0 - p1).Norm() / sqrt(getSubElements().size());
    }
    
    inline void Quadrant::updateSubElements(vector<vector <unsigned int> > &nsubs) {

        sub_elements = nsubs;
        return;

        /*cout << "updateSubElements: " << sub_elements.size();
        sub_elements.clear();
        sub_elements.reserve(nsubs.size());
        for (auto ne:nsubs) {
            sub_elements.push_back(ne);
        }
        cout << " -> " << sub_elements.size() << "\n";*/
    }

    inline void Quadrant::removeSubElement(unsigned int i) {

        sub_elements.erase(sub_elements.begin()+i );
    }

    inline double Quadrant::getMaxDistance() const {
		return max_dis;
	}
	
    inline void Quadrant::setSurface() {
		surface = true;
	}
	
    inline bool Quadrant::isSurface() const {
        return surface || !intersected_edges.empty();
	}
	
    inline void Quadrant::setIntersectedEdges(list<unsigned int> &iedges){
        intersected_edges = iedges;
	}
    
    inline const unsigned int&Quadrant::getIndex() {
        return q_id;
    }
	
    std::ostream& operator<<(ostream& o, const Quadrant &q);
}
#endif
