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
* @file MeshPoint.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#ifndef MeshPoint_h
#define MeshPoint_h 1

#include "Point3D.h"
#include "Utils.h"

#include <iostream>
#include <math.h>
#include <vector>
#include <list>
#include <limits>

using Clobscode::Point3D;
using std::list;
using std::vector;

namespace Clobscode
{

typedef enum {
    STATEMASK = 0x0000,
    INSIDE =  1,
    FEATURE =  2,
    PROJECTED =  4,
    OUTCHECKED =  8,
    ALL = 15
} MeshPointStateType;
typedef unsigned int MeshPointState;

	class MeshPoint{
		
	public:
		
		MeshPoint();
		
        MeshPoint(const Point3D &p);
		
        virtual ~MeshPoint();
		
        virtual void setPoint(const Point3D &p);
		
        //acces method:
        virtual Point3D &getPoint();
        virtual const Point3D &getPoint() const;

		virtual void addElement(unsigned int idx);
		
        virtual const list<unsigned int> &getElements() const;
		
        virtual void clearElements();
		
        virtual bool wasOutsideChecked() const;
		
		virtual void outsideChecked();
		
        //unconditionnal set
		virtual void setMaxDistance(double md);
		
        virtual double getMaxDistance() const;
		
        //set only if lower
        virtual void updateMaxDistance(double md);

        virtual void updateMaxDistanceByFactor(const double &per);
		
		virtual void setProjected();
		
        virtual bool wasProjected() const;
		
		//state methods
		virtual void setOutside();
		
		virtual void setInside();
		
		//returns true if node is inside any input mesh
        virtual bool isInside() const;
		
		//returns true if node is outside every input mesh
        virtual bool isOutside() const;
		
        virtual void featureProjected();
        
        virtual bool isFeature() const;
		
	protected:
		
		Point3D point;
		//this flag avoids to re-check if node is inside, 
		//which is an expensive task
//		bool outsidechecked, projected, feature;
		//inside is a flag to shrink elements to the surface.
		//it is a vector to know the state w.r.t. every input
		//geometry
//		bool inside;//, projected;
        MeshPointState state;
		list<unsigned int> elements;
		
		double maxdistance;
		
	};
	
	inline void MeshPoint::outsideChecked(){
        //outsidechecked = true;
        BITMASK_SET(state,OUTCHECKED);
	}
	
    inline bool MeshPoint::wasOutsideChecked() const{
        //return outsidechecked;
        return BITMASK_CHECK(state,OUTCHECKED);
	}
	
	inline void MeshPoint::setMaxDistance(double md){
		maxdistance = md;
	}
	
    inline double MeshPoint::getMaxDistance() const{
		return maxdistance;
	}
	
    inline void MeshPoint::updateMaxDistance(double md){
        if (maxdistance<md) {
            return;
        }
        maxdistance = md;
    }

    inline void MeshPoint::updateMaxDistanceByFactor(const double &per){
        maxdistance *= per;
    }

	inline void MeshPoint::setOutside(){
        //inside = false;
        BITMASK_CLEAR(state,INSIDE);
	}
	
	inline void MeshPoint::setInside(){
        //inside = true;
        BITMASK_SET(state,INSIDE);
    }
	
	//returns true if node is inside any input mesh
    inline bool MeshPoint::isInside() const {
        //return inside;
        return BITMASK_CHECK(state,INSIDE);
	}
	
	//returns true if node is outside every input mesh
    inline bool MeshPoint::isOutside() const{
        //return !inside;
        return !BITMASK_CHECK(state,INSIDE);
	}
	
	inline void MeshPoint::setProjected(){
        //projected = true;
        //inside = false;
        BITMASK_SET(state,PROJECTED);
        BITMASK_CLEAR(state,INSIDE);
	}
	
    inline bool MeshPoint::wasProjected() const {
        //return projected;
        return BITMASK_CHECK(state,PROJECTED);
	}
    
    inline void MeshPoint::featureProjected() {
        //feature = true;
        BITMASK_SET(state,PROJECTED);
    }
    
    inline bool MeshPoint::isFeature() const {
        //return feature;
        return BITMASK_CHECK(state,FEATURE);
    }
	
    inline Point3D &MeshPoint::getPoint(){
        return point;
    }
    inline const Point3D &MeshPoint::getPoint() const {
        return point;
    }

    inline void MeshPoint::setPoint(const Point3D &p){
		point = p;
	}
	
	inline void MeshPoint::addElement(unsigned int e){
		elements.push_back(e);
	}
	
	inline void MeshPoint::clearElements(){
		elements.clear();
	}
	
    inline const list<unsigned int> &MeshPoint::getElements() const {
		return elements;
	}
	
    std::ostream& operator<<(std::ostream& o, const MeshPoint &p);
	
}
#endif
