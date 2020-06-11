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
* @file FEMesh.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#ifndef FEMesh_h
#define FEMesh_h 1

#include <vector>
#include <array>
#include "Point3D.h"

using std::vector;

namespace Clobscode
{
	class FEMesh{
		
	public:
		
		FEMesh(){};
		
        virtual void setPoints(const vector<Point3D> &pts);
		
        virtual void setElements(const vector<vector<unsigned int> > &els);
		
        virtual void setOutsideNodes(const list<unsigned int> &outpts);
		
        virtual vector<Point3D> &getPoints();
        virtual const vector<Point3D> &getPoints() const ;

        virtual vector<vector<unsigned int> > &getElements();
        virtual const vector<vector<unsigned int> > &getElements() const;

        virtual list<unsigned int> &getOutsideNodes();
        virtual const list<unsigned int> &getOutsideNodes() const;
        
        //For debugging:
        virtual void setColoredCells(const vector<unsigned int> &colored);
        
        virtual const vector<unsigned int> &getColoredCells();

        //For decoration
        virtual const vector<unsigned short> &getRefLevels() const;
        virtual void setRefLevels(const vector<unsigned short> &rl);

        virtual const vector <double> &getMinAngles() const;
        virtual void setMinAngles(const vector<double> &ma);
        virtual const vector <double> &getMaxAngles() const;
        virtual void setMaxAngles(const vector<double> &ma);

        virtual const array <unsigned int,18> &getAnglesTriHistogram() const;
        virtual void setAnglesTriHistogram(const array<unsigned int,18> &ah);
        virtual const array <unsigned int,18> &getAnglesQuadHistogram() const;
        virtual void setAnglesQuadHistogram(const array<unsigned int,18> &ah);

        virtual const vector<unsigned short> &getSurfState() const;
        virtual void setSurfState(const vector<unsigned short> &surf);
        
        virtual const vector<unsigned short> &getDebugging() const;
        virtual void setDebugging(const vector<unsigned short> &deb);

	protected:
		
		vector<Point3D> points;
        vector<vector<unsigned int> > elements;
        vector <unsigned short> ref_levels, surf_state, deb_state;
        vector <double> min_angles, max_angles;
        array<unsigned int,18> angles_tri_histogram,  angles_quad_histogram;
        vector<unsigned int> color;
        list<unsigned int> outpts;
		
	};
    
    inline const vector<unsigned short> &FEMesh::getRefLevels() const {return ref_levels;}
    inline void FEMesh::setRefLevels(const vector<unsigned short> &rl) {ref_levels=rl;}

    inline const vector<double> &FEMesh::getMinAngles() const {return min_angles;}
    inline void FEMesh::setMinAngles(const vector<double> &ma) {min_angles=ma;}
    inline const vector<double> &FEMesh::getMaxAngles() const {return max_angles;}
    inline void FEMesh::setMaxAngles(const vector<double> &ma) {max_angles=ma;}

    inline const array<unsigned int,18> &FEMesh::getAnglesTriHistogram() const {return angles_tri_histogram;}
    inline void FEMesh::setAnglesTriHistogram(const array<unsigned int,18> &ah) {angles_tri_histogram=ah;}
    inline const array<unsigned int,18> &FEMesh::getAnglesQuadHistogram() const {return angles_quad_histogram;}
    inline void FEMesh::setAnglesQuadHistogram(const array<unsigned int,18> &ah) {angles_quad_histogram=ah;}

    inline const vector<unsigned short> &FEMesh::getSurfState() const {return surf_state;}
    inline void FEMesh::setSurfState(const vector<unsigned short> &surf) {surf_state=surf;}
    
    inline const vector<unsigned short> &FEMesh::getDebugging() const {return deb_state;}
    inline void FEMesh::setDebugging(const vector<unsigned short> &deb) {deb_state=deb;}

    inline void FEMesh::setColoredCells(const vector<unsigned int> &colored) {
        color = colored;
    }
    
    inline const vector<unsigned int> &FEMesh::getColoredCells() {
        return color;
    }
	
    inline void FEMesh::setOutsideNodes(const list<unsigned int> &outpts){
        this->outpts=outpts;
	}
	
    inline void FEMesh::setElements(const vector<vector<unsigned int> > &els){
        elements=els;
	}
	
    inline void FEMesh::setPoints(const vector<Point3D> &pts){
        points=pts;
	}
	
    inline vector<Point3D> &FEMesh::getPoints(){
        return points;
    }
    inline const vector<Point3D> &FEMesh::getPoints() const{
        return points;
    }

    inline vector<vector<unsigned int> > &FEMesh::getElements(){
        return elements;
    }
    inline const vector<vector<unsigned int> > &FEMesh::getElements() const {
        return elements;
    }

    inline list<unsigned int> &FEMesh::getOutsideNodes(){
        return outpts;
    }
    inline const list<unsigned int> &FEMesh::getOutsideNodes() const{
        return outpts;
    }
}
#endif
