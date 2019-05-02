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
* @file SplitVisitorReductionOpenMP.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#ifndef SplitVisitorReductionOpenMP_h
#define SplitVisitorReductionOpenMP_h 1

#include <set>
#include <vector>

#include "../MeshPoint.h"
#include "../Point3D.h"
#include "../QuadEdge.h"

#include "../Visitors/Visitor.h"
#include <cassert>


using Clobscode::MeshPoint;
using Clobscode::QuadEdge;
using Clobscode::Point3D;
using std::vector;
using std::set;


namespace Clobscode
{
    class SplitVisitorReductionOpenMP : public Visitor {

    public:
        
        SplitVisitorReductionOpenMP();

        bool visit(Quadrant *o) override;

        void setPoints(const vector<MeshPoint> &points);
        
        void setNewPts(vector<Point3D> &new_pts);
        
        void setCurrentEdges(const set<QuadEdge> &edges);
        
        void setNewEdges(set<QuadEdge> &edges);
        
        void setNewEles(vector<vector<unsigned int> > &new_eles);
        
        void setClipping(vector<vector<Point3D> > &clipping);

        //OpenMP
        void setThreadNum(unsigned int tn) { this->tn = tn; }

    protected:
        
        //references
        const vector<MeshPoint> *points;
        vector<Point3D> *new_pts;
        const set<QuadEdge> *current_edges; //reading only
        set<QuadEdge> *new_edges; //Filled
        vector<vector<unsigned int> > *new_eles;
        vector<vector<Point3D> > *clipping;

        //Openmp
        unsigned int tn = 0;

        bool splitEdge(unsigned int idx1,
                       unsigned int idx2,
                       unsigned int &c_n_pts,
                       unsigned int &mid_idx);

    };
}



#endif //MESHER_ROI_SPLITVISITORReductionOpenMP_H
