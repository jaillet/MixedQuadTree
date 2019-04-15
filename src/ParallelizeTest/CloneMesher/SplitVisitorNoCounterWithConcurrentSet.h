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
* @file SplitVisitorNoCounterWithConcurrentSet.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#ifndef SplitVisitorNoCounterWithConcurrentSet_h
#define SplitVisitorNoCounterWithConcurrentSet_h 1

#include <list>
#include <set>
#include <vector>
#include <tbb/concurrent_unordered_set.h>

#include "../../Visitors/Visitor.h"
#include "../../MeshPoint.h"
#include "../../QuadEdge.h"


using Clobscode::MeshPoint;
using Clobscode::QuadEdge;
using Clobscode::Point3D;
using std::vector;
using std::list;
using std::set;


namespace Clobscode
{
    class SplitVisitorNoCounterWithConcurrentSet : public Visitor {

    public:
        
        SplitVisitorNoCounterWithConcurrentSet();

        bool visit(Quadrant *o) override;

        void setPoints(const vector<MeshPoint> &points);
        
        void setNewPts(vector<Point3D> &new_pts);
        
        void setEdges(tbb::concurrent_unordered_set<QuadEdge, std::hash<QuadEdge>> &edges);
        
        void setNewEles(vector<vector<unsigned int> > &new_eles);
        
        void setClipping(vector<vector<Point3D> > &clipping);

    protected:
        
        //references
        const vector<MeshPoint> *points;
        vector<Point3D> *new_pts;
        tbb::concurrent_unordered_set<QuadEdge, std::hash<QuadEdge>> *edges;
        vector<vector<unsigned int> > *new_eles;
        vector<vector<Point3D> > *clipping;

        bool splitEdge(unsigned int idx1,
                       unsigned int idx2,
                       unsigned int &c_n_pts,
                       unsigned int &mid_idx);

    };
}



#endif //MESHER_ROI_SPLITVISITORNoCounterWithConcurrentSet_H
