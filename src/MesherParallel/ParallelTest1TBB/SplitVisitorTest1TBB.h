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
* @file SplitVisitorTest1.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#ifndef SplitVisitorTest1_h
#define SplitVisitorTest1_h 1

#include <list>
#include <set>
#include <vector>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_vector.h>
#include <mutex>


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
    class SplitVisitorTest1TBB : public Visitor {

    public:
        
        SplitVisitorTest1TBB();

        bool visit(Quadrant *o) override;

        void setPoints(const vector<MeshPoint> &points);
        
        void setNewPts(tbb::concurrent_vector<Point3D> &new_pts);
        
        void setEdges(tbb::concurrent_unordered_set<QuadEdge, std::hash<QuadEdge>> &edges);
        
        void setNewEles(vector<vector<unsigned int> > &new_eles);
        
        void setClipping(vector<vector<Point3D> > &clipping);

        void setCounterPoints(unsigned int counter);

    protected:
        
        //references
        const vector<MeshPoint> *points;
        tbb::concurrent_vector<Point3D> *new_pts;
        tbb::concurrent_unordered_set<QuadEdge, std::hash<QuadEdge>> *edges;
        vector<vector<unsigned int> > *new_eles;
        vector<vector<Point3D> > *clipping;
        unsigned int counter_points;

        bool splitEdge(unsigned int idx1,
                       unsigned int idx2,
                       unsigned int &mid_idx,
                       double idx5,
                       double idx6,
                       double idx7);

    };
}



#endif //MESHER_ROI_SplitVisitorTest1_H
