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
* @file CustomSplitVisitor.cpp
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "CustomSplitVisitor.h"
#include "../../Quadrant.h"

namespace Clobscode
{
//vector<MeshPoint> *points;
//list<Point3D> *new_pts;
//set<QuadEdge> *edges;
//vector<vector<unsigned int> > *new_eles;
//vector<vector<Point3D> > *clipping;

    CustomSplitVisitor::CustomSplitVisitor()
        :points(NULL),new_pts(NULL),edges(NULL),new_eles(NULL),clipping(NULL)
    { }

    void CustomSplitVisitor::setPoints(const vector<MeshPoint> &points) {
        this->points = &points;
    }
    
    void CustomSplitVisitor::setNewPts(vector<Point3D> &new_pts) {
        this->new_pts = &new_pts;
    }

    void CustomSplitVisitor::setNewEdges(set<QuadEdge> &new_edges) {
        this->new_edges = &new_edges;
    }
    
    void CustomSplitVisitor::setEdges(set<QuadEdge> &edges) {
        this->edges = &edges;
    }
    
    void CustomSplitVisitor::setNewEles(vector<vector<unsigned int> > &new_eles) {
        this->new_eles = &new_eles;
    }
    
    void CustomSplitVisitor::setClipping(vector<vector<Point3D> > &clipping) {
        this->clipping = &clipping;
    }

    bool CustomSplitVisitor::visit(Quadrant *o)
    {
        //getting variables for modification
        //preferably by reference, to avoid unnecessary copying
        const vector<unsigned int> &pi = o->getPointIndex();

        new_eles->reserve(4);

        unsigned int n_pts = points->size() + new_pts->size();
        //the vector containing all nodes of this Quadrant (and sons)
        vector<unsigned int> all_pts(9,0);

        //save the four nodes of this square first
        for (unsigned int i=0; i < pi.size(); i++) {
            all_pts[i] = pi[i];
        }

        const Point3D &min = points->at(pi[0]).getPoint();
        const Point3D &max = points->at(pi[2]).getPoint();
        const Point3D avg = (max-min)/2 + min;

        //inserting node 4 between nodes 0 and 1
        if (splitEdge(all_pts[0],all_pts[1],n_pts,all_pts[4])) {
            //the coordinates of node 8 must be computed and added to
            //new_pts list of points
            new_pts->push_back(Point3D (avg[0],min[1],avg[2]));
        }
        //inserting node 5 between nodes 1 and 2
        if (splitEdge(all_pts[1],all_pts[2],n_pts,all_pts[5])) {
            //the coordinates of node 9 must be computed and added to
            //new_pts list of points
            new_pts->push_back(Point3D (max[0],avg[1],avg[2]));
        }
        //inserting node 6 between nodes 2 and 3
        if (splitEdge(all_pts[2],all_pts[3],n_pts,all_pts[6])) {
            //the coordinates of node 10 must be computed and added to
            //new_pts list of points
            new_pts->push_back(Point3D (avg[0],max[1],avg[2]));
        }
        //inserting node 7 between nodes 3 and 0
        if (splitEdge(all_pts[0],all_pts[3],n_pts,all_pts[7])) {
            //the coordinates of node 11 must be computed and added to
            //new_pts list of points
            new_pts->push_back(Point3D (min[0],avg[1],avg[2]));
        }

        //of course all the intern edges and mid point were never inserted
        //before, so this task is performed without asking
        new_pts->push_back(Point3D (avg[0],avg[1],avg[2]));
        all_pts[8] = n_pts;

//        QuadEdge intern_edge1 (all_pts[4],all_pts[6]);
//        intern_edge1.updateMidPoint(all_pts[8]);
//        QuadEdge intern_edge2 (all_pts[5],all_pts[7]);
//        intern_edge2.updateMidPoint(all_pts[8]);
//        edges->insert(intern_edge1);
//        edges->insert(intern_edge2);

        new_edges->emplace(all_pts[4],all_pts[6],all_pts[8]);
        new_edges->emplace(all_pts[5],all_pts[7],all_pts[8]);

        new_edges->emplace(all_pts[4],all_pts[8]);
        new_edges->emplace(all_pts[6],all_pts[8]);
        new_edges->emplace(all_pts[5],all_pts[8]);
        new_edges->emplace(all_pts[7],all_pts[8]);
        
        //now that all edges were inserted, the elements can be easily built
        vector<unsigned int> son_element (4,0);
        son_element[0]=all_pts[0];
        son_element[1]=all_pts[4];
        son_element[2]=all_pts[8];
        son_element[3]=all_pts[7];

        new_eles->push_back(son_element);

        son_element[0]=all_pts[4];
        son_element[1]=all_pts[1];
        son_element[2]=all_pts[5];
        son_element[3]=all_pts[8];

        new_eles->push_back(son_element);

        son_element[0]=all_pts[8];
        son_element[1]=all_pts[5];
        son_element[2]=all_pts[2];
        son_element[3]=all_pts[6];

        new_eles->push_back(son_element);

        son_element[0]=all_pts[7];
        son_element[1]=all_pts[8];
        son_element[2]=all_pts[6];
        son_element[3]=all_pts[3];

        new_eles->push_back(son_element);

        //extreme nodes of each son to be used by clipping
        //method.
        clipping->reserve(4);
        vector<Point3D> extreme_nodes(2, Point3D ());

        //bottom/left son is defined by nodes 0 and 8
        extreme_nodes[0] = Point3D (min.X(),min.Y(),avg.Z());
        extreme_nodes[1] = Point3D (avg.X(),avg.Y(),avg.Z());
        clipping->push_back(extreme_nodes);

        //bottom/right son is defined by nodes 4 and 5
        extreme_nodes[0] = Point3D (avg.X(),min.Y(),avg.Z());
        extreme_nodes[1] = Point3D (max.X(),avg.Y(),avg.Z());
        clipping->push_back(extreme_nodes);

        //top/right son is defined by nodes 8 and 2
        extreme_nodes[0] = Point3D (avg.X(),avg.Y(),avg.Z());
        extreme_nodes[1] = Point3D (max.X(),max.Y(),avg.Z());
        clipping->push_back(extreme_nodes);

        //bottom/left son is defined by nodes 7 and 6
        extreme_nodes[0] = Point3D (min.X(),avg.Y(),avg.Z());
        extreme_nodes[1] = Point3D (avg.X(),max.Y(),avg.Z());
        clipping->push_back(extreme_nodes);

        return true;
    }

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------


    bool CustomSplitVisitor::splitEdge(unsigned int idx1, unsigned int idx2,
                                 unsigned int &c_n_pts, unsigned int &mid_idx){
        
        QuadEdge this_edge (idx1,idx2);

        auto found = new_edges->find(this_edge);

        if (found != new_edges->end()) {

            if ((*found)[2]!=0) {
                //if the edge was already split, then save its mid_point and
                //return false (the current process didn't split the edge)
                mid_idx = (*found)[2];
                return false;
            }

            this_edge.updateMidPoint(c_n_pts++);
            found = new_edges->erase(found);
            //PAUL : todo change ? peut poser pb
            new_edges->insert(found, this_edge); //using found as hint for insertion


            // bulding and inserting, with hint if possible
            new_edges->emplace_hint(found, this_edge[0],this_edge[2]);
            new_edges->emplace_hint(found, this_edge[2],this_edge[1]);

            mid_idx = this_edge[2];

        } else  {
            found = edges->find(this_edge);

            if (found != edges->end() && (*found)[2]!=0) {
                //if the edge was already split, then save its mid_point and
                //return false (the current process didn't split the edge)
                mid_idx = (*found)[2];
                return false;
            }

            this_edge.updateMidPoint(c_n_pts++);
            //found = new_edges->erase(found);
            //PAUL : todo change ? peut poser pb
            new_edges->insert(this_edge); //using found as hint for insertion


            // bulding and inserting, with hint if possible
            new_edges->emplace(this_edge[0],this_edge[2]);
            new_edges->emplace(this_edge[2],this_edge[1]);

            mid_idx = this_edge[2];
        }





        //this edge is about to be split. Note that no edge can have point index
        //0 as its mid_point. For this reason, we know that the edge was not
        //split by other Quadrant before. The current edge must be replaced in the
        //set by the same one plus a mid_point equal to c_n_pts (current number
        //of points). The coordinates of this new point will be inserted by the
        //split method above. The splitEdge method will only compute the index
        //for this new point and will insert the two new edges (derived from the
        //splitting process of the current edge). Note that c_n_pts will be
        //increased for next splitting process of another edge.

        // std::set doesn't permit to update an element
        // we should erase and reinsert...
        /*this_edge.updateMidPoint(c_n_pts++);
        found = new_edges->erase(found);
        //PAUL : todo change ? peut poser pb
        new_edges->insert(this_edge); //using found as hint for insertion


        // bulding and inserting, with hint if possible
        new_edges->emplace(this_edge[0],this_edge[2]);
        new_edges->emplace(this_edge[2],this_edge[1]);

        mid_idx = this_edge[2];*/
        return true;
    }
}


