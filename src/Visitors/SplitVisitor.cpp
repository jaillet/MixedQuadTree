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
* @file SplitVisitor.cpp
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "SplitVisitor.h"
#include "../Quadrant.h"

namespace Clobscode
{
//vector<MeshPoint> *points;
//list<Point3D> *new_pts;
//set<QuadEdge> *edges;
//vector<vector<unsigned int> > *new_eles;
//vector<vector<Point3D> > *clipping;

    SplitVisitor::SplitVisitor()
        :points(NULL),new_pts(NULL),MapEdges(NULL),new_eles(NULL),clipping(NULL)
    { idx = 0; }

    void SplitVisitor::setPoints(const vector<MeshPoint> &points) {
        this->points = &points;
    }
    
    void SplitVisitor::setNewPts(list<Point3D> &new_pts) {
        this->new_pts = &new_pts;
    }
    
    void SplitVisitor::setMapEdges(map<QuadEdge, EdgeInfo> &MapEdges) {
        this->MapEdges = &MapEdges;
    }
    
    void SplitVisitor::setNewEles(vector<vector<unsigned int> > &new_eles) {
        this->new_eles = &new_eles;
    }
    
    void SplitVisitor::setStartIndex(const unsigned int &sidx) {
        this->idx = sidx;
    }
    
    void SplitVisitor::setClipping(vector<vector<Point3D> > &clipping) {
        this->clipping = &clipping;
    }
    
    void SplitVisitor::setProcessedQuadVector(vector<Quadrant> &proQuadVec) {
        this->proQuadVec = &proQuadVec;
    }
    
    void SplitVisitor::setMapProcessed(map<unsigned int, unsigned int> &proQuadMap) {
        this->proQuadMap = &proQuadMap;
    }
    
    void SplitVisitor::setToBalanceList(list<pair<unsigned int,unsigned int> > &unBalanced) {
        this->unBalanced = &unBalanced;
    }

    bool SplitVisitor::visit(Quadrant *o)
    {
        //getting variables for modification
        //preferably by reference, to avoid unnecessary copying
        const vector<unsigned int> &pi = o->pointindex;

        new_eles->reserve(4);

        unsigned int n_pts = points->size() + new_pts->size();
        //the vector containing all nodes of this Quadrant (and sons)
        vector<unsigned int> all_pts(9,0);

        //save the four nodes of this square first
        for (unsigned int i=0; i< pi.size(); i++) {
            all_pts[i] = pi[i];
        }

        const Point3D &min = points->at(pi[0]).getPoint();
        const Point3D &max = points->at(pi[2]).getPoint();
        const Point3D avg = (max-min)/2 + min;
        
        unsigned int son_rl = o->getRefinementLevel() + 1;

        //inserting node 4 between nodes 0 and 1
        unsigned int nq_idx;
        if (splitEdge(all_pts[0],all_pts[1],idx,idx+1,1,n_pts,nq_idx,all_pts[4])) {
            
            //save the new mid_edge index
            all_pts[4] = n_pts;
            //create the Edges with parent neighbor
            MapEdges->emplace(QuadEdge (all_pts[0],n_pts,true), EdgeInfo (0,idx,nq_idx));
            MapEdges->emplace(QuadEdge (all_pts[1],n_pts,true), EdgeInfo (0,idx+1,nq_idx));
            
            //increase the index
            n_pts++;
            //the coordinates of node 4 must be computed and added to
            //new_pts list of points
            new_pts->push_back(Point3D (avg[0],min[1],avg[2]));
            //We have to check if neighbor Quad is unBalanced
            auto found = proQuadMap->find(nq_idx);
            if (found!=proQuadMap->end()) {
                unsigned int nrl = (proQuadVec->at(found->second)).getRefinementLevel();
                int diff = nrl - son_rl;
                
                if (diff>1 or diff<-1) {
                    unBalanced->push_back(*found);
                }
            }
        }
        
        //inserting node 5 between nodes 1 and 2
        if (splitEdge(all_pts[1],all_pts[2],idx+1,idx+2,1,n_pts,nq_idx,all_pts[5])) {
            //save the new mid_edge index
            all_pts[5] = n_pts;
            //create the Edges with parent neighbor
            MapEdges->emplace(QuadEdge (all_pts[1],n_pts,true), EdgeInfo (0,idx+1,nq_idx));
            MapEdges->emplace(QuadEdge (all_pts[2],n_pts,true), EdgeInfo (0,idx+2,nq_idx));
            
            //increase the index
            n_pts++;
            //the coordinates of node 5 must be computed and added to
            //new_pts list of points
            new_pts->push_back(Point3D (max[0],avg[1],avg[2]));
            //We have to check if neighbor Quad is unBalanced
            auto found = proQuadMap->find(nq_idx);
            if (found!=proQuadMap->end()) {
                unsigned int nrl = (proQuadVec->at(found->second)).getRefinementLevel();
                int diff = nrl - son_rl;
                
                if (diff>1 or diff<-1) {
                    unBalanced->push_back(*found);
                }
            }
        }
        
        //inserting node 6 between nodes 2 and 3
        if (splitEdge(all_pts[2],all_pts[3],idx+2,idx+3,2,n_pts,nq_idx,all_pts[6])) {
            //save the new mid_edge index
            all_pts[6] = n_pts;
            //create the Edges with parent neighbor
            MapEdges->emplace(QuadEdge (all_pts[2],n_pts,true), EdgeInfo (0,nq_idx,idx+2));
            MapEdges->emplace(QuadEdge (all_pts[3],n_pts,true), EdgeInfo (0,nq_idx,idx+3));
            
            //increase the index
            n_pts++;
            //the coordinates of node 6 must be computed and added to
            //new_pts list of points
            new_pts->push_back(Point3D (avg[0],max[1],avg[2]));
            //We have to check if neighbor Quad is unBalanced
            auto found = proQuadMap->find(nq_idx);
            if (found!=proQuadMap->end()) {
                unsigned int nrl = (proQuadVec->at(found->second)).getRefinementLevel();
                int diff = nrl - son_rl;
                
                if (diff>1 or diff<-1) {
                    unBalanced->push_back(*found);
                }
            }
        }
        
        //inserting node 7 between nodes 3 and 0
        if (splitEdge(all_pts[3],all_pts[0],idx+3,idx,2,n_pts,nq_idx,all_pts[7])) {
            //save the new mid_edge index
            all_pts[7] = n_pts;
            //create the Edges with parent neighbor
            MapEdges->emplace(QuadEdge (all_pts[3],n_pts,true), EdgeInfo (0,nq_idx,idx+3));
            MapEdges->emplace(QuadEdge (all_pts[0],n_pts,true), EdgeInfo (0,nq_idx,idx));
            
            //increase the index
            n_pts++;
            //the coordinates of node 7 must be computed and added to
            //new_pts list of points
            new_pts->push_back(Point3D (min[0],avg[1],avg[2]));
            //We have to check if neighbor Quad is unBalanced
            auto found = proQuadMap->find(nq_idx);
            if (found!=proQuadMap->end()) {
                unsigned int nrl = (proQuadVec->at(found->second)).getRefinementLevel();
                int diff = nrl - son_rl;
                
                if (diff>1 or diff<-1) {
                    unBalanced->push_back(*found);
                }
            }
        }

        //of course all the intern edges and mid point were never inserted
        //before, so this task is performed without asking
        new_pts->push_back(Point3D (avg[0],avg[1],avg[2]));
        all_pts[8] = n_pts;

        //NOTE: edge (4,6) and (5,7) with mid_point 8 (the center of the
        //Quad are not inserted because those edges ((4,6) and (5,7)) don't
        //belong to any element. These edges would only increase the use of
        //memory without adding necessary information.
        //edges->emplace(all_pts[4],all_pts[6],all_pts[8]);
        //edges->emplace(all_pts[5],all_pts[7],all_pts[8]);
        /*QuadEdge qeY(all_pts[4],all_pts[6]);
        MapEdges->emplace(qeY,all_pts[8]);
        
        QuadEdge qeX(all_pts[5],all_pts[7]);
        MapEdges->emplace(qeX,all_pts[8]);*/
        
        MapEdges->emplace(QuadEdge (all_pts[4],all_pts[8],true),EdgeInfo (0,idx,idx+1));
        MapEdges->emplace(QuadEdge (all_pts[6],all_pts[8],true),EdgeInfo (0,idx+3,idx+2));
        MapEdges->emplace(QuadEdge (all_pts[5],all_pts[8],true),EdgeInfo (0,idx+2,idx+1));
        MapEdges->emplace(QuadEdge (all_pts[7],all_pts[8],true),EdgeInfo (0,idx+3,idx));
        
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
     bool SplitVisitor::splitEdge(const unsigned int &idx1, const unsigned int &idx2,
                                  const unsigned int &q1, const unsigned int &q2,
                                  const unsigned int &pos, const unsigned int &nidx,
                                  unsigned int &nquad, unsigned int &mid_idx) {
         
         //If the Edge was split, update Quad info. Else: return the Quad index
         //of the other side of the edge
         auto help = MapEdges->find(QuadEdge (idx1,idx2));
         mid_idx = (help->second)[0];
         
         if (mid_idx!=0) {
             
             //update Quad neighbor info in EdgeInfo
             auto e1 = MapEdges->find(QuadEdge (idx1,mid_idx));
             e1->second[pos] = q1;
             
             auto e2 = MapEdges->find(QuadEdge (idx2,mid_idx));
             e2->second[pos] = q2;
             
             //Important note. It is possible that edges (idxX,mid_idx) have been already
             //split. In order to mantain updated info, it is necessary to check for
             //evetual edge sons and make them congruent.
             
             //auto sub1 = MapEdges->find(QuadEdge (idx1,mid_idx));
             //unsigned int sub_mid = (sub1->second)[0];
             unsigned int sub_mid = (e1->second)[0];
             if (sub_mid!=0) {
                 (MapEdges->find(QuadEdge (idx1,sub_mid))->second)[pos] = q1;
                 (MapEdges->find(QuadEdge (mid_idx,sub_mid))->second)[pos] = q1;
             }
             
             //auto sub2 = MapEdges->find(QuadEdge (idx2,mid_idx));
             //sub_mid = (sub2->second)[0];
             sub_mid = (e2->second)[0];
             if (sub_mid!=0) {
                 (MapEdges->find(QuadEdge (idx2,sub_mid))->second)[pos] = q2;
                 (MapEdges->find(QuadEdge (mid_idx,sub_mid))->second)[pos] = q2;
             }
             
             return false;
         }
         
         //Save in the current edge the new mid_index
         (help->second)[0] = nidx;
         
         //return neighbor quad index
         nquad = (help->second)[3-pos];
         return true;
     }

}


