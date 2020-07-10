/*
 <Mix-mesher: region type. This program generates a mixed-elements 2D mesh>
 
 Copyright (C) <2013,2020>  <Claudio Lobos> All rights reserved.
 
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
 * @file Mesher.cpp
 * @author Claudio Lobos, Fabrice Jaillet
 * @version 0.1
 * @brief
 **/

#include "Mesher.h"
#include <math.h>
#include <iomanip>

namespace Clobscode
{
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    Mesher::Mesher() {}
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    Mesher::~Mesher() {}
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    //Create a grid mesh regarding the Bounding Box of input mesh.
    //This will produce several cubes as roots of an octree structure.
    //Then split each initial element 4^rl times (where rl stands
    //for Refinement Level).
    std::shared_ptr<FEMesh> Mesher::refineMesh(Polyline &input, const unsigned short &rl,
                                               const string &name, list<unsigned int> &roctli,
                                               list<RefinementRegion *> &all_reg,
                                               GeometricTransform &gt,
                                               const unsigned short &minrl,
                                               const unsigned short &omaxrl,
                                               const bool &debugging, bool decoration){
        
        //Note: rotation are not enabled when refining an already produced mesh.
        bool rotated = !gt.Default();
        if(rotated) {
            /*cout << "rotating input surface mesh\n";
             cout << gt.getCentroid() << "\n";
             cout << gt.getXAxis() << " " << gt.getYAxis();
             cout << " " << gt.getZAxis() << "\n";*/
            gt.rotatePolyline(input);
        }
        
#if (VTKOUT==true)         //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> grid_octree=make_shared<FEMesh>();
            saveOutputMesh(grid_octree,points,Quadrants,debugging);
            string tmp_name = name + "_grid";
            Services::WriteVTK(tmp_name,grid_octree);
        }
#endif

        //split Quadrants until the refinement level (rl) is achieved.
        //The output will be a one-irregular mesh.
        splitQuadrants(rl,input,roctli,all_reg,name,minrl,omaxrl,debugging);
        
        //Save the Octant mesh for further refinement.
        Services::WriteQuadtreeMesh(name,points,Quadrants,MapEdges,gt);
        //Some Quads will be then removed due to proximity with the surface.
        //However we must preserve them if the oct mesh to avoid congruency
        //problems. For this reason we will keep track of removed quads
        //so we can easily link elements to quad index when reading an oct
        //mesh file.
        map<unsigned int, bool> removedquads;
        list<unsigned int> quadmeshidx;
        for (auto q: Quadrants) {
            removedquads[q.getIndex()] = true;
            quadmeshidx.push_back(q.getIndex());
        }
        
        //link element and node info for code optimization, also
        //detect Quadrants with features.
        detectFeatureQuadrants(input);
        linkElementsToNodes();
        detectInsideNodes(input);
        computeNodeMaxDist();
        
#if (VTKOUT==true)         //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> pure_octree=make_shared<FEMesh>();
            saveOutputMesh(pure_octree,points,Quadrants,debugging);
            string tmp_name = name + "_quads";
            Services::WriteVTK(tmp_name,pure_octree);
        }
#endif

        projectCloseToBoundaryNodes(input);

#if (VTKOUT==true) //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> closeto_octree=make_shared<FEMesh>();
            saveOutputMesh(closeto_octree,points,Quadrants);
            string tmp_name = name + "_closeto";
            Services::WriteVTK(tmp_name,closeto_octree);
        }
#endif
        
        removeOnSurfaceSafe(input);
        
        //update element and node info.
        linkElementsToNodes();
        
#if (VTKOUT==true)        //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> pure_octree=make_shared<FEMesh>();
            saveOutputMesh(pure_octree,points,Quadrants);
            string tmp_name = name + "_remSur";
            Services::WriteVTK(tmp_name,pure_octree);
        }
#endif

        //shrink outside nodes to the input domain boundary
        shrinkToBoundary(input);
        
#if (VTKOUT==true)        //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> shrink_octree=make_shared<FEMesh>();
            saveOutputMesh(shrink_octree,points,Quadrants);
            string tmp_name = name + "_shrink";
            Services::WriteVTK(tmp_name,shrink_octree);
        }
#endif

        //apply the surface Patterns
        applySurfacePatterns(input);
        //removeOnSurface(input);
        
        if (rotated) {
            // rotate the mesh
            for (unsigned int i=0; i<points.size(); i++) {
                gt.applyInverse(points[i].getPoint());
            }
            // rotate back the polyline as well
            //FJA no need for Refinement ??
            //for (unsigned int i=0; i<input.getPoints().size(); i++) {
            //    gt.applyInverse(input.getPoints()[i]);
            //}
        }
        
        //the almighty output mesh
        std::shared_ptr<FEMesh> mesh = std::make_shared<FEMesh>();
        
        //save the data of the mesh in its final state
        saveOutputMesh(mesh,decoration);
        
        //Write element-quad info the file
        Services::addOctElemntInfo(name,Quadrants,removedquads,quadmeshidx);
        
        return mesh;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    //Create a grid mesh regarding the Bounding Box of input mesh.
    //This will produce several cubes as roots of an octree structure.
    //Then split each initial element 4^rl times (where rl stands
    //for Refinement Level).
    std::shared_ptr<FEMesh> Mesher::generateMesh(Polyline &input, const unsigned short &rl,
                                                 const string &name,
                                                 list<RefinementRegion *> &all_reg,
                                                 const bool &debugging, bool decoration){
        
        //ATTENTION: geometric transform causes invalid input rotation when the
        //input is a cube.
        GeometricTransform gt;
        
        //rotate: This method is written below and its mostly commented because
        //it causes conflicts when the input is a cube. Must be checked.
        bool rotated = rotateGridMesh(input, all_reg, gt);
        
        //generate root Quadrants
        generateGridMesh(input);
        
#if (VTKOUT==true)         //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> grid_octree = make_shared<FEMesh>();
            saveOutputMesh(grid_octree,points,Quadrants,debugging);
            string tmp_name = name + "_grid";
            Services::WriteVTK(tmp_name,grid_octree);
        }
#endif
        
        //split Quadrants until the refinement level (rl) is achieved.
        //the last 0 correspond to min RL in the mesh. As this mesh
        //is starting from scratch, this value is set to 0. When refining
        //an existing mesh, this value may change.
        //The output will be a one-irregular mesh.
        generateQuadtreeMesh(rl,input,all_reg,name,0,debugging,Quadrants.size());
        
        Services::WriteQuadtreeMesh(name,points,Quadrants,MapEdges,gt);
        //Some Quads will be then removed due to proximity with the surface.
        //However we must preserve them if the oct mesh to avoid congruency
        //problems. For this reason we will keep track of removed quads
        //so we can easily link elements to quad index when reading an oct
        //mesh file.
        map<unsigned int, bool> removedquads;
        list<unsigned int> quadmeshidx;
        for (auto q: Quadrants) {
            removedquads[q.getIndex()] = true;
            quadmeshidx.push_back(q.getIndex());
        }
        
        
        //link element and node info for code optimization, also
        //detect Quadrants with features.
        detectFeatureQuadrants(input);
        linkElementsToNodes();
        detectInsideNodes(input);
        computeNodeMaxDist();
        
#if (VTKOUT==true)         //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> pure_octree = make_shared<FEMesh>();
            saveOutputMesh(pure_octree,points,Quadrants);
            string tmp_name = name + "_quads";
            Services::WriteVTK(tmp_name,pure_octree);
        }
#endif

        projectCloseToBoundaryNodes(input);
        
#if (VTKOUT==true)         //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> closeto_octree = make_shared<FEMesh>();
            saveOutputMesh(closeto_octree,points,Quadrants);
            string tmp_name = name + "_closeto";
            Services::WriteVTK(tmp_name,closeto_octree);
        }
#endif

        removeOnSurfaceSafe(input);
        
        //update element and node info.
        linkElementsToNodes();
        
#if (VTKOUT==true)        //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> pure_octree = make_shared<FEMesh>();
            saveOutputMesh(pure_octree,points,Quadrants);
            string tmp_name = name + "_remSur";
            Services::WriteVTK(tmp_name,pure_octree);
        }
#endif

        //shrink outside nodes to the input domain boundary
        shrinkToBoundary(input);
        
#if (VTKOUT==true)        //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> shrink_octree = make_shared<FEMesh>();
            saveOutputMesh(shrink_octree,points,Quadrants);
            string tmp_name = name + "_shrink";
            Services::WriteVTK(tmp_name,shrink_octree);
        }
#endif

        //apply the surface Patterns
        applySurfacePatterns(input);
        
        if (rotated) {
            // rotate the mesh
            for (unsigned int i=0; i<points.size(); i++) {
                gt.applyInverse(points[i].getPoint());
            }
            // rotate back the polyline as well
            for (unsigned int i=0; i<input.getPoints().size(); i++) {
                gt.applyInverse(input.getPoints()[i]);
            }
        }
        
        //the almighty output mesh
        std::shared_ptr<FEMesh> mesh = std::make_shared<FEMesh>();
        
        //save the data of the mesh in its final state
        saveOutputMesh(mesh,decoration);
        
        //Write element-quad info the file
        Services::addOctElemntInfo(name,Quadrants,removedquads,quadmeshidx);
        
        return mesh;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    bool Mesher::rotateGridMesh(Polyline &input, list<RefinementRegion *> &all_reg,
                                GeometricTransform &gt){
        
        list<RefinementRegion *>::const_iterator it, rrot;
        bool inputHasbeenRotated = false;
        
        for (it = all_reg.begin(); it!=all_reg.end(); it++) {
            //in case of input roi
            if((*it)->needsInputRotation()){
                if(!inputHasbeenRotated){
                    gt = (*it)->rotateWithinYou(input);
                    inputHasbeenRotated = true;
                    rrot=it;
                    break;
                }
            }
        }
        
        if (inputHasbeenRotated) {
            for (it = all_reg.begin(); it!=all_reg.end(); it++) {
                if (it!=rrot) {
                    if ((*it)->needsLocalRotation()) {
                        (*it)->rotate(gt);
                    }
                }
            }
        }
        
        return inputHasbeenRotated;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    void Mesher::generateGridMesh(Polyline &input){
        //vectors with each coordinate per axis
        vector<double> all_x, all_y;
        vector<vector<unsigned int> > elements;
        
        auto start_time = chrono::high_resolution_clock::now();
        
        GridMesher gm;
        gm.generatePoints(input.getBounds(),all_x,all_y);
        gm.generateMesh(all_x,all_y,points,elements);
        
        Quadrants.reserve(elements.size());
        
        //create the root Quadrants
        for (unsigned int i=0; i<elements.size(); i++) {
            Quadrant o (elements[i], 0, i);
            //Only when the Quadrant intersects the input
            //add it to the list of current Quadrants. As
            //This is the first time Quadrants are checked
            //for intersections they must be made w.r.t.
            //all input edges.
            IntersectionsVisitor iv(false);
            //if (o.checkIntersections(input,points)) {
            iv.setPolyline(input);
            iv.setPoints(points);
            if (o.accept(&iv)) {
                EdgeVisitor::insertEdges(&o,MapEdges);
                Quadrants.push_back(o);
            }
        }
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * generateGridMesh in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    void Mesher::setInitialState(vector<MeshPoint> &epts, vector<Quadrant> &eocts,
                                 map<QuadEdge, EdgeInfo> &medgs) {
        Quadrants.assign(make_move_iterator(eocts.begin()),
                         make_move_iterator(eocts.end()));
        
        //Erase previous quadrants to save memory
        eocts.erase(eocts.begin(),eocts.end());
        
        points.assign(make_move_iterator(epts.begin()),make_move_iterator(epts.end()));
        epts.erase(epts.begin(),epts.end());
        
        MapEdges = medgs;
        //recompute Octant indexes and update edges with those indexes
        for (unsigned int i=0;i<Quadrants.size();i++) {
            vector<unsigned int> qpts = Quadrants[i].getPointIndex();
            for (unsigned int j=0; j<4; j++) {
                
                auto qe = MapEdges.find(QuadEdge (qpts[j],qpts[(j+1)%4]));
                
                if (qe==MapEdges.end()) {
                    cout << "  Error at Mesher::setInitialState edge not found\n";
                }
                
                unsigned int mide = (qe->second)[0];
                if (mide!=0) {
                    auto qes1 = MapEdges.find(QuadEdge (qpts[j],mide));
                    auto qes2 = MapEdges.find(QuadEdge (qpts[(j+1)%4],mide));
                    if (j==0 || j==1) {
                        (qes1->second)[1] = i;
                        (qes2->second)[1] = i;
                    }
                    else {
                        (qes1->second)[2] = i;
                        (qes2->second)[2] = i;
                    }
                }
                else {
                    if (j==0 || j==1) {
                        (qe->second)[1] = i;
                    }
                    else {
                        (qe->second)[2] = i;
                    }
                }
            }
        }
        cout << "  Mesh initialized\n";

//        auto end_time = chrono::high_resolution_clock::now();
//        cout << "    * generateGridMesh in "
//        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
//        cout << " ms"<< endl;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    void Mesher::splitQuadrants(const unsigned short &rl, Polyline &input,
                                list<unsigned int> &roctli,
                                list<RefinementRegion *> &all_reg, const string &name,
                                const unsigned short &minrl,
                                const unsigned short &maxrl,
                                const bool &debugging){
        
        //The list of candidate quads to refine and the tmp version of
        //adding those how are still candidates for the next iteration.
        list<Quadrant> candidates, new_candidates, clean_processed, refine_tmp;
        
        //The Quads that don't need further refinement.
        vector<Quadrant> processed;
        
        //A map that connect the index of processed Quads with their position in the
        //vector of processed Quads.
        map<unsigned int, unsigned int> idx_pos_map;
        
        //list of the points added at this refinement iteration:
        list<Point3D> new_pts;
        
        //list<Quadrant>::iterator iter;
        
        //A set containing the index of Quads to be refined to
        //maintain balancing
        list<pair<unsigned int,unsigned int> > toBalance;
        
        //initialising the vector and the map
        candidates.assign(make_move_iterator(Quadrants.begin()),
                          make_move_iterator(Quadrants.end()));
        
        //The starting point for assignaiting indexes
        unsigned int new_q_idx = candidates.size();
        unsigned int future_q_idx = new_q_idx + 4*(unsigned int)roctli.size();
        
        //Erase previous quadrants to save memory
        Quadrants.erase(Quadrants.begin(),Quadrants.end());
        
        //create visitors and give them variables
        SplitVisitor sv;
        sv.setPoints(points);
        sv.setMapEdges(MapEdges);
        sv.setNewPts(new_pts);
        sv.setProcessedQuadVector(processed);
        sv.setMapProcessed(idx_pos_map);
        sv.setToBalanceList(toBalance);
        
        auto start_refine_quad_time = chrono::high_resolution_clock::now();
        
        bool listref = false, localmax = false;
        if (!roctli.empty()) {
            listref = true;
        }
        
        
        //----------------------------------------------------------
        //refine once each Quadrant in the list
        //----------------------------------------------------------
        
        //note: roctli is an index list sorted when readed.
        //First we need to build the map with actual index
        //of quads in the processed quad vector. The quads
        //needing refinement will go into another list.
        if (!roctli.empty()) {
            
            list<unsigned int>::iterator octidx = roctli.begin();
            
            unsigned int qua_pos = 0;
            while (!candidates.empty()) {
                
                Quadrant quad = *(candidates.begin());
                candidates.pop_front();
                
                if (octidx == roctli.end() || qua_pos!=*octidx) {
                    idx_pos_map[quad.getIndex()] = processed.size();
                    processed.push_back(quad);
                }
                else {
                    //we advance to next quadrant in the list for the next iteration.
                    octidx++;
                    refine_tmp.push_back(quad);
                }
                qua_pos++;
            }
            
            for (auto quad: refine_tmp) {
                
                //start refinement process for current quadrant.
                list<unsigned int> &inter_edges = quad.getIntersectedEdges();
                unsigned short qrl = quad.getRefinementLevel();
                
                if (qrl==maxrl) {
                    localmax = true;
                }
                
                vector<vector<Point3D> > clipping_coords;
                sv.setClipping(clipping_coords);
                
                vector<vector<unsigned int> > split_elements;
                sv.setNewEles(split_elements);
                sv.setStartIndex(new_q_idx);
                
                quad.accept(&sv);
                
                if (inter_edges.empty()) {
                    for (unsigned int j=0; j<split_elements.size(); j++) {
                        Quadrant o (split_elements[j], qrl+1, new_q_idx++);
                        new_candidates.push_back(o);
                    }
                }
                else {
                    for (unsigned int j=0; j<split_elements.size(); j++) {
                        Quadrant o (split_elements[j],qrl+1,new_q_idx++);
                        //the new points are inserted in bash at the end of this
                        //iteration. For this reason, the coordinates must be passed
                        //"manually" at this point (clipping_coords).
                        
                        //select_faces = true
                        IntersectionsVisitor iv(true);
                        //if (o.checkIntersections(input,inter_edges,clipping_coords[j]))
                        iv.setPolyline(input);
                        iv.setEdges(inter_edges);
                        iv.setCoords(clipping_coords[j]);
                        
                        if (o.accept(&iv)) {
                            new_candidates.push_back(o);
                        }
                        else {
                            //The element doesn't intersect any input edge.
                            //It must be checked if it's inside or outside.
                            //Only in the first case add it to new_Quadrants.
                            //Test this with parent Quadrant faces only.
                            if (isItIn(input,inter_edges,clipping_coords[j])) {
                                new_candidates.push_back(o);
                            }
                            else {
                                //we must update neighbor information at the edges
                                EdgeVisitor::removeEdges(&o, MapEdges);
                            }
                        }
                    }
                }
            }
        
            //Erase the list to refine
            refine_tmp.erase(refine_tmp.begin(),refine_tmp.end());
            
            //cout << "To balance list has " << toBalance.size() << " quads\n";
        
            while (!toBalance.empty()) {
                
                list<pair<unsigned int, unsigned int> > tmp_toBalance;
                std::swap(toBalance,tmp_toBalance);
                tmp_toBalance.sort();
                tmp_toBalance.unique();
                
                while (!tmp_toBalance.empty()) {
                    unsigned int key = tmp_toBalance.begin()->first;
                    unsigned int val = tmp_toBalance.begin()->second;
                    
                    Quadrant quad = processed[val];
                    tmp_toBalance.pop_front();
                    list<unsigned int> &inter_edges = quad.getIntersectedEdges();
                    unsigned short qrl = quad.getRefinementLevel();
                    
                    vector<vector<Point3D> > clipping_coords;
                    sv.setClipping(clipping_coords);
                    
                    vector<vector<unsigned int> > split_elements;
                    sv.setNewEles(split_elements);
                    sv.setStartIndex(new_q_idx);
                    
                    quad.accept(&sv);
                    
                    if (inter_edges.empty()) {
                        for (unsigned int j=0; j<split_elements.size(); j++) {
                            Quadrant o (split_elements[j], qrl+1, new_q_idx++);
                            idx_pos_map[o.getIndex()] = processed.size();
                            processed.push_back(o);
                        }
                    }
                    else {
                        for (unsigned int j=0; j<split_elements.size(); j++) {
                            Quadrant o (split_elements[j],qrl+1,new_q_idx++);
                            //the new points are inserted in bash at the end of this
                            //iteration. For this reason, the coordinates must be passed
                            //"manually" at this point (clipping_coords).
                            
                            //select_faces = true
                            IntersectionsVisitor iv(true);
                            //if (o.checkIntersections(input,inter_edges,clipping_coords[j]))
                            iv.setPolyline(input);
                            iv.setEdges(inter_edges);
                            iv.setCoords(clipping_coords[j]);
                            
                            if (o.accept(&iv)) {
                                idx_pos_map[o.getIndex()] = processed.size();
                                processed.push_back(o);
                                
                            }
                            else {
                                //The element doesn't intersect any input edge.
                                //It must be checked if it's inside or outside.
                                //Only in the first case add it to new_Quadrants.
                                //Test this with parent Quadrant faces only.
                                if (isItIn(input,inter_edges,clipping_coords[j])) {
                                    idx_pos_map[o.getIndex()] = processed.size();
                                    processed.push_back(o);
                                }
                                else {
                                    //we must update neighbor information at the edges
                                    EdgeVisitor::removeEdges(&o, MapEdges);
                                }
                            }
                        }
                    }
                    //To mantain congruency in the map, we must erase all
                    //Quadrants (index) that have been split due to balancing.
                    auto delquad = idx_pos_map.find(key);
                    idx_pos_map.erase(delquad);
                }
            }
            
            // don't forget to update list
            std::swap(candidates,new_candidates);
            
            if (!new_pts.empty()) {
                //add the new points to the vector
                points.reserve(points.size() + new_pts.size());
                points.insert(points.end(),new_pts.begin(),new_pts.end());
                
                //cout << "new points inserted\n";
                //cout << "processed size " << processed.size() << endl;
                //cout << "map size " << idx_pos_map.size() << endl;
            }
            
            //unsigned int pro_quads = 0;
            
            //clean non used Quads.
            for (auto used_quad: idx_pos_map) {
                clean_processed.push_back(processed[used_quad.second]);
                //pro_quads++;
            }

            processed.erase(processed.begin(),processed.end());
            
        }
        
        //If there are no more refinement regions, we must apply transition patterns
        //at this moment and then finish the process.

        //----------------------------------------------------------
        // apply transition patterns
        //----------------------------------------------------------
        
        /*auto end_refine_quad_time = chrono::high_resolution_clock::now();
        
        //TransitionPatternVisitor section
        TransitionPatternVisitor tpv;
        tpv.setMapEdges(MapEdges);
        if (localmax) {
            tpv.setMaxRefLevel(maxrl+1);
        }
        else {
            tpv.setMaxRefLevel(maxrl);
        }
        new_pts.clear();

        //Apply transition patterns to remaining Quads
        unsigned mixedn = 0;
        for (auto &tq: clean_processed) {
            unsigned int sen = tq.getSubElements().size();
            if (!tq.accept(&tpv)) {
                std::cerr << "Error at Mesher::generateQuadtreeMesh";
                std::cerr << " Transition Pattern not found\n";
            }
            if (tq.getSubElements().size()!=sen) {
                mixedn++;
            }
        }
        
        //if no points were added at this iteration, it is no longer
        //necessary to continue the refinement.
        if (!new_pts.empty()) {
            //add the new points to the vector
            points.reserve(points.size() + new_pts.size());
            points.insert(points.end(),new_pts.begin(),new_pts.end());
        }*/
        
        //insert will reserve space as well
        Quadrants.insert(Quadrants.end(),make_move_iterator(candidates.begin()),make_move_iterator(candidates.end()));
        // better to erase as let in a indeterminate state by move
        candidates.erase(candidates.begin(),candidates.end());
        Quadrants.insert(Quadrants.end(),make_move_iterator(clean_processed.begin()),make_move_iterator(clean_processed.end()));
        clean_processed.erase(clean_processed.begin(),clean_processed.end());
        
#if (VTKOUT==true)            //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> transition_octree=make_shared<FEMesh>();
            saveOutputMesh(transition_octree,points,Quadrants);
            string tmp_name = name + "_transition";
            Services::WriteVTK(tmp_name,transition_octree);
        }
#endif
        
        if (localmax) {
            generateQuadtreeMesh(rl,input,all_reg,name,minrl,maxrl+1,debugging,future_q_idx);
        }
        else {
            generateQuadtreeMesh(rl,input,all_reg,name,minrl,maxrl,debugging,future_q_idx);
        }
        
        //If there are more refinement regions, continue with the process
        //were the quads positions will change in the final vector due to
        //map indexing that allows to optimize the research of neighbors.
        /*if (all_reg.size()>1) {
            
            auto end_time = chrono::high_resolution_clock::now();
            cout << "       * List refinement in "
            << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_refine_quad_time).count();
            cout << " ms"<< endl;
            
            if (listref) {
                
                //CL Debbuging
                {
                    //save pure octree mesh
                    std::shared_ptr<FEMesh> refined_octree=make_shared<FEMesh>();
                    saveOutputMesh(refined_octree,points,Quadrants);
                    string tmp_name = name + "_listRefinement";
                    Services::WriteVTK(tmp_name,refined_octree);
                }
            }
            //Continue with the rest of the refinement and apply transition patterns
            generateQuadtreeMesh(rl,input,all_reg,name,minrl,maxrl,debugging);
            return;
        }
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "       * Transition Patterns in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-end_refine_quad_time).count();
        cout << " ms"<< endl;
        cout << "    * generateQuadtreeMesh in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_refine_quad_time).count();
        cout << " ms"<< endl;*/
        
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    void Mesher::generateQuadtreeMesh(const unsigned short &rl, Polyline &input,
                                      const list<RefinementRegion *> &all_reg,
                                      const string &name, const unsigned short &minrl,
                                      const unsigned short &givenmaxrl,
                                      const bool &debugging,
                                      unsigned int new_q_idx) {
        
        
        auto start_time = chrono::high_resolution_clock::now();
        
        //The list of candidate quads to refine and the tmp version of
        //adding those how are still candidates for the next iteration.
        list<Quadrant> candidates, new_candidates, refine_tmp, clean_processed;
        
        //The Quads that don't need further refinement.
        vector<Quadrant> processed;
        
        //A map that connect the index of processed Quads with their position in the
        //vector of processed Quads. 
        map<unsigned int, unsigned int> idx_pos_map;
        
        //list of the points added at this refinement iteration:
        list<Point3D> new_pts;
        
        //A set containing the index of Quads to be refined to
        //maintain balancing
        list<pair<unsigned int,unsigned int> > toBalance;
        
        //initialising the vector and the map
        candidates.assign(make_move_iterator(Quadrants.begin()),
                          make_move_iterator(Quadrants.end()));
        
        //Erase previous quadrants to save memory
        Quadrants.clear();
        
        //an iterator for the regions
        list<RefinementRegion *>::const_iterator reg_iter;
        
        //create visitors and give them variables
        SplitVisitor sv;
        sv.setPoints(points);
        sv.setMapEdges(MapEdges);
        sv.setNewPts(new_pts);
        sv.setProcessedQuadVector(processed);
        sv.setMapProcessed(idx_pos_map);
        sv.setToBalanceList(toBalance);
        
        auto start_refine_quad_time = chrono::high_resolution_clock::now();
        
        //------------------------------------------------------------
        //refine each Quadrant until the Boundary is correctly handled
        //------------------------------------------------------------
        
        auto start_refine_rl_time = chrono::high_resolution_clock::now();
        
        //        list<RefinementRegion *>::const_iterator reg_iter=all_reg.begin();;
        unsigned int i=0; //current quad level
        /*do { // until no new quads are created
            new_pts.clear();
            
            //split the Quadrants as needed
            while (!tmp_Quadrants.empty()) {
                iter=tmp_Quadrants.begin();
                
                bool to_refine = false;
                
                iter->computeMaxDistance(points); //TODO, avoid recompute if already checked
                if ((*all_reg.begin())->intersectsQuadrant(points,*iter)) {
                    to_refine = true;
                }
                
                //now if refinement is not needed, we add the Quadrant as it was.
                if (!to_refine) {
                    new_Quadrants.push_back(*iter);
                    tmp_Quadrants.pop_front();
                    continue;
                }
                else {
                    list<unsigned int> &inter_edges = iter->getIntersectedEdges();
                    
                    vector<vector<Point3D> > clipping_coords;
                    sv.setClipping(clipping_coords);
                    
                    vector<vector<unsigned int> > split_elements;
                    sv.setNewEles(split_elements);
                    
                    sv.setStartIndex(new_Quadrants.size());
                    
                    //iter->split(points,new_pts,QuadEdges,split_elements,clipping_coords);
                    //cout << "Accept" << endl;
                    iter->accept(&sv);
                    
                    if (inter_edges.empty()) { //inner quad
                        for (unsigned int j=0; j<split_elements.size(); j++) {
                            Quadrant o (split_elements[j],i+1);
                            new_Quadrants.push_back(o);
                        }
                    }
                    else {
                        for (unsigned int j=0; j<split_elements.size(); j++) {
                            Quadrant o (split_elements[j],i+1);
                            //the new points are inserted in bash at the end of this
                            //iteration. For this reason, the coordinates must be passed
                            //"manually" at this point (clipping_coords).
                            
                            //select_faces = true
                            IntersectionsVisitor iv(true);
                            //if (o.checkIntersections(input,inter_edges,clipping_coords[j]))
                            iv.setPolyline(input);
                            iv.setEdges(inter_edges);
                            iv.setCoords(clipping_coords[j]);
                            
                            if (o.accept(&iv)) {
                                new_Quadrants.push_back(o);
                            }
                            else {
                                //The element doesn't intersect any input face.
                                //It must be checked if it's inside or outside.
                                //Only in the first case add it to new_Quadrants.
                                //Test this with parent Quadrant faces only.
                                
                                //Comment the following lines of this 'else' if
                                //only intersecting Quadrants are meant to be
                                //displayed.
                                
                                //note: inter_edges is quite enough to check if
                                //element is inside input, no Quadrant needed,
                                //so i moved the method to mesher  --setriva
                                
                                if (isItIn(input,inter_edges,clipping_coords[j])) {
                                    new_Quadrants.push_back(o);
                                }
                            }
                        }
                    }
                }
                // remove yet processed Quad
                tmp_Quadrants.pop_front();
                
            } // while
            
            // don't forget to update list
            std::swap(tmp_Quadrants,new_Quadrants);
            
            //if no points were added at this iteration, it is no longer
            //necessary to continue the refinement.
            if (new_pts.empty()) {
                cout << "warning at Mesher::generateQuadtreeMesh no new points!!!\n";
                break;
            }
            //add the new points to the vector
            points.reserve(points.size() + new_pts.size());
            points.insert(points.end(),new_pts.begin(),new_pts.end());
            
            
            ++i;
        } while (!new_pts.empty());
        
        #if (VTKOUT==true) //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> bound_octree=make_shared<FEMesh>();
            saveOutputMesh(bound_octree,points,new_candidates);
            string tmp_name = name + "_bound";
            Services::WriteVTK(tmp_name,bound_octree);
        }
        #endif */

        unsigned short max_rl;
        if (givenmaxrl==0) {
            max_rl = i;
        }
        else {
            max_rl = givenmaxrl;
        }
        
        //cout << "rl : [" << minrl << "," << max_rl << "]\n";
        
        
        auto end_refine_rl_time = chrono::high_resolution_clock::now();
        cout << "         * boundary max " << i << " in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_refine_rl_time-start_refine_rl_time).count();
        cout << " ms"<< endl;
        
        //----------------------------------------------------------
        //refine each Quadrant until the Refinement Level is reached
        //----------------------------------------------------------
        
        //when producing a new mesh, minrl will be always 0. But as this
        //method can also be called from "refineMesh" in this last case
        //the starting refinement level will not be 0, but the min rl
        //among the quadrants in the starting mesh.
        for (unsigned short i=minrl; i<rl; i++) {
            auto start_refine_rl_time = chrono::high_resolution_clock::now();
            
            //the new_pts is a list that holds the coordinates of
            //new points inserted at this iteration. At the end of
            //this bucle, they are inserted in the point vector
            new_pts.clear();
            
            //First we need to build the map with actual index
            //of quads in the processed quad vector. The quads
            //needing refinement will go into another list.
            while (!candidates.empty()) {
                Quadrant quad = *(candidates.begin());
                candidates.pop_front();
                
                bool to_refine = false;
                
                for (reg_iter=all_reg.begin(),reg_iter++; reg_iter!=all_reg.end(); ++reg_iter) {
                    
                    unsigned short region_rl = (*reg_iter)->getRefinementLevel();
                    if (region_rl<i) {
                        continue;
                    }
                    
                    //If the Quadrant has a greater RL than the region needs, continue
                    if (region_rl<= quad.getRefinementLevel()) {
                        continue;
                    }
                    
                    //in the case of a Region Refinement (a surface) this will change
                    //the state of inregion = true.
                    if ((*reg_iter)->intersectsQuadrant(points,quad)) {
                        to_refine = true;
                    }
                }
                
                //now if refinement is not needed, we add the Quadrant as it was.
                if (!to_refine) {
                    idx_pos_map[quad.getIndex()] = processed.size();
                    processed.push_back(quad);
                    continue;
                }
                else {
                    refine_tmp.push_back(quad);
                }
            }
            
            //now we can start to refine those needing it.
            for (auto quad: refine_tmp) {
                list<unsigned int> &inter_edges = quad.getIntersectedEdges();
                unsigned short qrl = quad.getRefinementLevel();
                
                vector<vector<Point3D> > clipping_coords;
                sv.setClipping(clipping_coords);
                
                vector<vector<unsigned int> > split_elements;
                sv.setNewEles(split_elements);
                sv.setStartIndex(new_q_idx);
                
                quad.accept(&sv);
                
                //If quad has not intersected faces directly add sons to new
                //candidates
                if (inter_edges.empty()) {
                    for (unsigned int j=0; j<split_elements.size(); j++) {
                        Quadrant o (split_elements[j], qrl+1, new_q_idx++);
                        //o.setInRegionState(reg_state);
                        new_candidates.push_back(o);
                    }
                }
                else {
                    for (unsigned int j=0; j<split_elements.size(); j++) {
                        Quadrant o (split_elements[j],qrl+1,new_q_idx++);
                        //the new points are inserted in bash at the end of this
                        //iteration. For this reason, the coordinates must be passed
                        //"manually" at this point (clipping_coords).
                        
                        //select_faces = true
                        IntersectionsVisitor iv(true);
                        //if (o.checkIntersections(input,inter_edges,clipping_coords[j]))
                        iv.setPolyline(input);
                        iv.setEdges(inter_edges);
                        iv.setCoords(clipping_coords[j]);
                        
                        if (o.accept(&iv)) {
                            //o.setInRegionState(reg_state);
                            new_candidates.push_back(o);
                        }
                        else {
                            //The element doesn't intersect any input edge.
                            //It must be checked if it's inside or outside.
                            //Only in the first case add it to new_Quadrants.
                            //Test this with parent Quadrant faces only.
                            if (isItIn(input,inter_edges,clipping_coords[j])) {
                                //o.setInRegionState(reg_state);
                                new_candidates.push_back(o);
                            }
                            else {
                                //we must update neighbor information at the edges
                                EdgeVisitor::removeEdges(&o, MapEdges);
                            }
                        }
                    }
                }
            }
            
            auto end_refine_rl_time = chrono::high_resolution_clock::now();
            cout << "         * level " << i << " in "
            << std::chrono::duration_cast<chrono::milliseconds>(end_refine_rl_time-start_refine_rl_time).count();
            cout << " ms";
            
            //Erase the list to refine
            refine_tmp.erase(refine_tmp.begin(),refine_tmp.end());
        
            //Refine non balanced Quads
            while (!toBalance.empty()) {
                
                list<pair<unsigned int, unsigned int> > tmp_toBalance;
                std::swap(toBalance,tmp_toBalance);
                tmp_toBalance.sort();
                tmp_toBalance.unique();
                
                while (!tmp_toBalance.empty()) {
                    unsigned int key = tmp_toBalance.begin()->first;
                    unsigned int val = tmp_toBalance.begin()->second;
                    tmp_toBalance.pop_front();

                    auto delquad = idx_pos_map.find(key);
                    if (delquad == idx_pos_map.end()) {
                        //Note: let us say that quad N is in the tmp_to_balance list at current
                        //iteration. Before arriving to it, a neighbor of N is refined to such
                        //level that adds N to list to_balance (for the next iteration). However
                        //as N is in the current tmp_to_balance it will be removed from the
                        //map and in next iteration it will be attempt to be split and removed
                        //once more in the mesh. Therefore, this if avoids this case.
                        continue;
                    }
                    
                    Quadrant quad = processed[val];
                    list<unsigned int> &inter_edges = quad.getIntersectedEdges();
                    unsigned short qrl = quad.getRefinementLevel();
                    
                    vector<vector<Point3D> > clipping_coords;
                    sv.setClipping(clipping_coords);
                    
                    vector<vector<unsigned int> > split_elements;
                    sv.setNewEles(split_elements);
                    sv.setStartIndex(new_q_idx);
                    
                    quad.accept(&sv);
                    
                    if (inter_edges.empty()) {
                        for (unsigned int j=0; j<split_elements.size(); j++) {
                            Quadrant o (split_elements[j], qrl+1, new_q_idx++);
                            idx_pos_map[o.getIndex()] = processed.size();
                            processed.push_back(o);
                        }
                    }
                    else {
                        for (unsigned int j=0; j<split_elements.size(); j++) {
                            Quadrant o (split_elements[j],qrl+1,new_q_idx++);
                            //the new points are inserted in bash at the end of this
                            //iteration. For this reason, the coordinates must be passed
                            //"manually" at this point (clipping_coords).
                            
                            //select_faces = true
                            IntersectionsVisitor iv(true);
                            //if (o.checkIntersections(input,inter_edges,clipping_coords[j]))
                            iv.setPolyline(input);
                            iv.setEdges(inter_edges);
                            iv.setCoords(clipping_coords[j]);
                            
                            if (o.accept(&iv)) {
                                idx_pos_map[o.getIndex()] = processed.size();
                                processed.push_back(o);
                                
                            }
                            else {
                                //The element doesn't intersect any input edge.
                                //It must be checked if it's inside or outside.
                                //Only in the first case add it to new_Quadrants.
                                //Test this with parent Quadrant faces only.
                                if (isItIn(input,inter_edges,clipping_coords[j])) {
                                    idx_pos_map[o.getIndex()] = processed.size();
                                    processed.push_back(o);
                                }
                                else {
                                    //we must update neighbor information at the edges
                                    EdgeVisitor::removeEdges(&o, MapEdges);
                                }
                            }
                        }
                    }
                    //To mantain congruency in the map, we must erase all
                    //Quadrants (index) that have been split due to balancing.
                    if (delquad != idx_pos_map.end()) {
                        idx_pos_map.erase(delquad);
                    }
                }
            }
            
            // don't forget to update list
            std::swap(candidates,new_candidates);
            
            //if no points were added at this iteration, it is no longer
            //necessary to continue the refinement.
            if (new_pts.empty()) {
                cerr << "warning at Mesher::generateQuadtreeMesh no new points!!!\n";
                auto end_balanced_time = chrono::high_resolution_clock::now();
                cout << " - Balanced in "
                << std::chrono::duration_cast<chrono::milliseconds>(end_balanced_time-end_refine_rl_time).count();
                cout << " ms"<< endl;
                break;
            }
            
            //add the new points to the vector
            points.reserve(points.size() + new_pts.size());
            points.insert(points.end(),new_pts.begin(),new_pts.end());
            
            auto end_balanced_time = chrono::high_resolution_clock::now();
            cout << " - Balanced in "
            << std::chrono::duration_cast<chrono::milliseconds>(end_balanced_time-end_refine_rl_time).count();
            cout << " ms"<< endl;
        }
        
        //clean non used Quads.
        for (auto used_quad: idx_pos_map) {
            clean_processed.push_back(processed[used_quad.second]);
        }
        
        //save space erasing processed quads.
        processed.erase(processed.begin(),processed.end());
        
#if (VTKOUT==true)        //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> refined_octree=make_shared<FEMesh>();
            if (debugging) {
                saveOutputMesh(refined_octree,points,clean_processed,true);
            }
            else {
                saveOutputMesh(refined_octree,points,clean_processed);
            }
            string tmp_name = name + "_refined";
            Services::WriteVTK(tmp_name,refined_octree);
        }
#endif

        auto end_refine_quad_time = chrono::high_resolution_clock::now();
        cout << "       * Refine Quad in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_refine_quad_time-start_refine_quad_time).count();
        cout << " ms"<< endl;

        
        //----------------------------------------------------------
        // apply transition patterns
        //----------------------------------------------------------
        
        //TransitionPatternVisitor section
        TransitionPatternVisitor tpv;
        tpv.setMapEdges(MapEdges);
        tpv.setMaxRefLevel(max_rl);
        new_pts.clear();
        
        //Apply transition patterns to remaining Quads
        unsigned mixedn = 0;
        for (auto &tq: clean_processed) {
            unsigned int sen = tq.getSubElements().size();
            if (!tq.accept(&tpv)) {
                std::cerr << "Error at Mesher::generateQuadtreeMesh";
                std::cerr << " Transition Pattern not found\n";
            }
            if (tq.getSubElements().size()!=sen) {
                mixedn++;
            }
        }
        
        //if no points were added at this iteration, it is no longer
        //necessary to continue the refinement.
        if (!new_pts.empty()) {
            //add the new points to the vector
            points.reserve(points.size() + new_pts.size());
            points.insert(points.end(),new_pts.begin(),new_pts.end());
        }
        
        //insert will reserve space as well
        Quadrants.insert(Quadrants.end(),make_move_iterator(candidates.begin()),make_move_iterator(candidates.end()));
        // better to erase as let in a indeterminate state by move
        candidates.erase(candidates.begin(),candidates.end());
        Quadrants.insert(Quadrants.end(),make_move_iterator(clean_processed.begin()),make_move_iterator(clean_processed.end()));
        clean_processed.erase(clean_processed.begin(),clean_processed.end());
        
#if (VTKOUT==true) //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> transition_octree=make_shared<FEMesh>();
            saveOutputMesh(transition_octree,points,Quadrants);
            string tmp_name = name + "_transition";
            Services::WriteVTK(tmp_name,transition_octree);
        }
#endif

        auto end_time = chrono::high_resolution_clock::now();
        cout << "       * Transition Patterns in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-end_refine_quad_time).count();
        cout << " ms"<< endl;
        cout << "    * generateQuadtreeMesh in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    bool Mesher::isItIn(const Polyline &mesh, const list<unsigned int> &faces, const vector<Point3D> &coords) const {
        //this method is meant to be used by Quadrants that don't
        //intersect input domains. If they are inside of at least
        //one input mesh, then they must remain in the output mesh.
        
        bool first = mesh.pointIsInMesh(coords[0],faces);
        bool second = mesh.pointIsInMesh(coords[1],faces);
        if (first==second) {
            return first;
        }
        
        //cout << "one inconsistency detected -> hard test\n";
        //return mesh.pointIsInMesh(coords[0],faces);
        return mesh.pointIsInMesh(coords[0]);
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    unsigned int Mesher::saveOutputMesh(const shared_ptr<FEMesh> &mesh, bool decoration){
        auto start_time = chrono::high_resolution_clock::now();
        
        vector<vector<unsigned int> > out_els;
        vector<unsigned short > out_els_ref_level, out_els_surf, out_els_deb;
        vector<double > out_els_min_angle, out_els_max_angle;
        array <unsigned int,180> out_els_angle_tri_histogram = {0},
                out_els_angle_quad_histogram = {0};

        //even if we don't know the quantity of elements, it will be at least the
        //number of quadrants, so we reserve for the vectors of VTK decoration.
        out_els_min_angle.reserve(Quadrants.size());
        out_els_max_angle.reserve(Quadrants.size());
        out_els_ref_level.reserve(Quadrants.size());
        out_els_surf.reserve(Quadrants.size());
        
        //new_idxs will hold the index of used nodes in the outside vector for points.
        //If the a node is not used by any element, its index will be 0 in this vector,
        //therefore the actual index is shiffted in 1. In other words, node 0 is node 1,
        //and node n is node n+1.
        vector<unsigned int> new_idxs (points.size(),0);
        unsigned int out_node_count = 0;
        vector<Point3D> out_pts;
 
        //recompute node indexes and update elements with them.
        for (unsigned int i=0; i<Quadrants.size(); i++) {
            const vector<vector<unsigned int> > &sub_els= Quadrants[i].getSubElements();
            for (unsigned int j=0; j<sub_els.size(); j++) {
                
                vector<unsigned int> sub_ele_new_idxs = sub_els[j];
                for (unsigned int k=0; k<sub_ele_new_idxs.size();k++) {
                    
                    unsigned int p_idx = sub_ele_new_idxs[k];
                    
                    if (new_idxs[p_idx]==0) {
                        sub_ele_new_idxs[k] = out_node_count++;
                        new_idxs[p_idx]=out_node_count;
                        out_pts.push_back(points[p_idx].getPoint());
                    }
                    else {
                        sub_ele_new_idxs[k] = new_idxs[p_idx]-1;
                    }
                }
                if (decoration) {
                    //surface regarding quad
                    if (Quadrants[i].isSurface()) {
                        out_els_surf.push_back(1);
                    }
                    else {
                        out_els_surf.push_back(0);
                    }
                    
                    if (Quadrants[i].isDebugging()) {
                        out_els_deb.push_back(1);
                    }
                    else {
                        out_els_deb.push_back(0);
                    }
                    
                    //refinment level herited from quad
                    out_els_ref_level.push_back(Quadrants[i].getRefinementLevel());


                    //compute minAngle and maxAngle (only for elements on the surface)
                    unsigned int np=sub_ele_new_idxs.size(); //nb points of the element
                    double minAngle=std::numeric_limits<double>::max();
                    double maxAngle=std::numeric_limits<double>::min();
                    if (Quadrants[i].isSurface()) {
                        for (unsigned int k=0; k<np; ++k) {

                            const Point3D &P0 = out_pts[sub_ele_new_idxs[(k-1+np)%np]];
                            const Point3D &P1 = out_pts[sub_ele_new_idxs[k]];
                            const Point3D &P2 = out_pts[sub_ele_new_idxs[(k+1)%np]];

                            double angle=P1.angle3Points(P0,P2);
                            minAngle=std::min(minAngle, angle);
                            maxAngle=std::max(maxAngle, angle);
                            if (np==3)
                                out_els_angle_tri_histogram[ min(179,((int)( round(angle)/*/10.*/)) %180) ]++;
                            else if (angle<89.5 || angle>91.5)
                                out_els_angle_quad_histogram[ min(179,((int)(round(angle)/*/10.*/)) %180) ]++;
                        }
                    }
                    out_els_min_angle.push_back(minAngle);
                    out_els_max_angle.push_back(maxAngle);
                }
                out_els.push_back(sub_ele_new_idxs);
            }
        }
        
        mesh->setPoints(out_pts);
        mesh->setElements(out_els);
        mesh->setRefLevels(out_els_ref_level);
        mesh->setMinAngles(out_els_min_angle);
        mesh->setMaxAngles(out_els_max_angle);
        mesh->setAnglesTriHistogram(out_els_angle_tri_histogram);
        mesh->setAnglesQuadHistogram(out_els_angle_quad_histogram);
        mesh->setSurfState(out_els_surf);
        mesh->setDebugging(out_els_deb);
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * SaveOutputMesh in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
        
        return out_els.size();
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    unsigned int Mesher::saveOutputMesh(const shared_ptr<FEMesh> &mesh,
                                        vector<MeshPoint> &tmp_points,
                                        list<Quadrant> &tmp_Quadrants,
                                        const bool &debugging,
                                        const list<Point3D> &extra_pts) {
        
        auto start_time = chrono::high_resolution_clock::now();
        
        vector<Point3D> out_pts;
        list<vector<unsigned int> > tmp_elements;
        vector<vector<unsigned int> > out_els;
        vector<unsigned short> deb_els;
        
        if (debugging) {
            deb_els.reserve(tmp_Quadrants.size());
        }
        
        unsigned int n = tmp_points.size();
        out_pts.reserve(n + extra_pts.size());
        for (unsigned int i=0; i<n; i++) {
            out_pts.push_back(points[i].getPoint());
        }
        
        //copy tmp points if available
        for (auto p: extra_pts) {
            out_pts.push_back(p);
        }
        
        OneIrregularVisitor oiv;
        if (debugging) {
            oiv.setEdges(MapEdges);
        }
        
        list<Quadrant>::iterator o_iter;
        
        bool errors = false;
        if (debugging) {
            for (o_iter=tmp_Quadrants.begin(); o_iter!=tmp_Quadrants.end(); ++o_iter) {
                
                vector<vector<unsigned int> > sub_els= o_iter->getSubElements();
                for (unsigned int j=0; j<sub_els.size(); j++) {
                    tmp_elements.push_back(sub_els[j]);
                }
                //Check balancing
                bool result = o_iter->accept(&oiv);
                if (result) {
                    deb_els.push_back(0);
                }
                else {
                    errors = true;
                    deb_els.push_back(1);
                }
            }
            mesh->setDebugging(deb_els);
            if (errors) {
                cerr << "->Warning mesh not balanced at Mesher::saveOutputMesh";
                cerr << " with debugging option on.\n";
            }
            else {
                cout << "    * No balancing problems detected\n";
            }
        }
        else {
            for (o_iter=tmp_Quadrants.begin(); o_iter!=tmp_Quadrants.end(); ++o_iter) {
                
                vector<vector<unsigned int> > sub_els= o_iter->getSubElements();
                for (unsigned int j=0; j<sub_els.size(); j++) {
                    tmp_elements.push_back(sub_els[j]);
                }
            }
        }
        
        out_els.reserve(tmp_elements.size());
        list<vector<unsigned int> >::const_iterator e_iter;
        
        for (e_iter=tmp_elements.begin(); e_iter!=tmp_elements.end(); ++e_iter) {
            out_els.push_back(*e_iter);
        }
        
        mesh->setPoints(out_pts);
        mesh->setElements(out_els);
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * SaveOutputMesh in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
        
        return out_els.size();
    }
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    unsigned int Mesher::saveOutputMesh(const shared_ptr<FEMesh> &mesh,
                                        vector<MeshPoint> &tmp_points,
                                        vector<Quadrant> &tmp_Quadrants,
                                        const bool &debugging,
                                        const list<Point3D> &extra_pts) {
        
        auto start_time = chrono::high_resolution_clock::now();
        
        vector<Point3D> out_pts;
        list<vector<unsigned int> > tmp_elements;
        vector<vector<unsigned int> > out_els;
        vector<unsigned short> deb_els;
        
        if (debugging) {
            deb_els.reserve(tmp_Quadrants.size());
        }
        
        unsigned int n = tmp_points.size();
        out_pts.reserve(n+extra_pts.size());
        for (unsigned int i=0; i<n; i++) {
            out_pts.push_back(points[i].getPoint());
        }
        
        OneIrregularVisitor oiv;
        if (debugging) {
            oiv.setEdges(MapEdges);
        }
        
        
        //copy tmp points if available
        for (auto p: extra_pts) {
            out_pts.push_back(p);
        }
        
        unsigned int frl = tmp_Quadrants[0].getRefinementLevel();
        
        bool errors = false;
        
        if (debugging) {
            for (unsigned int i=0; i<tmp_Quadrants.size(); i++) {
                vector<vector<unsigned int> > sub_els= tmp_Quadrants[i].getSubElements();
                for (unsigned int j=0; j<sub_els.size(); j++) {
                    tmp_elements.push_back(sub_els[j]);
                }
                //Check balancing
                bool result = tmp_Quadrants[i].accept(&oiv);
                if (result) {
                    deb_els.push_back(0);
                }
                else {
                    errors = true;
                    deb_els.push_back(1);
                }
            }
            mesh->setDebugging(deb_els);
            if (errors) {
                cerr << "->Warning mesh not balanced at Mesher::saveOutputMesh";
                cerr << " with debugging option on.\n";
            }
            else {
                cout << "    * No balancing problems detected\n";
            }
        }
        else {
            for (unsigned int i=0; i<tmp_Quadrants.size(); i++) {
                vector<vector<unsigned int> > sub_els= tmp_Quadrants[i].getSubElements();
                for (unsigned int j=0; j<sub_els.size(); j++) {
                    tmp_elements.push_back(sub_els[j]);
                }
            }
        }
        
        out_els.reserve(tmp_elements.size());
        list<vector<unsigned int> >::const_iterator e_iter;
        
        for (e_iter=tmp_elements.begin(); e_iter!=tmp_elements.end(); ++e_iter) {
            out_els.push_back(*e_iter);
        }
        
        mesh->setPoints(out_pts);
        mesh->setElements(out_els);
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * SaveOutputMesh in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
        
        return out_els.size();
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void Mesher::detectFeatureQuadrants(Polyline &input) {
        
        auto start_time = chrono::high_resolution_clock::now();
        
        unsigned int featCount = 0;
        for (unsigned int i=0; i<Quadrants.size(); i++) {
            if (input.getNbFeatures(Quadrants[i],points)>0) {
                //Quadrants[i].setFeature(); now set by Polyline::getNbFeatures()
                featCount++;
            }
        }
        //cout << "number of quadrants with features: " << featCount << "\n";
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * detectFeatureQuadrants in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
        
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void Mesher::linkElementsToNodes(){
        
        auto start_time = chrono::high_resolution_clock::now();
        
        //clear previous information
        for (unsigned int i=0; i<points.size(); i++) {
            points[i].clearElements();
        }
        
        //link element info to nodes
        for (unsigned int i=0; i<Quadrants.size(); i++) {
            
            const vector <unsigned int> &q_indpts = Quadrants[i].getSubPointIndex();
            
            for (unsigned int j=0; j<q_indpts.size(); j++) {
                points.at(q_indpts[j]).addElement(i);
            }
        }
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * linkElementsToNodes in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    //WARNING: linkElementsToNodes() must be called before
    void Mesher::computeNodeMaxDist() {
        for (auto &q:Quadrants) {
            q.computeMaxDistance(points);
        }
        for (auto& p:points) {
            for (const auto &e:p.getElements()) {
                p.updateMaxDistance(Quadrants[e].getMaxDistance());
            }
        }
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    void Mesher::detectInsideNodes(Polyline &input){
        auto start_time = chrono::high_resolution_clock::now();
        
        for (unsigned int i=0; i<points.size(); i++) {
            if (points[i].wasOutsideChecked()) {
                continue;
            }
            
            list<unsigned int> p_eles = points[i].getElements(), p_edges;
            points[i].outsideChecked();
            if (p_eles.empty()) {
                continue;
            }
            list<unsigned int>::const_iterator iter;
            for (iter=p_eles.begin(); iter!=p_eles.end(); ++iter) {
                const list<unsigned int> &qedges= Quadrants[*iter].getIntersectedEdges();
                //                list<unsigned int>::const_iterator qe_iter;
                //                if (qedges.empty()) {
                //                    continue;
                //                }
                // append qedges to p_edges
                p_edges.insert(p_edges.end(),qedges.begin(),qedges.end());
            }
            
            p_edges.sort();
            p_edges.unique();
            
            // p_edges: edges intersected
            if (p_edges.empty() || input.pointIsInMesh(points[i].getPoint(),p_edges)) {
                points[i].setInside();
            }
        }
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * detectInsideNodes in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    void Mesher::removeOnSurface(Polyline &input){
        auto start_time = chrono::high_resolution_clock::now();
        
        list<Quadrant> newele,removed;
        RemoveSubElementsVisitor rsv;
        rsv.setPoints(points);
        //remove elements without an inside node.
        for (unsigned int i=0; i<Quadrants.size(); i++) {
            if (Quadrants[i].isInside()) {
                newele.push_back(Quadrants[i]);
                continue;
            }
            
            if (Quadrants[i].hasIntersectedFeatures()) {
                if (Quadrants[i].accept(&rsv)) {
                    //Quadrant with feature to be removed
                    //only if average node is outside the
                    //the edges it intersects.
                    Point3D avg;
                    for (auto quaNoIdx:Quadrants[i].getPointIndex()) {
                        avg+=points[quaNoIdx].getPoint();
                    }
                    avg/=Quadrants[i].getPointIndex().size();
                    if (input.pointIsInMesh(avg,Quadrants[i].getIntersectedEdges())) {
                        newele.push_back(Quadrants[i]);
                    }
                    else {
                        removed.push_back(Quadrants[i]);
                    }
                }
                else {
                    newele.push_back(Quadrants[i]);
                }
            }
            else { //FJA add a "else" here as some quadrants are inserted twice
                
                //if (Quadrants[i].removeOutsideSubElements(points)) {
                if (Quadrants[i].accept(&rsv)) {
                    removed.push_back(Quadrants[i]);
                }
                else {
                    newele.push_back(Quadrants[i]);
                }
            }
        }
        
        if (removed.empty()) {
            return;
        }
        
        //clear removed elements
        removed.clear();
        //now element std::list from Surface mesh can be cleared, as all remaining
        //elements are still in use and attached to newele std::list.
        Quadrants.clear();
        Quadrants.assign(make_move_iterator(newele.begin()),make_move_iterator(newele.end()));
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * RemoveOnSurface in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void Mesher::removeOnSurfaceSafe(Polyline &input){
        auto start_time = chrono::high_resolution_clock::now();
        
        list<Quadrant> newele,removed;
        
        //Keep all inside quads and all of th
        for (auto q:Quadrants) {
            if (q.isInside()) {
                newele.push_back(q);
                continue;
            }
            
            bool onein = false, feaProj = false;
            //Check inside nodes, including internal ones
            for (auto quaNoIdx:q.getSubPointIndex()) {
                if (points[quaNoIdx].isInside()) {
                    onein = true;
                    break;
                }
                if (points[quaNoIdx].isFeature()) {
                    feaProj = true;
                }
            }
            
            //if it has at least one node inside
            //we don't remove it.
            if (onein) {
                newele.push_back(q);
                continue;
            }
            
            //Compute average node of the Quadrant
            //only with corner nodes
            Point3D avg;
            for (auto quaNoIdx:q.getPointIndex()) {
                avg+=points[quaNoIdx].getPoint();
            }
            avg/=q.getPointIndex().size();
            //Sometimes the avg node is just over domain's boundary,
            //despite the fact that a (big) portion of it's still
            //inside the domain and no other Quad will cover this
            //section. For this reason the entire element is shrink
            //to 99% of it. If one node of this shrink element is still
            //inside the domain, the Octant remains. Recall that a
            //node projected onto the surface is treated as an outside
            //node.
            bool accepted = false;
            for (auto quaNoIdx:q.getPointIndex()) {
                //Point3D test = (avg + points[quaNoIdx].getPoint())/2;
                Point3D test = points[quaNoIdx].getPoint() - avg;
                test*=0.99;
                test+=avg;
                if (input.pointIsInMesh(test,q.getIntersectedEdges())) {
                    newele.push_back(q);
                    accepted = true;
                    break;
                }
            }
            if (!accepted) {
                removed.push_back(q);
            }
            
            //the most expensive case: regarding
            //the original intersected edges, does
            //this quad still have a node inside?
            //if yes we keep it. It means that all
            //of its nodes are outside or projected
            //onto the boundary so if we erase it we
            //will loose boundary representation.
            /*bool border = false;
             for (auto qp:q.getPointIndex()) {
             if (input.pointIsInMesh(points[qp].getPoint(),q.getIntersectedEdges())) {
             newele.push_back(q);
             border = true;
             break;
             }
             }
             if (!border) {
             removed.push_back(q);
             }*/
        }
        
        if (removed.empty()) {
            return;
        }
        
        //clear removed elements
        removed.clear();
        //now element std::list from Vomule mesh can be cleared, as all remaining
        //elements are still in use and attached to newele std::list.
        Quadrants.clear();
        Quadrants.assign(std::make_move_iterator(newele.begin()),std::make_move_iterator(newele.end()));
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * RemoveOnSurfaceSafe in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    void Mesher::applySurfacePatterns(Polyline &input){
        auto start_time = chrono::high_resolution_clock::now();
        
        //apply patterns to avoid flat, invalid and
        //poor quality elements.
        SurfaceTemplatesVisitor stv;
        stv.setPoints(points);
        stv.setPolyline(input);
        
        
        for (auto &q:Quadrants) {
            
            if (q.getPointIndex().size()!=4) {
                continue;
            }
            
            if (q.isSurface()) {
                if (!q.accept(&stv)) {
                    cout << "Error in Mesher::applySurfacePatterns: coultd't apply";
                    cout << " a surface pattern\n";
                    cout << q << "\n";
                    continue;
                }
            }
            
        }
        
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * ApplySurfacePatterns in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }
    
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    //shrink elements intersecting the envelope defined by all
    //input surfaces
    
    void Mesher::shrinkToBoundary(Polyline &input){
        auto start_time = chrono::high_resolution_clock::now();
        
        //Slow element removed (but works): from elements intersecting the
        //input domain, detect inner nodes. Project this nodes onto the
        //surface. If after all is done, if an element counts only with "on
        //surface" and "outside" nodes, remove it.
        list<unsigned int> out_nodes;
        list<Quadrant>::iterator oiter;
        
        
        //Manage Quadrants with Features first
        for (const auto &q:Quadrants) {
            if (q.isInside()) {
                continue;
            }
            
            //Save all outside nodes and manage them after
            //dealing with all features quadrants.
            //if a node was projected to a feature
            //it will be skipped during the rest of
            //node projection.
            for (auto pIdx:q.getSubPointIndex()) {
                if (points[pIdx].isOutside()) {
                    out_nodes.push_back(pIdx);
                }
            }
            
            /*
             All features were already managed.
             
             if (!q.hasFeature()) {
             continue;
             }
             
             list<unsigned int> fs = input.getFeatureProjection(q,points);
             
             if (fs.empty()) {
             cerr << "Error at Mesher::shrinkToBoundary";
             cerr << " Quadrant labeled with feature has none\n";
             continue;
             }
             
             const vector<unsigned int> &epts = q.getPointIndex();
             
             unsigned int fsNum = fs.size(), outNo = 0;
             for (auto pIdx:epts) {
             if (points[pIdx].isOutside() && !points[pIdx].wasProjected()) {
             outNo++;
             }
             }
             
             list<unsigned int>::const_iterator iter;
             
             for (iter=fs.begin(); iter!=fs.end(); ++iter) {
             
             if (outNo==0) {
             break;
             }
             
             double best = std::numeric_limits<double>::infinity();
             Point3D projected = input.getPoints()[*iter];
             unsigned int pos = 0;
             bool push = false;
             
             for (auto pIdx:epts) {
             if (points[pIdx].isOutside() && !points[pIdx].wasProjected()) {
             
             const Point3D &current = points[pIdx].getPoint();
             double dis = (current - projected).Norm();
             
             if(best>dis){
             best = dis;
             pos = pIdx;
             push = true;
             }
             }
             }
             
             if (push) {
             points[pos].setProjected();
             points[pos].setPoint(projected);
             //Feature projected flag will be used later to
             //apply surface patterns.
             points[pos].featureProjected();
             for (auto pe:points[pos].getElements()) {
             
             //this should be studied further.
             if (Quadrants.at(pe).intersectsSurface()) {
             Quadrants[pe].setSurface();
             }
             }
             outNo--;
             }
             }*/
        }
        
        //Manage non Feature Quadrants.
        out_nodes.sort();
        out_nodes.unique();
        
        for (auto p:out_nodes) {
            
            if (points.at(p).wasProjected()) {
                continue;
            }
            
            //get the faces of Quadrants sharing this node
            list<unsigned int> p_qInterEdges;
            
            for (auto pe:points[p].getElements()) { //elements containing p
                //append this to the list of edges
                p_qInterEdges.insert(p_qInterEdges.end(),
                                     Quadrants[pe].getIntersectedEdges().begin(),
                                     Quadrants[pe].getIntersectedEdges().end());
            }
            
            p_qInterEdges.sort();
            p_qInterEdges.unique();
            
            if (p_qInterEdges.empty()) {
                cout << "\nWarning at Mesher::shrinkToBoundary";
                cout << " no faces to project an outside node\n";
                cout << p << " n_els " << points.at(p).getElements().size() << ":";
                for (auto pe:points.at(p).getElements()) {
                    cout << " " << pe;
                }
                cout << "\n";
                continue;
            }
            
            const Point3D &current = points[p].getPoint();
            Point3D projected = input.getProjection(current,p_qInterEdges);
            
            //if ( current.distance(projected)< points[p].getMaxDistance() )
            {
                points[p].setPoint(projected);
                points[p].setProjected();
            }
        }
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * ShrinkToBoundary in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }
    
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    void Mesher::projectCloseToBoundaryNodes(Polyline &input){
        auto start_time = chrono::high_resolution_clock::now();
        
        //Slow element removed (but works): from elements intersecting the
        //input domain, detect inner nodes. Project this nodes onto the surface.
        list<unsigned int> in_nodes;
        
        for (auto &q:Quadrants) {
            
            //if (!q.isSurface()) {
            if (q.isInside()) {
                continue;
            }
            //    continue;
            //}
            
            //Save the intern nodes to a list to manage them later
            //and continue the process just for feature Quadrants first.
            const vector<unsigned int> subpointindex(q.getSubPointIndex());
            for (auto pIdx:subpointindex) {
                if (points[pIdx].isInside()) {
                    in_nodes.push_back(pIdx);
                }
            }
            if (!q.hasIntersectedFeatures()) {
                continue;
            }
            
            //Manage Quadrants with Features
            list<unsigned int> fs = input.getFeatureProjection(q,points);
            
            unsigned int fsNum = fs.size();
            //we use a list and interator to erase the projected node
            //from the list of features.
            list<unsigned int>::const_iterator itFs;
            for (itFs=fs.begin(); itFs!=fs.end(); ++itFs) {
                
                //                if (fsNum==0) {
                //                    break;
                //                }
                if ( input.featureHasProjectedPt(*itFs)) {
                    break;
                }
                
                const Point3D &featProjected = input.getPoints()[*itFs];
                
                double best = std::numeric_limits<double>::infinity();
                unsigned int candidate;
                
                bool push = false;
                for (auto pIdx:subpointindex) {
                    
                    if (points[pIdx].wasProjected()) {
                        /* no more needed now we store a map Feature-ProjectedPt !!!
                         // special case: when a Feature is right on an edge/node
                         // and as already been treated in the neighboor quadrant
                         //FJA: work only if ONLY ONE feature per quad
                         // otherwise, should store the Feature number in the MeshPoint
                         // but memory consuming...
                         if (points[pIdx].isFeature()) {
                         push=false;
                         break;
                         } */
                        continue;
                    }
                    
                    const Point3D &current = points[pIdx].getPoint();
                    double dis = (current - featProjected).Norm();
                    
                    //if(points[pIdx].getMaxDistance()>dis && best>dis){
                    if(best>dis) {
                        best = dis;
                        candidate = pIdx;
                        push = true;
                    }
                }
                
                if (push) {
                    points[candidate].setProjected();
                    points[candidate].setPoint(featProjected);
                    //Feature projected flag will be used later to
                    //apply surface patterns.
                    points[candidate].setFeature();
                    input.setFeatureIndexProjectedPt(*itFs,candidate);
                    
                    for (auto pe:points[candidate].getElements()) {
                        //all the quads associated to a projected node
                        //will be labeled as surface, even if they don't
                        //intersect a boundary. This is to enable surfacePatterns
                        //for each quadrant that must be treated (note that
                        //a surface quad is not the same as an inside quad).
                        Quadrants[pe].setSurface();
                        
                        //Also, if quadrant was completely inside, then we must
                        //set as intersected edges the edges of its neighbor.
                        //This information will be used to detect if sub-elements
                        //are still inside the domain and for that, the subset
                        //of intersected edges is necessary.
                        if (Quadrants[pe].getIntersectedEdges().empty()) {
                            Quadrants[pe].setIntersectedEdges(q.getIntersectedEdges());
                        }
                    }
                    
                    fsNum--;
                }
            }
        }
        
        in_nodes.sort();
        in_nodes.unique(); //need to be sorted to call this function
        
        //move (when possible) all inner points to surface
        for (auto p:in_nodes) {
            
            if (points[p].wasProjected()) {
                //projected due to a feature => already managed.
                continue;
            }
            
            //get the faces of Quadrants sharing this node
            list<unsigned int> p_qInterEdges;
            
            //elements containing p
            for (auto pe:points.at(p).getElements()) {
                //append this to the list of edges
                p_qInterEdges.insert(p_qInterEdges.end(),
                                     Quadrants[pe].getIntersectedEdges().begin(),
                                     Quadrants[pe].getIntersectedEdges().end());
            }
            
            p_qInterEdges.sort();
            p_qInterEdges.unique();
            
            const Point3D &current = points[p].getPoint();
            Point3D projected = input.getProjection(current,p_qInterEdges);
            double dis = (current - projected).Norm();
            
            if(dis<points[p].getMaxDistance()){
                //this node has been moved to boundary, thus every element
                //sharing this node must be set as a border element in order
                //to avoid topological problems.
                points[p].setProjected();
                points[p].setPoint(projected);
                for (auto pe:points[p].getElements()) {
                    
                    //intersectsSurface() returns true when the Quad previously
                    //had non empty list of intersected edges.
                    if (Quadrants[pe].intersectsSurface()) {
                        //setSurface() will change the flags of "inside" wich
                        Quadrants[pe].setSurface();
                    }
                }
            }
        }
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * ProjectCloseToBoundary in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
        
    }
}
